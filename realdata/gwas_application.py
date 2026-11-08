import argparse, gzip, hashlib, json, os, time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from bed_reader import open_bed
from scipy.stats import norm, chi2, t as tdist
from xcal.core import panel_moments, jackknife_moments, satterthwaite_df, ld_scores, estimate_c, sd_from_eaf

ap = argparse.ArgumentParser(description="Apply standard TWAS, adapted EUR variance control and XCAL-t to a real GWAS (BioBank Japan height)")
ap.add_argument("--sumstats", required=True, help="BBJ BOLT file: Variants CHR POS REF ALT ALT_freq MAF Rsq BETA SE P_INF P_BOLT (hg19)")
ap.add_argument("--geno", required=True, help="1000 Genomes PLINK prefix template with {} for super-population")
ap.add_argument("--models", required=True, help="directory with mesa_{AFA,CAU,HIS}_weights.csv.gz and mesa_*_extra.csv")
ap.add_argument("--target", default="EAS"); ap.add_argument("--n-gwas", type=int, default=159095)
ap.add_argument("--outdir", required=True); ap.add_argument("--threads", type=int, default=2)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
SETS = ["AFA", "CAU", "HIS"]
WIN = 1_000_000
t0 = time.time()
log = open(f"{a.outdir}/application_log.txt", "w", encoding="utf-8")


def say(*x):
    s = " ".join(str(i) for i in x)
    print(s, flush=True); log.write(s + "\n"); log.flush()


bim = pd.read_csv(a.geno.format(a.target) + ".bim", sep="\t", header=None, names=["chr", "rsid", "cm", "pos", "a1", "a2"])
bim["idx"] = np.arange(len(bim))
bim["key"] = bim.chr.astype(str) + ":" + bim.pos.astype(str)
keys = set(bim.key)
parts = []
for ch in pd.read_csv(a.sumstats, sep=r"\s+", usecols=["CHR", "POS", "REF", "ALT", "ALT_freq", "Rsq", "BETA", "SE"], chunksize=2_000_000,
                      dtype={"CHR": str}):
    ch["key"] = ch.CHR.astype(str) + ":" + ch.POS.astype(str)
    parts.append(ch[ch.key.isin(keys) & (ch.Rsq >= 0.8)])
ss = pd.concat(parts, ignore_index=True)
ss = ss.merge(bim[["key", "rsid", "a1", "a2", "idx"]], on="key")
same = (ss.ALT == ss.a1) & (ss.REF == ss.a2)
flip = (ss.ALT == ss.a2) & (ss.REF == ss.a1)
ss = ss[same | flip].copy()
ss["z"] = np.where(flip[same | flip], -1, 1) * ss.BETA / ss.SE
ss["eaf_a1"] = np.where(flip[same | flip], 1 - ss.ALT_freq, ss.ALT_freq)
ss = ss[(ss.eaf_a1 > 0.01) & (ss.eaf_a1 < 0.99) & np.isfinite(ss.z)]
ss = ss.drop_duplicates("rsid", keep=False)
say(f"matched SNPs: {len(ss):,} of {len(bim):,} panel SNPs ({time.time() - t0:.0f}s)")

with open_bed(a.geno.format(a.target) + ".bed", count_A1=True, num_threads=a.threads) as bed:
    n_target = bed.iid_count
ldf = f"{a.outdir}/ld_scores_{a.target}.npy"
if os.path.exists(ldf):
    ldl = np.load(ldf)
else:
    ldl = np.full(len(bim), np.nan)
    for chrom in range(1, 23):
        idx = bim.index[bim.chr == chrom].values
        with open_bed(a.geno.format(a.target) + ".bed", count_A1=True, num_threads=a.threads) as bed:
            X = bed.read(index=np.s_[:, idx], dtype="float32")
        X = np.where(X < 0, np.nan, X)
        X = np.where(np.isnan(X), np.nanmean(X, 0), X)
        ldl[idx] = ld_scores(X, bim.pos.values[idx], window_bp=WIN)
        say(f"LD scores chr{chrom} ({time.time() - t0:.0f}s)")
    np.save(ldf + ".tmp.npy", ldl); os.replace(ldf + ".tmp.npy", ldf)
ss["ld"] = ldl[ss.idx.values]
fitf = f"{a.outdir}/ldsc_fits.json"
if os.path.exists(fitf):
    fit_free, fit_fixed = json.load(open(fitf))
else:
    fit_free = estimate_c(ss.z.values, ss.ld.values, n_blocks=200)
    fit_fixed = estimate_c(ss.z.values, ss.ld.values, n_blocks=200, fix_intercept=1.0)
    json.dump([fit_free, fit_fixed], open(fitf, "w"))
M = int(np.isfinite(ss.ld).sum())
say("LDSC free intercept:", json.dumps(fit_free), "\nLDSC fixed intercept:", json.dumps(fit_fixed))
say(f"implied N*h2 = c*M = {fit_free['c'] * M:,.0f} (M = {M:,}); h2_SNP ~ {fit_free['c'] * M / a.n_gwas:.3f}")
loof = f"{a.outdir}/ldsc_loco.json"
if os.path.exists(loof):
    loco = {int(k): v for k, v in json.load(open(loof)).items()}
else:
    chr_ss = ss.CHR.astype(int).values
    loco = {}
    for chrom in range(1, 23):
        keep = chr_ss != chrom
        loco[chrom] = estimate_c(ss.z.values[keep], ss.ld.values[keep], n_blocks=200)
        say(f"LOCO LDSC without chr{chrom}: c = {loco[chrom]['c']:.4f}, a = {loco[chrom]['intercept']:.3f}")
    json.dump(loco, open(loof + ".tmp", "w")); os.replace(loof + ".tmp", loof)

zmap = pd.Series(ss.z.values, index=ss.idx.values)
sdmap = pd.Series(sd_from_eaf(ss.eaf_a1.values), index=ss.idx.values)
W = pd.concat([pd.read_csv(f"{a.models}/mesa_{s}_weights.csv.gz").assign(set=s) for s in SETS])
E = pd.concat([pd.read_csv(f"{a.models}/mesa_{s}_extra.csv").assign(set=s) for s in SETS])
keep = E[E["pred.perf.pval"] < 0.05][["set", "gene", "genename"]]
W = W.merge(keep, on=["set", "gene"]).merge(bim[["rsid", "chr", "pos", "a1", "a2", "idx"]], on="rsid")
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
s1 = (W.eff_allele == W.a1) & (W.ref_allele == W.a2) | (W.eff_allele.map(COMP) == W.a1) & (W.ref_allele.map(COMP) == W.a2)
s2 = (W.eff_allele == W.a2) & (W.ref_allele == W.a1) | (W.eff_allele.map(COMP) == W.a2) & (W.ref_allele.map(COMP) == W.a1)
amb = W.eff_allele.map(COMP) == W.ref_allele
W = W[(s1 | s2) & ~amb].copy()
W["w"] = np.where(s1[(s1 | s2) & ~amb], W.weight, -W.weight)
W = W[W.idx.isin(zmap.index)]
say(f"models with >=1 usable SNP: {W.groupby(['set', 'gene']).ngroups:,}")


def std(X):
    X = X.astype(np.float32)
    X[X < 0] = np.nan
    mu = np.nanmean(X, 0)
    bad = np.isnan(X)
    if bad.any():
        X[bad] = np.take(mu, np.nonzero(bad)[1])
    sd = X.std(0); ok = sd > 0
    X -= X.mean(0)
    X[:, ok] /= sd[ok]
    X[:, ~ok] = 0.0
    return X, sd.astype(np.float64), ok


CHUNK, PAD = 20_000_000, 3_000_000
rows = []
os.makedirs(f"{a.outdir}/per_chr", exist_ok=True)
for chrom in range(1, 23):
    pcf = f"{a.outdir}/per_chr/chr{chrom}.csv"
    if os.path.exists(pcf):
        if os.path.getsize(pcf) > 1:
            rows.extend(pd.read_csv(pcf).to_dict("records"))
        continue
    r0 = len(rows)
    fc = loco[chrom]
    c, aint, se_c, se_a = fc["c"], max(fc["intercept"], 1.0), fc["se_c"], fc["se_intercept"]
    bc_all = bim[bim.chr == chrom]
    Wc_all = W[W.chr == chrom]
    gmin = Wc_all.groupby(["set", "gene"]).pos.min()
    for cs in range(0, int(bc_all.pos.max()) + 1, CHUNK):
        genes_here = gmin[(gmin >= cs) & (gmin < cs + CHUNK)].index
        if len(genes_here) == 0:
            continue
        bc = bc_all[(bc_all.pos >= cs - PAD) & (bc_all.pos < cs + CHUNK + PAD)]
        with open_bed(a.geno.format(a.target) + ".bed", count_A1=True, num_threads=a.threads) as bed:
            ZT, sT, okT = std(bed.read(index=np.s_[:, bc.idx.values], dtype="int8"))
        with open_bed(a.geno.format("EUR") + ".bed", count_A1=True, num_threads=a.threads) as bed:
            ZE, sE, okE = std(bed.read(index=np.s_[:, bc.idx.values], dtype="int8"))
        pos = bc.pos.values; base = bc.idx.values[0]
        blocks = np.array_split(np.random.default_rng(chrom * 1000 + cs // CHUNK).permutation(ZT.shape[0]), 20)
        Wc = Wc_all.set_index(["set", "gene"]).loc[genes_here].reset_index()
        for (s, gene), g in Wc.groupby(["set", "gene"]):
            loc = g.idx.values - base
            lo = np.searchsorted(pos, pos[loc].min() - WIN); hi = np.searchsorted(pos, pos[loc].max() + WIN, side="right")
            sl = loc - lo
            w = np.zeros(hi - lo); w[sl] = g.w.values
            zw = np.zeros(hi - lo); zw[sl] = zmap[g.idx.values].values
            gsd = np.zeros(hi - lo); gsd[sl] = sdmap[g.idx.values].values
            vG, vT, vE = np.where(okT[lo:hi], w * gsd, 0.0), w * sT[lo:hi], w * sE[lo:hi]
            ZTw, ZEw = ZT[:, lo:hi].astype(np.float64), ZE[:, lo:hi].astype(np.float64)
            mG = panel_moments(ZTw, vG, okT[lo:hi]); mT = panel_moments(ZTw, vT, okT[lo:hi])
            mE = panel_moments(ZEw, vE, okE[lo:hi])
            if min(mG["L"], mT["L"], mE["L"]) <= 0:
                continue
            Lb, Qb = jackknife_moments(ZTw, vG, blocks=blocks)
            D = aint * mG["L"] + c * mG["Q"]
            df = float(satterthwaite_df(aint * Lb + c * Qb, D, extra_var=(mG["Q"] * se_c) ** 2 + (mG["L"] * se_a) ** 2))
            qE = mE["Q"] / mE["L"]
            zG, zT, zE = vG @ zw, vT @ zw, vE @ zw
            r = dict(set=s, gene=gene, genename=g.genename.iloc[0], chr=chrom, n_snps=len(g),
                     z_standard_target=zT / np.sqrt(mT["L"]), z_standard_eur=zE / np.sqrt(mE["L"]),
                     z_vc_eur_eurld=zE / np.sqrt(mE["L"] * (1 + c * qE)), z_vc_eur_targetld=zT / np.sqrt(mT["L"] * (1 + c * qE)),
                     z_xcal=zG / np.sqrt(D), df_xcalt=df, c_used=c, a_used=aint, q_target=mG["Q"] / mG["L"], q_eur=qE, floor=mG["floor_active"])
            r["p_xcalt"] = 2 * tdist.sf(abs(r["z_xcal"]), df)
            for k in ("z_standard_target", "z_standard_eur", "z_vc_eur_eurld", "z_vc_eur_targetld", "z_xcal"):
                r[k.replace("z_", "p_")] = 2 * norm.sf(abs(r[k]))
            rows.append(r)
        del ZT, ZE
    pd.DataFrame(rows[r0:]).to_csv(pcf + ".tmp", index=False); os.replace(pcf + ".tmp", pcf)
    say(f"chr{chrom}: {len(rows):,} models ({time.time() - t0:.0f}s)")
R = pd.DataFrame(rows)
R.to_csv(f"{a.outdir}/application_per_model.csv.gz", index=False)
n = len(R); bonf = 0.05 / n
methods = {"Standard TWAS, EAS LD": "p_standard_target", "Standard TWAS, EUR LD": "p_standard_eur",
           "Adapted variance control (EUR q), EUR LD": "p_vc_eur_eurld", "Adapted variance control (EUR q), EAS LD": "p_vc_eur_targetld",
           "XCAL (normal reference)": "p_xcal", "XCAL-t (GWAS SDs; recommended)": "p_xcalt"}
S = []
for lab, col in methods.items():
    p = R[col].values
    S.append(dict(method=lab, n_models=n, lambda_median=float(np.median(chi2.isf(p, 1)) / chi2.isf(0.5, 1)),
                  n_sig_bonferroni=int((p < bonf).sum()), n_sig_genomewide_2p5e6=int((p < 2.5e-6).sum()),
                  n_distinct_genes_bonf=int(R.loc[p < bonf, "genename"].nunique())))
S = pd.DataFrame(S)
S.to_csv(f"{a.outdir}/application_summary.csv", index=False)
ratio = (R.q_target / R.q_eur).describe(percentiles=[.05, .5, .95])
ratio.to_csv(f"{a.outdir}/application_q_ratio.csv")
json.dump(dict(fit_free=fit_free, fit_fixed=fit_fixed, loco={str(k): dict(c=v["c"], se_c=v["se_c"], intercept=v["intercept"], se_intercept=v["se_intercept"]) for k, v in loco.items()}, M=M, n_gwas=a.n_gwas, matched_snps=int(len(ss)), n_models=n,
               sumstats_md5=hashlib.md5(open(a.sumstats, "rb").read()).hexdigest()), open(f"{a.outdir}/application_ldsc.json", "w"), indent=1)

fig, ax = plt.subplots(figsize=(6.2, 5.6))
exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
for lab, col, c_ in [("Standard TWAS, EAS LD", "p_standard_target", "#2a78d6"), ("Adapted variance control (EUR q), EUR LD", "p_vc_eur_eurld", "#1baf7a"),
                     ("XCAL-t (GWAS SDs; recommended)", "p_xcalt", "#e34948")]:
    ax.plot(exp, np.sort(-np.log10(np.maximum(R[col].values, 1e-300)))[::-1], ".", ms=2.5, color=c_, label=lab)
ax.plot([0, exp.max()], [0, exp.max()], color="#555", lw=1); ax.axhline(-np.log10(bonf), color="#999", ls="--", lw=1)
ax.set_xlabel("Expected −log10 p"); ax.set_ylabel("Observed −log10 p"); ax.legend(frameon=False, fontsize=9)
ax.set_ylim(0, min(60, ax.get_ylim()[1]))
fig.tight_layout(); fig.savefig(f"{a.outdir}/fig_application_qq.png", dpi=300); plt.close(fig)
say(S.to_string()); say(ratio.round(3).to_string())
