import argparse, glob, hashlib, json, os, time
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import norm, t as tdist
from joblib import Parallel, delayed
from xcal.calibration import Panel
from xcal.core import satterthwaite_df

ap = argparse.ArgumentParser(description="Null benchmark: counted false positives from simulated GWAS z-scores with exactly known covariance")
ap.add_argument("genodir"); ap.add_argument("weightdir"); ap.add_argument("outdir")
ap.add_argument("--reps", type=int, default=30000)
ap.add_argument("--jobs", type=int, default=2)
ap.add_argument("--max-regions", type=int, default=None)
ap.add_argument("--gc-lambda", default=None, help="CSV of genomic-control lambdas per cell (from 05); runs the exact GC counting pass")
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

POPS = ("YRI", "CEU", "CHB")
C_DENSE = [0.0, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8]
SPARSE = [(pi, c) for pi in (0.1, 0.01, 0.001) for c in (0.1, 0.4)]
CELLS = [(1.0, c) for c in C_DENSE] + SPARSE
ALPHAS = np.array([0.05, 1e-3, 1e-4, 2.5e-6])
ZA = norm.isf(ALPHAS / 2)
# numerator columns per gene: 0 = v_G (GWAS SDs), 1 = v_T (target-panel SDs), 2 = v_E (EUR-panel SDs)
METHODS = ["ORACLE", "STD_M", "STD_E", "VC_E_E", "VC_E_M", "XCAL", "XCALT", "XCALT_PSD"]
GRID = np.logspace(-3, 3.5, 2000)
SCRIPT_SHA = hashlib.sha256(open(__file__, "rb").read()).hexdigest()
LAMBDA = None
if a.gc_lambda:
    lt = pd.read_csv(a.gc_lambda)
    LAMBDA = {(r.gwas_pop, float(r.pi), float(r.c)): float(r.lambda_gc) for r in lt.itertuples()}


def sha(f):
    return hashlib.sha256(open(f, "rb").read()).hexdigest()


def psd_sqrt(C):
    lam, U = np.linalg.eigh((C + C.T) / 2)
    return U * np.sqrt(np.clip(lam, 0, None))


def config_hash(rid, wf):
    cfg = dict(script=SCRIPT_SHA, reps=a.reps, cells=CELLS, alphas=ALPHAS.tolist(), geno=sha(f"{a.genodir}/region_{rid:04d}.npz"),
               weights=sha(wf), gc=sorted((k[0], k[1], k[2], v) for k, v in LAMBDA.items()) if LAMBDA else None)
    return hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()[:16]


def region(wf):
    rid = int(os.path.basename(wf)[8:12])
    tag = "gc_" if LAMBDA else ""
    of = f"{a.outdir}/{tag}region_{rid:04d}.csv.gz"
    h = config_hash(rid, wf)
    if os.path.exists(of) and os.path.exists(of + ".hash") and open(of + ".hash").read() == h:
        return of
    recs = [r for r in json.loads(str(np.load(wf)["recs"])) if r["n_w"] > 0]
    d = np.load(f"{a.genodir}/region_{rid:04d}.npz")
    p = int(d["n_snps"])
    T = {g: Panel(d[f"{g}_ldref"]) for g in POPS}
    rng = np.random.default_rng(7000 + rid)
    rows, grids = [], {}
    for g in POPS:
        Rg = d[f"{g}_R_gwas"].astype(np.float64)
        sG = d[f"{g}_sd_gwas"]
        Pt, Pe = T[g], T["CEU"]
        genes = []
        for r in recs:
            w = np.zeros(p); w[r["idx"]] = r["w"]
            vG, vT, vE = w * sG, Pt.v(w), Pe.v(w)
            mG, mT, mE = Pt.moments(vG), Pt.moments(vT), Pe.moments(vE)
            if min(mG["L"], mT["L"], mE["L"]) <= 1e-12:
                continue
            LbG, QbG = Pt.jackknife(vG, K=20, seed=rid)
            LbT, QbT = Pt.jackknife(vT, K=20, seed=rid)
            aG = Rg @ vG
            genes.append(dict(r=r, V=(vG, vT, vE), mG=mG, mT=mT, mE=mE, jkG=(LbG, QbG), jkT=(LbT, QbT),
                              L_true=float(vG @ aG), Q_true=float(aG @ aG)))
        if not genes:
            continue
        m = len(genes)
        Vm = np.column_stack([x for gg in genes for x in gg["V"]])
        A = Rg @ Vm
        C0 = Vm.T @ A
        C1 = A.T @ A
        for pi, c in CELLS:
            if pi == 1.0:
                num = rng.standard_normal((a.reps, 3 * m)) @ psd_sqrt(C0 + c * C1).T
            else:
                noise = rng.standard_normal((a.reps, 3 * m)) @ psd_sqrt(C0).T
                mask = rng.random((a.reps, p)) < pi
                rr, cc = np.nonzero(mask)
                vals = rng.normal(0.0, np.sqrt(c / pi), len(rr))
                num = noise + sparse.csr_matrix((vals, (rr, cc)), shape=(a.reps, p)) @ A
            gridcount = np.zeros(len(GRID))
            for k, gg in enumerate(genes):
                i0 = 3 * k
                var = [C0[i0 + j, i0 + j] + c * C1[i0 + j, i0 + j] for j in range(3)]
                LG, QG = gg["mG"]["L"], gg["mG"]["Q"]
                LT, QT = gg["mT"]["L"], gg["mT"]["Q"]
                LE, qE = gg["mE"]["L"], gg["mE"]["Q"] / gg["mE"]["L"]
                dfG = float(satterthwaite_df(gg["jkG"][0] + c * gg["jkG"][1], LG + c * QG))
                dfT = float(satterthwaite_df(gg["jkT"][0] + c * gg["jkT"][1], LT + c * QT))
                spec = {"ORACLE": (0, np.sqrt(var[0]), None), "STD_M": (1, np.sqrt(LT), None), "STD_E": (2, np.sqrt(LE), None),
                        "VC_E_E": (2, np.sqrt(LE * (1 + c * qE)), None), "VC_E_M": (1, np.sqrt(LT * (1 + c * qE)), None),
                        "XCAL": (0, np.sqrt(LG + c * QG), None), "XCALT": (0, np.sqrt(LG + c * QG), dfG),
                        "XCALT_PSD": (1, np.sqrt(LT + c * QT), dfT)}
                r = gg["r"]
                row = dict(gene=r["gene"], region=rid, gwas_pop=g, pi=pi, c=c, reps=a.reps)
                if LAMBDA:
                    lam = LAMBDA[(g, pi, c)]
                    stat = np.abs(num[:, i0 + 1]) / np.sqrt(LT)
                    for ai, al in enumerate(ALPHAS):
                        row[f"n_GC_{al:g}"] = int(np.count_nonzero(stat > ZA[ai] * np.sqrt(lam)))
                    rows.append(row)
                    continue
                row.update(train_pop=r["train_pop"], n_eqtl=r["n_eqtl"], n_w=r["n_w"], r2_cv=r.get("r2_cv"), r2_target=r.get(f"r2_{g}"),
                           L_G=LG, Q_G=QG, Q_G_untrunc=gg["mG"]["Q_untruncated"], Q_G_raw=gg["mG"]["Q_raw"], floor_G=gg["mG"]["floor_active"],
                           L_true=gg["L_true"], Q_true=gg["Q_true"], q_true=gg["Q_true"] / gg["L_true"],
                           qG_panel=QG / LG, qT=QT / LT, qE=qE, floor_T=gg["mT"]["floor_active"], floor_E=gg["mE"]["floor_active"],
                           df=dfG, df_psd=dfT)
                for meth in METHODS:
                    j, sd, df = spec[meth]
                    stat = np.abs(num[:, i0 + j]) / sd
                    V = var[j] / sd ** 2
                    thr = tdist.isf(ALPHAS / 2, df) if df is not None else ZA
                    row[f"V_{meth}"] = V
                    for ai, al in enumerate(ALPHAS):
                        row[f"n_{meth}_{al:g}"] = int(np.count_nonzero(stat > thr[ai]))
                        row[f"pred_{meth}_{al:g}"] = float(2 * norm.sf(thr[ai] / np.sqrt(V)))
                    if meth == "STD_M":
                        s2 = np.sort(stat ** 2)
                        gridcount += len(s2) - np.searchsorted(s2, GRID, side="right")
                rows.append(row)
            if not LAMBDA:
                grids[f"{g}|{pi}|{c}"] = gridcount
    pd.DataFrame(rows).to_csv(of + ".tmp", index=False, compression="gzip")
    if not LAMBDA:
        np.savez_compressed(f"{a.outdir}/grid_{rid:04d}.npz", **grids)
    os.replace(of + ".tmp", of)
    open(of + ".hash", "w").write(h)
    return of


files = sorted(glob.glob(f"{a.weightdir}/weights_*.npz"))[: a.max_regions]
t = time.time()
Parallel(n_jobs=a.jobs)(delayed(region)(f) for f in files)
print(("GC pass: " if LAMBDA else "") + f"benchmarked {len(files)} regions in {round(time.time() - t)} s")
