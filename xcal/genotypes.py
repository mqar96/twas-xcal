import numpy as np
import stdpopsim
import warnings

POPS = ("YRI", "CEU", "CHB")
SPLITS = {"train": 1000, "ldref": 500, "gwas": 1100}


def simulate_region(seed, length=1_000_000, n_per_pop=2600, max_snps=400, maf_any=0.05, maf_all=0.01):
    warnings.filterwarnings("ignore")
    sp = stdpopsim.get_species("HomSap")
    model = sp.get_demographic_model("OutOfAfrica_3G09")
    contig = sp.get_contig(length=length)
    eng = stdpopsim.get_engine("msprime")
    ts = eng.simulate(model, contig, {p: n_per_pop for p in POPS}, seed=seed)
    G = ts.genotype_matrix().astype(np.int8)
    keep = np.array([len(v.alleles) == 2 for v in ts.variants()]) if False else (G.max(axis=1) <= 1)
    G = G[keep]
    pop_of_sample = np.array([ts.node(n).population for n in ts.samples()])
    names = [ts.population(i).metadata.get("name", str(i)) for i in range(ts.num_populations)]
    out = {}
    for p in POPS:
        pid = names.index(p)
        idx = np.where(pop_of_sample == pid)[0]
        h = G[:, idx]
        dip = (h[:, 0::2] + h[:, 1::2]).T
        out[p] = dip
    freqs = np.vstack([out[p].mean(0) / 2 for p in POPS])
    maf = np.minimum(freqs, 1 - freqs)
    ok = (maf.max(0) >= maf_any) & (maf.min(0) >= maf_all)
    sel = np.where(ok)[0]
    rng = np.random.default_rng(seed + 7)
    if len(sel) > max_snps:
        sel = np.sort(rng.choice(sel, max_snps, replace=False))
    res = {}
    for p in POPS:
        X = out[p][:, sel]
        o = 0
        for k, n in SPLITS.items():
            res[f"{p}_{k}"] = X[o:o + n].astype(np.int8)
            o += n
    res["n_snps"] = np.array(len(sel))
    return res


def standardize(X):
    X = X.astype(np.float64)
    mu = X.mean(0)
    sd = X.std(0)
    sd[sd == 0] = 1.0
    return (X - mu) / sd, mu, sd


def ld(X):
    Z, _, sd = standardize(X)
    R = Z.T @ Z / Z.shape[0]
    mono = X.std(0) == 0
    R[mono, :] = 0
    R[:, mono] = 0
    R[mono, mono] = 1.0
    return R
