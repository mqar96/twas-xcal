import warnings
import numpy as np

POPS = ("YRI", "CEU", "CHB")
DEFAULT_SPLITS = {"train": 1000, "ldref": 500, "gwas": 20000}


def simulate_region(seed, splits=None, length=1_000_000, max_snps=400, maf_any=0.05, maf_all=0.01):
    import stdpopsim
    warnings.filterwarnings("ignore")
    splits = dict(DEFAULT_SPLITS if splits is None else splits)
    n_per_pop = sum(splits.values())
    sp = stdpopsim.get_species("HomSap")
    model = sp.get_demographic_model("OutOfAfrica_3G09")
    contig = sp.get_contig(length=length)
    ts = stdpopsim.get_engine("msprime").simulate(model, contig, {p: n_per_pop for p in POPS}, seed=seed)
    pop_of_sample = np.array([ts.node(n).population for n in ts.samples()])
    names = [ts.population(i).metadata.get("name", str(i)) for i in range(ts.num_populations)]
    idx = {p: np.where(pop_of_sample == names.index(p))[0] for p in POPS}
    ok_sites = []
    for k, var in enumerate(ts.variants()):
        g = var.genotypes
        if g.max() > 1:
            continue
        f = np.array([g[idx[p]].mean() for p in POPS])
        maf = np.minimum(f, 1 - f)
        if maf.max() >= maf_any and maf.min() >= maf_all:
            ok_sites.append(k)
    rng = np.random.default_rng(seed + 7)
    if len(ok_sites) > max_snps:
        ok_sites = sorted(rng.choice(ok_sites, max_snps, replace=False).tolist())
    want = set(ok_sites)
    H = np.zeros((len(ok_sites), ts.num_samples), np.int8)
    pos = np.zeros(len(ok_sites))
    j = 0
    for k, var in enumerate(ts.variants()):
        if k in want:
            H[j] = var.genotypes
            pos[j] = var.site.position
            j += 1
    res = {"pos": pos}
    for p in POPS:
        h = H[:, idx[p]]
        dip = (h[:, 0::2] + h[:, 1::2]).T
        o = 0
        for name, n in splits.items():
            res[f"{p}_{name}"] = dip[o:o + n]
            o += n
    res["n_snps"] = np.array(len(ok_sites))
    return res
