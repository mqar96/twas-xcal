import argparse, glob, json, os, time
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import norm, t as tdist
from joblib import Parallel, delayed
from xcal.calibration import Panel
from xcal.core import satterthwaite_df

ap = argparse.ArgumentParser(description="Null benchmark: counted false positives from simulated GWAS z-scores with exact covariance")
ap.add_argument("genodir"); ap.add_argument("weightdir"); ap.add_argument("outdir")
ap.add_argument("--reps", type=int, default=30000)
ap.add_argument("--jobs", type=int, default=2)
ap.add_argument("--max-regions", type=int, default=None)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

POPS = ("YRI", "CEU", "CHB")
C_DENSE = [0.0, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8]
SPARSE = [(pi, c) for pi in (0.1, 0.01, 0.001) for c in (0.1, 0.4)]
CELLS = [(1.0, c) for c in C_DENSE] + SPARSE
ALPHAS = np.array([0.05, 1e-3, 1e-4, 2.5e-6])
ZA = norm.isf(ALPHAS / 2)
METHODS = ["ORACLE", "STD_M", "STD_E", "VC_E_E", "VC_E_M", "XCAL", "XCALT"]
GRID = np.logspace(-3, 3.5, 2000)


def psd_sqrt(C):
    lam, U = np.linalg.eigh((C + C.T) / 2)
    return U * np.sqrt(np.clip(lam, 0, None))


def region(wf):
    rid = int(os.path.basename(wf)[8:12])
    of = f"{a.outdir}/region_{rid:04d}.csv.gz"
    if os.path.exists(of):
        return of
    recs = [r for r in json.loads(str(np.load(wf)["recs"])) if r["n_w"] > 0]
    if not recs:
        pd.DataFrame().to_csv(of, index=False)
        return of
    d = np.load(f"{a.genodir}/region_{rid:04d}.npz")
    p = int(d["n_snps"])
    T = {g: Panel(d[f"{g}_ldref"]) for g in POPS}
    rng = np.random.default_rng(7000 + rid)
    rows, grids = [], {}
    for g in POPS:
        Rg = d[f"{g}_R_gwas"].astype(np.float64)
        Pt, Pe = T[g], T["CEU"]
        genes = []
        for r in recs:
            w = np.zeros(p); w[r["idx"]] = r["w"]
            vT, vE = Pt.v(w), Pe.v(w)
            mT, mE = Pt.moments(w), Pe.moments(w)
            if mT["L"] <= 1e-12 or mE["L"] <= 1e-12:
                continue
            Lb, Qb = Pt.jackknife(w, K=20, seed=rid)
            aT = Rg @ vT
            genes.append(dict(r=r, vT=vT, vE=vE, mT=mT, mE=mE, Lb=Lb, Qb=Qb, L_true=float(vT @ aT), Q_true=float(aT @ aT),
                              q_true=float(aT @ aT) / float(vT @ aT)))
        if not genes:
            continue
        m = len(genes)
        Vm = np.column_stack([x for gg in genes for x in (gg["vT"], gg["vE"])])
        A = Rg @ Vm
        C0 = Vm.T @ A
        C1 = A.T @ A
        for pi, c in CELLS:
            den = {}
            for k, gg in enumerate(genes):
                LT, QT, LE = gg["mT"]["L"], gg["mT"]["Q"], gg["mE"]["L"]
                qE = gg["mE"]["Q"] / LE
                varT = C0[2 * k, 2 * k] + c * C1[2 * k, 2 * k]
                varE = C0[2 * k + 1, 2 * k + 1] + c * C1[2 * k + 1, 2 * k + 1]
                df = float(satterthwaite_df(gg["Lb"] + c * gg["Qb"], LT + c * QT))
                den[k] = dict(ORACLE=(2 * k, np.sqrt(varT)), STD_M=(2 * k, np.sqrt(LT)), STD_E=(2 * k + 1, np.sqrt(LE)),
                              VC_E_E=(2 * k + 1, np.sqrt(LE * (1 + c * qE))), VC_E_M=(2 * k, np.sqrt(LT * (1 + c * qE))),
                              XCAL=(2 * k, np.sqrt(LT + c * QT)), XCALT=(2 * k, np.sqrt(LT + c * QT)),
                              _df=df, _var=(varT, varT, varE, varE, varT, varT, varT))
            if pi == 1.0:
                num = rng.standard_normal((a.reps, 2 * m)) @ psd_sqrt(C0 + c * C1).T
            else:
                noise = rng.standard_normal((a.reps, 2 * m)) @ psd_sqrt(C0).T
                mask = rng.random((a.reps, p)) < pi
                rr, cc = np.nonzero(mask)
                vals = rng.normal(0.0, np.sqrt(c / pi), len(rr))
                U = sparse.csr_matrix((vals, (rr, cc)), shape=(a.reps, p))
                num = noise + U @ A
            gridcount = np.zeros(len(GRID))
            for k, gg in enumerate(genes):
                r = gg["r"]; dk = den[k]; df = dk["_df"]
                row = dict(gene=r["gene"], region=rid, train_pop=r["train_pop"], n_eqtl=r["n_eqtl"], n_w=r["n_w"],
                           r2_cv=r.get("r2_cv"), r2_target=r.get(f"r2_{g}"), gwas_pop=g, pi=pi, c=c, reps=a.reps,
                           L_T=gg["mT"]["L"], Q_T=gg["mT"]["Q"], Q_T_untrunc=gg["mT"]["Q_untruncated"],
                           floor_T=gg["mT"]["floor_active"], floor_E=gg["mE"]["floor_active"],
                           qT=gg["mT"]["Q"] / gg["mT"]["L"], qE=gg["mE"]["Q"] / gg["mE"]["L"], q_true=gg["q_true"], L_true=gg["L_true"], Q_true=gg["Q_true"],
                           Q_T_raw=gg["mT"]["Q_raw"], df=df)
                for mi, meth in enumerate(METHODS):
                    col, sd = dk[meth]
                    stat = np.abs(num[:, col]) / sd
                    V = dk["_var"][mi] / sd ** 2
                    thr = tdist.isf(ALPHAS / 2, df) if meth == "XCALT" else ZA
                    row[f"V_{meth}"] = V
                    for ai, al in enumerate(ALPHAS):
                        row[f"n_{meth}_{al:g}"] = int(np.count_nonzero(stat > thr[ai]))
                        row[f"pred_{meth}_{al:g}"] = float(2 * norm.sf(thr[ai] / np.sqrt(V)))
                    if meth == "STD_M":
                        s2 = np.sort(stat ** 2)
                        gridcount += len(s2) - np.searchsorted(s2, GRID, side="right")
                rows.append(row)
            grids[f"{g}|{pi}|{c}"] = gridcount
    pd.DataFrame(rows).to_csv(of + ".tmp", index=False, compression="gzip")
    np.savez_compressed(f"{a.outdir}/grid_{rid:04d}.npz", **{k: v for k, v in grids.items()})
    os.replace(of + ".tmp", of)
    return of


files = sorted(glob.glob(f"{a.weightdir}/weights_*.npz"))[: a.max_regions]
t = time.time()
Parallel(n_jobs=a.jobs)(delayed(region)(f) for f in files)
print("benchmarked", len(files), "regions in", round(time.time() - t), "s")
