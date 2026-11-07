import sys, os, glob, json, time, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from xcal.calibration import panel_stats, v_std, poly_q, v_xcal, jackknife_panels, LQ_from_panels, satterthwaite_df, quad_R2_unbiased
from joblib import Parallel, delayed

gdir, wdir, out = sys.argv[1], sys.argv[2], sys.argv[3]
POPS = ("YRI", "CEU", "CHB")
CS = [0.0, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8]


def region(wf):
    rid = int(os.path.basename(wf)[8:12])
    recs = [r for r in json.loads(str(np.load(wf)["recs"])) if r["n_w"] > 0]
    if not recs:
        return []
    d = np.load(f"{gdir}/region_{rid:04d}.npz")
    P = {f"{p}_{k}": panel_stats(d[f"{p}_{k}"]) for p in POPS for k in ("ldref", "gwas")}
    p = d["YRI_gwas"].shape[1]
    JK = {g: jackknife_panels(d[f"{g}_ldref"], K=20, seed=rid) for g in POPS}
    rows = []
    for r in recs:
        w = np.zeros(p); w[r["idx"]] = r["w"]
        for g in POPS:
            Pg = P[f"{g}_gwas"]; Pt = P[f"{g}_ldref"]; Pe = P["CEU_ldref"]
            qE = poly_q(w, P["CEU_ldref"]); qT = poly_q(w, Pt)
            Lk, Qk = LQ_from_panels(w, JK[g])
            vt = w * Pt["s"]; Lf = vt @ Pt["R"] @ vt; Qf = quad_R2_unbiased(Pt, vt)
            for c in CS:
                Vm = v_std(w, Pt, Pg, c)
                Ve = v_std(w, Pe, Pg, c)
                rows.append(dict(gene=r["gene"], region=rid, train_pop=r["train_pop"], n_eqtl=r["n_eqtl"], h2=r["h2"], k=r["k"],
                                 n_w=r["n_w"], r2_cv=r.get("r2_cv"), r2_target=r.get(f"r2_{g}"), gwas_pop=g, c=c,
                                 V_unadj_matchedLD=Vm, V_unadj_EURLD=Ve,
                                 V_VCeur_matchedLD=Vm / (1 + c * qE), V_VCeur_EURLD=Ve / (1 + c * qE),
                                 V_VCtarget_matchedLD=Vm / (1 + c * qT),
                                 V_XCAL=v_xcal(w, Pt, Pg, c), df_XCAL=satterthwaite_df(Lk, Qk, c, Lf + c * Qf), qE=qE, qT=qT))
    return rows

files = sorted(glob.glob(f"{wdir}/weights_*.npz"))
t = time.time()
res = Parallel(n_jobs=2)(delayed(region)(f) for f in files)
df = pd.DataFrame([x for r in res for x in r])
df.to_parquet(out) if out.endswith(".parquet") else df.to_csv(out, index=False)
print("rows", len(df), "genes", df.gene.nunique(), "time", round(time.time() - t, 1))
