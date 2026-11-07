import sqlite3, sys, pandas as pd
base = sys.argv[1]
rows = []
for p in ["AFA", "CAU", "HIS", "ALL", "AFHI"]:
    c = sqlite3.connect(f"file:{base}/unfiltered_dbs/{p}_imputed_10_peer_3_pcs_2.db?mode=ro", uri=True)
    tabs = [r[0] for r in c.execute("select name from sqlite_master where type='table'")]
    w = pd.read_sql("select * from weights", c)
    e = pd.read_sql("select * from extra", c)
    rows.append(dict(model=p, tables=",".join(tabs), genes=w.gene.nunique(), weight_rows=len(w),
                     median_snps_per_gene=w.groupby("gene").size().median(), extra_cols=",".join(e.columns)))
    w.to_csv(f"mesa_{p}_weights.csv.gz", index=False); e.to_csv(f"mesa_{p}_extra.csv", index=False)
print(pd.DataFrame(rows).to_string())
