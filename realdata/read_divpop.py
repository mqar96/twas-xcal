import argparse, hashlib, os, sqlite3
import pandas as pd

DIVPOP_COMMIT = "85100c2a41e2022e1e684e165980c3e9868db925"
GIT_BLOB_SHA1 = {"AFA": "9371c3ab71ec431a539b7ea410d36bba87144730",
                 "CAU": "04aabe2e4001c3b2a2b168138feeb1f4c60c04ad",
                 "HIS": "1a71a6e5f0f4f01ed54316951115fae65e09d8b4"}

ap = argparse.ArgumentParser(description="Export MESA DivPop elastic-net models (pinned commit) to CSV")
ap.add_argument("divpop_dir", help=f"clone of https://github.com/WheelerLab/DivPop at commit {DIVPOP_COMMIT}")
ap.add_argument("outdir")
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
for p, want in GIT_BLOB_SHA1.items():
    f = f"{a.divpop_dir}/unfiltered_dbs/{p}_imputed_10_peer_3_pcs_2.db"
    data = open(f, "rb").read()
    got = hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()
    if got != want:
        raise SystemExit(f"{f}: git blob SHA-1 {got} != expected {want}; check out DivPop commit {DIVPOP_COMMIT}")
    c = sqlite3.connect(f"file:{f}?mode=ro", uri=True)
    w = pd.read_sql("select * from weights", c)
    e = pd.read_sql("select * from extra", c)
    w.to_csv(f"{a.outdir}/mesa_{p}_weights.csv.gz", index=False)
    e.to_csv(f"{a.outdir}/mesa_{p}_extra.csv", index=False)
    print(p, "genes", w.gene.nunique(), "weights", len(w), "verified")
