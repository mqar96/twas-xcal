import argparse, hashlib, os, sys

ap = argparse.ArgumentParser(description="Verify 1000 Genomes PLINK files against the upstream Zenodo MD5 list")
ap.add_argument("geno_dir")
a = ap.parse_args()
bad = 0
for line in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "zenodo_md5.txt")):
    md5, name = line.split()
    h = hashlib.md5()
    with open(os.path.join(a.geno_dir, name), "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 24), b""):
            h.update(chunk)
    ok = h.hexdigest() == md5
    bad += not ok
    print(("OK  " if ok else "BAD ") + name)
sys.exit(1 if bad else 0)
