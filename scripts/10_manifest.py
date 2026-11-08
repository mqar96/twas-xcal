import hashlib, json, os, platform, sys
import importlib.metadata as md

outdir = sys.argv[1]
files = {}
for root, _, names in os.walk(outdir):
    for n in sorted(names):
        if n == "MANIFEST.json":
            continue
        p = os.path.join(root, n)
        files[os.path.relpath(p, outdir)] = hashlib.sha256(open(p, "rb").read()).hexdigest()
pk = {}
for k in ["numpy", "scipy", "pandas", "scikit-learn", "msprime", "stdpopsim", "matplotlib", "joblib", "bed-reader", "twas-xcal"]:
    try:
        pk[k] = md.version(k)
    except md.PackageNotFoundError:
        pk[k] = None
json.dump(dict(python=platform.python_version(), packages=pk, files=files), open(os.path.join(outdir, "MANIFEST.json"), "w"), indent=1)
print("manifest:", len(files), "files")
