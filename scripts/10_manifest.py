"""Write results/MANIFEST.json: SHA-256 of every output plus the provenance records of the steps that made them."""
import glob, hashlib, json, os, sys

outdir = sys.argv[1]
files = {}
for root, _, names in os.walk(outdir):
    for n in sorted(names):
        if n == "MANIFEST.json":
            continue
        p = os.path.join(root, n)
        files[os.path.relpath(p, outdir).replace(os.sep, "/")] = hashlib.sha256(open(p, "rb").read()).hexdigest()
steps = {}
for f in sorted(glob.glob(os.path.join(outdir, "provenance", "*.json"))):
    r = json.load(open(f, encoding="utf-8"))
    steps[r["step"]] = {k: r[k] for k in ("command", "start", "end", "exit_code", "git_revision", "git_dirty", "xcal_module",
                                          "xcal_version", "installed_twas_xcal", "python") if k in r}
json.dump(dict(steps=steps, files=files), open(os.path.join(outdir, "MANIFEST.json"), "w", encoding="utf-8"), indent=1)
print("manifest:", len(files), "files,", len(steps), "steps")
