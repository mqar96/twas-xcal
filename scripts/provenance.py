"""Run one pipeline step and record its provenance.

Usage: python scripts/provenance.py WORKDIR STEP_NAME -- command [args...]
Writes WORKDIR/provenance/STEP_NAME.json with: command, start/end time, exit code, git revision and dirty state
of the repository, the imported xcal module path and version, SHA-256 of every source file, Python version and
the installed package list.
"""
import datetime, glob, hashlib, importlib.metadata as md, json, os, platform, subprocess, sys

workdir, step = sys.argv[1], sys.argv[2]
cmd = sys.argv[sys.argv.index("--") + 1:]
repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def git(*a):
    try:
        return subprocess.run(["git", "-C", repo, *a], capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return None


import xcal
src = {}
for pat in ("xcal/*.py", "scripts/*.py", "scripts/*.sh", "realdata/*.py"):
    for f in sorted(glob.glob(os.path.join(repo, pat))):
        src[os.path.relpath(f, repo)] = hashlib.sha256(open(f, "rb").read()).hexdigest()
rec = dict(step=step, command=cmd, cwd=os.getcwd(), start=datetime.datetime.now(datetime.timezone.utc).isoformat(),
           git_revision=git("rev-parse", "HEAD"), git_dirty=bool(git("status", "--porcelain", "--untracked-files=no")),
           xcal_module=os.path.abspath(xcal.__file__), xcal_version=xcal.__version__,
           installed_twas_xcal=md.version("twas-xcal") if "twas-xcal" in {d.metadata["Name"] for d in md.distributions()} else None,
           python=platform.python_version(), platform=platform.platform(),
           packages=sorted(f"{d.metadata['Name']}=={d.version}" for d in md.distributions()), source_sha256=src)
r = subprocess.run(cmd)
rec.update(end=datetime.datetime.now(datetime.timezone.utc).isoformat(), exit_code=r.returncode)
os.makedirs(os.path.join(workdir, "provenance"), exist_ok=True)
with open(os.path.join(workdir, "provenance", f"{step}.json"), "w", encoding="utf-8") as fh:
    json.dump(rec, fh, indent=1)
sys.exit(r.returncode)
