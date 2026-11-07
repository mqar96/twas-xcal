# Real-data pipeline (phase 2)

1. **Download 1000 Genomes phase 3 PLINK files** (Zenodo 10.5281/zenodo.6614170, CC BY 4.0): `1000G_{AFR,EUR,EAS,SAS,AMR}.{bed,bim,fam}`. Verify them against `zenodo_md5.txt`.
2. **Get the MESA models:** `git clone https://github.com/WheelerLab/DivPop`, then from a separate directory run `python -I read_divpop.py /path/to/DivPop`.
3. **Run the calibration:** `python real_calibration.py out/real` (about 40 min on 2 cores; peak RAM about 2 GB; resumable per chromosome). Edit the `G` path at the top first.
4. **Analyse:** `python analyze_real.py out/real results/real`.
