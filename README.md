# zephyr-taxprofile

## Install

```bash
# core only (QC + align + coverage + confidence + figures on local FASTAs)
pip install -r requirements.txt

# add discovery tools (kraken2, sourmash) + optional CLI binaries
mamba env create -f environment.yml && mamba activate zephyr-taxprofile
```

## Quickstart

```bash
# 1. Download the reference panel from NCBI (one file per virus -> references/raw/),
#    then assemble it. Panel matches SecureBio's 14 "viruses detected" rows.
bash references/fetch_references.sh      # needs ncbi-datasets-cli; or download manually
python references/build_reference.py

# 2. Put your pool FASTAs (.fasta or .fasta.gz) in a directory, e.g. data/
#    e.g. 250107-BoDT-NAS-P1.respiratory.fasta.gz

# 3. Point config.yaml at your kraken2 / sourmash DBs (optional; blank = align-only)

# 4. Run the batch
python scripts/run_batch.py --input data/ --config config.yaml

# 5. (optional) open notebooks/analysis.ipynb for the write-up
```

Verify the install with the bundled synthetic test (no external data):

```bash
pytest tests/test_coverage.py -q
```

