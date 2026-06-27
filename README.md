# zephyr-taxprofile

## Install

```bash
pip install -r requirements.txt
```

## Quickstart

```bash
# 1. Download the reference panel from NCBI, then assemble it.
bash references/fetch_references.sh      
python references/build_reference.py

# 2. Run the batch
python scripts/run_batch.py --input data/ --config config.yaml

