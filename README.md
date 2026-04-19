# IDPC_Reproduction

Reproduction repository for the paper:

**Intersection-Defined Phase Coordinates Reveal Localized Selection and a Non-Closed Observational Structure**  
https://doi.org/10.5281/zenodo.19628769

This repository provides the full pipeline required to reproduce the experiments and figures reported in the paper.


## Quick Demo (Browser)

Run the Binder demo:

https://bids.mybinder.org/v2/gh/SATORU-SIEL/IDPC_Reproduction/HEAD?urlpath=%2Fdoc%2Ftree%2FIDPC_Repro_Demo.ipynb

## Repository Structure
```
IDPC_Reproduction/

├── IDPC_Repro_Demo.ipynb        # Binder demo

├── IDPC_Reproduction.ipynb      # Full reproduction notebook

├── Quantum Measurement.ipynb    # Quantum measurement notebook

├── requirements.txt             # Dependencies

├── scripts/                     # Reproduction and control-analysis scripts

├── reports/                     # Reproduction reports, metrics, and figures

├── IDPC_Reproduction/           # CSV files (demo)

├── IDPC_Reproduction_ricci/     # CSV files (demo)
```

## Key Reproduction Results

The bundled demo CSV files reproduce the main reported metrics:

| Metric | Reproduced value |
| --- | ---: |
| switch_gain (best intersection model) | 0.809793 |
| block-permutation null mean | 0.564516 |
| empirical permutation p | 0.004975 |
| LOSO state-classification AUC | 0.791875 |
| LOSO state-classification accuracy | 0.695122 |
| mean absolute Pearson | 0.419253 |

See `reports/reproduction_analysis_report.md` for the basic reproduction report.

Additional control analyses are included:

- `reports/mismatch_pipeline_report_en.md`: full-pipeline mismatched-pairing control for the Section 7.2 switch_gain claim
- `reports/sic_prediction_report_en.md`: SIC boundary-event prediction experiment
- `reports/random_noise_experiment_report_en.md`: random-noise control experiment

## Quick Local Reproduction

Install dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt jupyter
```

Run the demo notebook:

```bash
jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=3600 \
  --output IDPC_Repro_Demo_executed.ipynb IDPC_Repro_Demo.ipynb
```

Run the included script-based analyses:

```bash
python scripts/sic_prediction.py
python scripts/sic_prediction_figures.py
python scripts/mismatch_pipeline_collect.py
python scripts/mismatch_pipeline_figure.py
```

## Full Reproduction

Raw EEG dataset:

**Raw EEG Data (26 Sessions) for Reproduction of Neural–Quantum Structural Analysis in IDPC**

https://doi.org/10.5281/zenodo.19624924

To reproduce all results:

1. Download the dataset from Zenodo  
2. Run `IDPC_Reproduction.ipynb`


## Notes

- CSV files included in this repository are for demo execution
- Full reproduction requires the raw EEG dataset
- Quantum Measurement.ipynb is required only when reproducing the quantum experiment
