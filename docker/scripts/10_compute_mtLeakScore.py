import argparse
import json
from pathlib import Path
import yaml
import numpy as np
import pandas as pd
from scipy.stats.mstats import winsorize


EPS = 1e-6


def compute_mtLeakScore(mtUMI, oxphos):
    mt_term = np.log2(mtUMI + 1 + EPS)
    ox_term = (oxphos - np.nanmean(oxphos)) / (np.nanstd(oxphos) + EPS)
    score = mt_term * (-ox_term)
    score = winsorize(score, limits=[0.005, 0.005])  # 0.5% trimming
    return score


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    base = Path(cfg["base_dir"])
    proc_dir = base / cfg["processed_data_dir"]
    out_dir = base / "results/mtLeakScore"
    out_dir.mkdir(parents=True, exist_ok=True)

    # -------------------------------------
    # Load needed matrices
    # -------------------------------------
    mtUMI_path = proc_dir / "mtUMI.tsv"
    oxphos_path = proc_dir / "OXPHOS_score.tsv"

    mtUMI = pd.read_csv(mtUMI_path, sep="\t", index_col=0)
    oxphos = pd.read_csv(oxphos_path, sep="\t", index_col=0)

    # Align indices
    mtUMI, oxphos = mtUMI.align(oxphos, join="inner", axis=0)

    score = compute_mtLeakScore(mtUMI.values.flatten(),
                                oxphos.values.flatten())

    df = pd.DataFrame({
        "sample": mtUMI.index,
        "mtLeakScore": score
    })
    df.to_csv(out_dir / "mtLeakScore.tsv", sep="\t", index=False)

    # metadata
    meta = {
        "formula": "log2(mtUMI+1) * -(scaled OXPHOS)",
        "winsorize": "0.5% each side",
        "EPS": EPS,
    }
    with open(out_dir / "score_method.json", "w") as f:
        json.dump(meta, f, indent=2)

    print("[OK] mtLeakScore saved.")


if __name__ == "__main__":
    main()
