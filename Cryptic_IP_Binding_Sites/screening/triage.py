import sys
import argparse
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns


def triage_results(organism: str):
    """Triages hit results and generates a summary report/plots."""
    results_dir = Path(f"data/results/{organism}")
    csv_path = results_dir / f"{organism}_master.csv"

    if not csv_path.exists():
        print(f"Error: {csv_path} not found.")
        sys.exit(1)

    df = pd.read_csv(csv_path)

    # Filter for candidates
    hits = df[df["ManualReviewFlag"]].sort_values(by="CompositeScore", ascending=False)

    print(f"Total proteins screened: {df['UniProtID'].nunique()}")
    print(f"Total pockets analyzed: {len(df)}")
    print(f"Total flagged candidates: {len(hits)}")

    # Save top hits
    hits.to_csv(results_dir / f"{organism}_top_candidates.csv", index=False)

    # Generate basic plots
    if len(df) > 0:
        plt.figure(figsize=(10, 6))
        sns.histplot(df["CompositeScore"], bins=50)
        plt.axvline(x=0.7, color="r", linestyle="--", label="Candidate Threshold")
        plt.title(f"{organism.capitalize()} Proteome: Composite Score Distribution")
        plt.xlabel("Composite Score")
        plt.legend()
        plt.savefig(results_dir / f"{organism}_score_distribution.png")

    if len(hits) > 0:
        plt.figure(figsize=(10, 6))
        sns.scatterplot(
            data=hits,
            x="MeanSASA",
            y="CompositeScore",
            hue="MeanPLDDT",
            palette="viridis",
        )
        plt.title(f"{organism.capitalize()} Top Candidates: Score vs SASA")
        plt.savefig(results_dir / f"{organism}_hits_scatter.png")

    print("Triage complete. Files saved.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--organism", required=True)
    args = parser.parse_args()

    triage_results(args.organism)
