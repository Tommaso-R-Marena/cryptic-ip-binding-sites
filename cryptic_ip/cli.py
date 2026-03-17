"""
Command-line interface for cryptic IP site detection.
"""

import click
import json
from pathlib import Path

from .analysis import ProteinAnalyzer
from .validation import validate_adar2, ValidationSuite, StructureValidator, ResultsValidator
from .validation.md_validation import OpenMMMDValidationPipeline
from .database import ProteomeDownloader, ProteomeManager, DatabaseIntegrityChecker
from .reproducibility import deterministic_sort_dataframe, set_global_seed


@click.group()
@click.version_option(version="0.1.0")
def main():
    """
    Cryptic IP Binding Site Detection Pipeline

    Identify buried inositol phosphate binding sites in protein structures.
    """
    pass


@main.command()
@click.argument("pdb_file", type=click.Path(exists=True))
@click.option("--output", "-o", type=click.Path(), help="Output file for results")
@click.option("--score-threshold", "-t", default=0.60, help="Minimum score threshold")
@click.option("--use-ml-model", is_flag=True, help="Use trained ML classifier if available")
@click.option("--seed", type=int, default=42, show_default=True, help="Random seed")
@click.option(
    "--model-path",
    type=click.Path(exists=False),
    default="models/cryptic_ip_classifier_v1.pkl",
    help="Path to serialized ML model",
)
def analyze(pdb_file, output, score_threshold, use_ml_model, seed, model_path):
    """
    Analyze a single protein structure for cryptic IP binding sites.
    """
    click.echo(f"Analyzing {pdb_file}...\n")

    set_global_seed(seed)

    # Create analyzer
    analyzer = ProteinAnalyzer(pdb_file, use_ml_model=use_ml_model, model_path=model_path)

    # Detect and score pockets
    click.echo("Detecting pockets...")
    pockets = analyzer.detect_pockets()
    click.echo(f"Found {len(pockets)} pockets\n")

    click.echo("Scoring pockets...")
    scored = analyzer.score_all_pockets()

    # Filter by threshold
    candidates = scored[scored["composite_score"] >= score_threshold]

    click.echo(f"\nFound {len(candidates)} candidates above threshold {score_threshold}\n")

    if len(candidates) > 0:
        click.echo("Top candidates:")
        click.echo(
            candidates[["pocket_id", "composite_score", "volume", "sasa", "basic_residues"]]
            .head(10)
            .to_string()
        )

    # Save results
    if output:
        scored.to_csv(output, index=False)
        click.echo(f"\nResults saved to {output}")


@main.command()
@click.option(
    "--structure", type=click.Path(exists=True), help="Validate a single structure file (PDB/mmCIF)"
)
@click.option(
    "--results",
    "results_file",
    type=click.Path(exists=True),
    help="Validate a results CSV/JSON file",
)
@click.option(
    "--database",
    "database_path",
    type=click.Path(exists=True),
    help="Validate SQLite cache/database",
)
@click.option(
    "--all", "validate_all", is_flag=True, help="Run all selected validations in one pass"
)
@click.option(
    "--schema-file",
    multiple=True,
    type=click.Path(exists=True),
    help="Additional result files to compare schema consistency",
)
@click.option(
    "--checksums",
    type=click.Path(exists=True),
    help="JSON mapping of file path to expected SHA256 checksum",
)
@click.option(
    "--alphafold/--crystal",
    default=True,
    help="Use AlphaFold or crystal structure for legacy ADAR2 validation",
)
def validate(
    structure, results_file, database_path, validate_all, schema_file, checksums, alphafold
):
    """
    Validate structures, outputs, and database integrity.
    """
    run_legacy = not any([structure, results_file, database_path, validate_all])
    checksums_map = {}
    if checksums:
        checksums_map = json.loads(Path(checksums).read_text())

    if run_legacy:
        click.echo("Running ADAR2 validation...\n")
        legacy_results = validate_adar2(use_alphafold=alphafold)
        if legacy_results["validation_passed"]:
            click.secho("\n✓ Validation PASSED", fg="green", bold=True)
            return
        click.secho("\n✗ Validation FAILED", fg="red", bold=True)
        raise click.exceptions.Exit(1)

    errors_found = False

    def emit_report(name, report):
        nonlocal errors_found
        click.echo(f"\n{name} validation:")
        if report.metrics:
            click.echo(f"  Metrics: {report.metrics}")
        if report.issues:
            for issue in report.issues:
                color = "red" if issue.level == "error" else "yellow"
                click.secho(f"  [{issue.level.upper()}] {issue.check}: {issue.message}", fg=color)
                click.echo(f"    Suggested fix: {issue.suggested_fix}")
        else:
            click.secho("  No issues detected.", fg="green")
        errors_found = errors_found or any(issue.level == "error" for issue in report.issues)

    if structure or validate_all:
        if not structure:
            raise click.BadOptionUsage(
                "--all", "Provide --structure when using --all for structure validation."
            )
        struct_report = StructureValidator().validate(structure)
        emit_report("Structure", struct_report)

    if results_file or validate_all:
        if not results_file:
            raise click.BadOptionUsage(
                "--all", "Provide --results when using --all for results validation."
            )
        results_validator = ResultsValidator()
        res_report = results_validator.validate(results_file)
        emit_report("Results", res_report)
        if schema_file:
            schema_report = results_validator.validate_schema_consistency(
                (results_file, *schema_file)
            )
            emit_report("Schema consistency", schema_report)

    if database_path or validate_all:
        if not database_path:
            raise click.BadOptionUsage(
                "--all", "Provide --database when using --all for database validation."
            )
        db_report = DatabaseIntegrityChecker().validate(
            database_path, expected_checksums=checksums_map
        )
        emit_report("Database", db_report)

    if errors_found:
        raise click.exceptions.Exit(1)

    click.secho("\n✓ Validation checks completed without blocking errors", fg="green", bold=True)


@main.command()
def validate_suite():
    """
    Run full validation suite (positive and negative controls).
    """
    click.echo("Running complete validation suite...\n")
    suite = ValidationSuite()
    results = suite.run_full_validation()


@main.command()
@click.argument("organism", type=click.Choice(["yeast", "human", "dictyostelium"]))
@click.option("--data-dir", "-d", default="data/structures", help="Data directory")
def download(organism, data_dir):
    """
    Download AlphaFold proteome structures.
    """
    downloader = ProteomeDownloader(data_dir)

    # Show info
    info = downloader.get_info(organism)
    click.echo(f"\nDownloading {info['organism']}")
    click.echo(f"Proteins: {info['proteins']}")
    click.echo(f"Estimated size: {info['size_gb']} GB\n")

    if click.confirm("Continue with download?"):
        proteome_dir = downloader.download_proteome(organism)
        click.secho(f"\n✓ Download complete: {proteome_dir}", fg="green")


@main.command()
@click.argument("proteome_dir", type=click.Path(exists=True))
@click.option("--output", "-o", default="screening_results.csv", help="Output CSV file")
@click.option("--score-threshold", "-t", default=0.60, help="Minimum score threshold")
@click.option("--max-structures", "-n", type=int, help="Maximum structures to process")
@click.option("--use-ml-model", is_flag=True, help="Use trained ML classifier if available")
@click.option("--seed", type=int, default=42, show_default=True, help="Random seed")
@click.option(
    "--model-path",
    type=click.Path(exists=False),
    default="models/cryptic_ip_classifier_v1.pkl",
    help="Path to serialized ML model",
)
def screen(proteome_dir, output, score_threshold, max_structures, use_ml_model, seed, model_path):
    """
    Screen entire proteome for cryptic IP binding sites.
    """
    click.echo(f"Screening proteome: {proteome_dir}\n")
    set_global_seed(seed)

    # Initialize manager
    manager = ProteomeManager(proteome_dir)
    catalog = manager.build_catalog()

    if max_structures:
        catalog = catalog.head(max_structures)

    click.echo(f"Processing {len(catalog)} structures...\n")

    all_results = []

    with click.progressbar(catalog.iterrows(), length=len(catalog)) as bar:
        for idx, row in bar:
            try:
                pdb_path = row["filepath"]
                analyzer = ProteinAnalyzer(
                    pdb_path, use_ml_model=use_ml_model, model_path=model_path
                )
                scored = analyzer.score_all_pockets()

                # Add metadata
                scored["uniprot_id"] = row["uniprot_id"]
                scored["protein_file"] = row["filename"]

                # Filter by threshold
                candidates = scored[scored["composite_score"] >= score_threshold]
                if len(candidates) > 0:
                    all_results.append(candidates)

            except Exception as e:
                click.echo(f"\nError processing {row['uniprot_id']}: {e}")
                continue

    # Combine and save results
    if all_results:
        import pandas as pd

        final_results = pd.concat(all_results, ignore_index=True)
        final_results = deterministic_sort_dataframe(
            final_results,
            ["composite_score", "uniprot_id", "pocket_id"],
            ascending=[False, True, True],
        )
        final_results.to_csv(output, index=False)

        click.secho(f"\n✓ Found {len(final_results)} candidate sites", fg="green")
        click.secho(f"Results saved to {output}", fg="green")
    else:
        click.secho("\nNo candidates found above threshold", fg="yellow")


@main.command("md-validate")
@click.argument("candidates_csv", type=click.Path(exists=True))
@click.option("--output-dir", "-o", default="results/md_validation", help="Output directory")
@click.option("--top-n", default=20, show_default=True, help="Top candidates to validate")
def md_validate(candidates_csv, output_dir, top_n):
    """Validate top candidates with OpenMM MD and pocket dynamics analysis."""
    click.echo(f"Running MD validation for top {top_n} candidates from {candidates_csv}...")

    pipeline = OpenMMMDValidationPipeline(output_dir=output_dir)
    report = pipeline.validate_top_candidates(candidates_csv=candidates_csv, top_n=top_n)

    report_path = Path(output_dir) / "md_validation_report.csv"
    stable_count = (report["classification"] == "stably buried").sum()

    click.secho(f"✓ MD validation complete: {report_path}", fg="green")
    click.echo(f"Stable pockets retained: {stable_count}/{len(report)}")


if __name__ == "__main__":
    main()
