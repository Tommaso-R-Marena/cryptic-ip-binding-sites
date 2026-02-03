"""
Command-line interface for cryptic IP site detection.
"""

import click
from pathlib import Path

from .analysis import ProteinAnalyzer
from .validation import validate_adar2, ValidationSuite
from .database import ProteomeDownloader, ProteomeManager


@click.group()
@click.version_option(version='0.1.0')
def main():
    """
    Cryptic IP Binding Site Detection Pipeline
    
    Identify buried inositol phosphate binding sites in protein structures.
    """
    pass


@main.command()
@click.argument('pdb_file', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), help='Output file for results')
@click.option('--score-threshold', '-t', default=0.60, help='Minimum score threshold')
def analyze(pdb_file, output, score_threshold):
    """
    Analyze a single protein structure for cryptic IP binding sites.
    """
    click.echo(f"Analyzing {pdb_file}...\n")
    
    # Create analyzer
    analyzer = ProteinAnalyzer(pdb_file)
    
    # Detect and score pockets
    click.echo("Detecting pockets...")
    pockets = analyzer.detect_pockets()
    click.echo(f"Found {len(pockets)} pockets\n")
    
    click.echo("Scoring pockets...")
    scored = analyzer.score_all_pockets()
    
    # Filter by threshold
    candidates = scored[scored['composite_score'] >= score_threshold]
    
    click.echo(f"\nFound {len(candidates)} candidates above threshold {score_threshold}\n")
    
    if len(candidates) > 0:
        click.echo("Top candidates:")
        click.echo(candidates[['pocket_id', 'composite_score', 'volume', 'sasa', 'basic_residues']].head(10).to_string())
    
    # Save results
    if output:
        scored.to_csv(output, index=False)
        click.echo(f"\nResults saved to {output}")


@main.command()
@click.option('--alphafold/--crystal', default=True, help='Use AlphaFold or crystal structure')
def validate(alphafold):
    """
    Validate pipeline on ADAR2 (gold standard).
    """
    click.echo("Running ADAR2 validation...\n")
    results = validate_adar2(use_alphafold=alphafold)
    
    if results['validation_passed']:
        click.secho("\n✓ Validation PASSED", fg='green', bold=True)
    else:
        click.secho("\n✗ Validation FAILED", fg='red', bold=True)


@main.command()
def validate_suite():
    """
    Run full validation suite (positive and negative controls).
    """
    click.echo("Running complete validation suite...\n")
    suite = ValidationSuite()
    results = suite.run_full_validation()


@main.command()
@click.argument('organism', type=click.Choice(['yeast', 'human', 'dictyostelium']))
@click.option('--data-dir', '-d', default='data/structures', help='Data directory')
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
    
    if click.confirm('Continue with download?'):
        proteome_dir = downloader.download_proteome(organism)
        click.secho(f"\n✓ Download complete: {proteome_dir}", fg='green')


@main.command()
@click.argument('proteome_dir', type=click.Path(exists=True))
@click.option('--output', '-o', default='screening_results.csv', help='Output CSV file')
@click.option('--score-threshold', '-t', default=0.60, help='Minimum score threshold')
@click.option('--max-structures', '-n', type=int, help='Maximum structures to process')
def screen(proteome_dir, output, score_threshold, max_structures):
    """
    Screen entire proteome for cryptic IP binding sites.
    """
    click.echo(f"Screening proteome: {proteome_dir}\n")
    
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
                pdb_path = row['filepath']
                analyzer = ProteinAnalyzer(pdb_path)
                scored = analyzer.score_all_pockets()
                
                # Add metadata
                scored['uniprot_id'] = row['uniprot_id']
                scored['protein_file'] = row['filename']
                
                # Filter by threshold
                candidates = scored[scored['composite_score'] >= score_threshold]
                if len(candidates) > 0:
                    all_results.append(candidates)
                    
            except Exception as e:
                click.echo(f"\nError processing {row['uniprot_id']}: {e}")
                continue
    
    # Combine and save results
    if all_results:
        import pandas as pd
        final_results = pd.concat(all_results, ignore_index=True)
        final_results = final_results.sort_values('composite_score', ascending=False)
        final_results.to_csv(output, index=False)
        
        click.secho(f"\n✓ Found {len(final_results)} candidate sites", fg='green')
        click.secho(f"Results saved to {output}", fg='green')
    else:
        click.secho("\nNo candidates found above threshold", fg='yellow')


if __name__ == '__main__':
    main()
