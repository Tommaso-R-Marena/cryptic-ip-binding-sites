"""Utility functions for Jupyter notebooks.

Provides robust download functions with fallback mechanisms.
"""

import requests
import gzip
from pathlib import Path
from typing import Optional


def download_alphafold_structure(
    uniprot_id: str,
    output_file: Path,
    version: int = 4,
    timeout: int = 30
) -> Path:
    """Download AlphaFold structure with API+FTP fallback.
    
    Uses the same robust logic as Notebook 01:
    1. Try API metadata to get download URL
    2. Fall back to direct FTP download
    
    Args:
        uniprot_id: UniProt accession (e.g., 'P78563')
        output_file: Where to save the PDB file
        version: AlphaFold model version (default: 4)
        timeout: Request timeout in seconds
        
    Returns:
        Path to downloaded file
        
    Raises:
        requests.HTTPError: If download fails from both API and FTP
    """
    if output_file.exists():
        print(f'✓ Using cached: {output_file}')
        return output_file
    
    print(f'Downloading AlphaFold structure for {uniprot_id}...')
    
    # Try API method first
    api_url = f'https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}'
    print(f'  Fetching metadata from API...')
    
    try:
        response = requests.get(api_url, timeout=10)
        response.raise_for_status()
        data = response.json()
        
        # Extract PDB file URL from API response
        entry = data[0] if isinstance(data, list) and len(data) > 0 else data
        pdb_url = entry.get('pdbUrl') or entry.get('cifUrl')
        
        if not pdb_url:
            # Construct URL from entry ID
            entry_id = entry.get('entryId', f'AF-{uniprot_id}-F1')
            version_num = entry.get('latestVersion', version)
            pdb_url = f'https://alphafold.ebi.ac.uk/files/{entry_id}-model_v{version_num}.pdb'
        
        print(f'  Download URL: {pdb_url}')
        pdb_response = requests.get(pdb_url, timeout=timeout)
        pdb_response.raise_for_status()
        
        output_file.write_bytes(pdb_response.content)
        print(f'✓ Downloaded: {output_file}')
        return output_file
        
    except requests.HTTPError as e:
        print(f'\n✗ API download failed: {e}')
        print(f'Trying FTP fallback...')
        
        # Fallback: Direct FTP download
        ftp_url = f'https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/AF-{uniprot_id}-F1-model_v{version}.pdb.gz'
        print(f'  FTP URL: {ftp_url}')
        
        ftp_response = requests.get(ftp_url, timeout=timeout)
        ftp_response.raise_for_status()
        
        # Decompress and save
        decompressed = gzip.decompress(ftp_response.content)
        output_file.write_bytes(decompressed)
        print(f'✓ Downloaded from FTP: {output_file}')
        return output_file


def download_pdb_structure(
    pdb_id: str,
    output_file: Path,
    timeout: int = 30
) -> Path:
    """Download crystal structure from RCSB PDB.
    
    Args:
        pdb_id: PDB identifier (e.g., '1ZY7')
        output_file: Where to save the PDB file
        timeout: Request timeout in seconds
        
    Returns:
        Path to downloaded file
    """
    if output_file.exists():
        print(f'✓ Using cached: {output_file}')
        return output_file
    
    print(f'Downloading PDB structure {pdb_id}...')
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    
    response = requests.get(url, timeout=timeout)
    response.raise_for_status()
    
    output_file.write_bytes(response.content)
    print(f'✓ Downloaded: {output_file}')
    return output_file


def setup_notebook_environment():
    """Setup environment for notebooks (Colab detection, path setup).
    
    Returns:
        Tuple of (IN_COLAB: bool, data_dir: Path)
    """
    import sys
    import os
    from pathlib import Path
    
    IN_COLAB = 'google.colab' in sys.modules
    
    if IN_COLAB:
        print('Running in Google Colab - installing dependencies...')
        # Dependencies installed via !pip in notebook cells
        
        # Clone repository if needed
        if not Path('cryptic-ip-binding-sites').exists():
            import subprocess
            subprocess.run(
                ['git', 'clone', 'https://github.com/Tommaso-R-Marena/cryptic-ip-binding-sites.git'],
                check=True
            )
            os.chdir('cryptic-ip-binding-sites')
        
        sys.path.insert(0, str(Path.cwd()))
    else:
        # Local or Binder
        sys.path.insert(0, str(Path.cwd().parent))
    
    # Create data directory
    data_dir = Path('notebook_data')
    data_dir.mkdir(exist_ok=True)
    
    print('Setup complete!')
    return IN_COLAB, data_dir


if __name__ == '__main__':
    # Test the functions
    print('Testing notebook utilities...')
    
    from pathlib import Path
    import tempfile
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        # Test AlphaFold download
        print('\n1. Testing AlphaFold download:')
        af_file = tmpdir / 'test_adar2.pdb'
        try:
            download_alphafold_structure('P78563', af_file)
            print(f'   File size: {af_file.stat().st_size / 1024:.1f} KB')
            print('   ✓ AlphaFold download OK')
        except Exception as e:
            print(f'   ✗ Failed: {e}')
        
        # Test PDB download
        print('\n2. Testing PDB download:')
        pdb_file = tmpdir / 'test_1zy7.pdb'
        try:
            download_pdb_structure('1ZY7', pdb_file)
            print(f'   File size: {pdb_file.stat().st_size / 1024:.1f} KB')
            print('   ✓ PDB download OK')
        except Exception as e:
            print(f'   ✗ Failed: {e}')
    
    print('\nTesting complete!')
