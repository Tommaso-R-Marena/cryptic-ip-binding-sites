"""
Download and organize AlphaFold proteome structures.
"""

import urllib.request
import tarfile
from pathlib import Path
from typing import Dict, Optional
from tqdm import tqdm


class ProteomeDownloader:
    """
    Download AlphaFold proteome structures.
    
    Supported organisms:
    - Saccharomyces cerevisiae (yeast): UP000002311
    - Homo sapiens (human): UP000005640
    - Dictyostelium discoideum: UP000002195
    """
    
    PROTEOMES = {
        'yeast': {
            'uniprot_id': 'UP000002311',
            'organism': 'Saccharomyces cerevisiae',
            'taxon': '559292',
            'proteins': 6049,
            'size_gb': 15
        },
        'human': {
            'uniprot_id': 'UP000005640',
            'organism': 'Homo sapiens',
            'taxon': '9606',
            'proteins': 23391,
            'size_gb': 50
        },
        'dictyostelium': {
            'uniprot_id': 'UP000002195',
            'organism': 'Dictyostelium discoideum',
            'taxon': '44689',
            'proteins': 12622,
            'size_gb': 30
        }
    }
    
    ALPHAFOLD_BASE = "https://ftp.ebi.ac.uk/pub/databases/alphafold/latest"
    
    def __init__(self, data_dir: str = "data/structures"):
        """
        Initialize downloader.
        
        Args:
            data_dir: Base directory for proteome structures
        """
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
    
    def download_proteome(self, organism: str, force: bool = False) -> Path:
        """
        Download complete proteome from AlphaFold.
        
        Args:
            organism: Organism key ('yeast', 'human', 'dictyostelium')
            force: Re-download even if exists
            
        Returns:
            Path to proteome directory
        """
        if organism not in self.PROTEOMES:
            raise ValueError(f"Unknown organism: {organism}")
        
        info = self.PROTEOMES[organism]
        proteome_dir = self.data_dir / organism
        
        # Check if already downloaded
        if proteome_dir.exists() and not force:
            print(f"{organism} proteome already exists at {proteome_dir}")
            return proteome_dir
        
        proteome_dir.mkdir(parents=True, exist_ok=True)
        
        # Construct download URL
        uniprot_id = info['uniprot_id']
        taxon = info['taxon']
        filename = f"{uniprot_id}_{taxon}_{organism.upper()}_v4.tar"
        url = f"{self.ALPHAFOLD_BASE}/{filename}"
        
        tar_path = self.data_dir / filename
        
        # Download
        print(f"\nDownloading {organism} proteome...")
        print(f"Size: ~{info['size_gb']} GB, {info['proteins']} proteins")
        print(f"This may take a while...\n")
        
        try:
            with tqdm(unit='B', unit_scale=True, desc=filename) as pbar:
                def report(block_num, block_size, total_size):
                    pbar.total = total_size
                    pbar.update(block_size)
                
                urllib.request.urlretrieve(url, tar_path, reporthook=report)
            
            # Extract
            print(f"\nExtracting to {proteome_dir}...")
            with tarfile.open(tar_path, 'r') as tar:
                tar.extractall(proteome_dir)
            
            # Clean up tar file
            tar_path.unlink()
            
            print(f"\nDownload complete: {proteome_dir}")
            return proteome_dir
            
        except Exception as e:
            print(f"\nDownload failed: {e}")
            # Clean up partial downloads
            if tar_path.exists():
                tar_path.unlink()
            raise
    
    def get_info(self, organism: str) -> Dict:
        """
        Get information about a proteome.
        
        Args:
            organism: Organism key
            
        Returns:
            Proteome information dictionary
        """
        if organism not in self.PROTEOMES:
            raise ValueError(f"Unknown organism: {organism}")
        return self.PROTEOMES[organism]
    
    def list_available(self) -> Dict:
        """
        List all available proteomes.
        
        Returns:
            Dictionary of proteome information
        """
        return self.PROTEOMES
