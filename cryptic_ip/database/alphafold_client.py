"""AlphaFold Database API client for fetching protein structures."""

import logging
import requests
from pathlib import Path
from typing import Dict, Optional, List
import time

logger = logging.getLogger(__name__)


class AlphaFoldClient:
    """Client for AlphaFold Protein Structure Database.
    
    Uses the official EBI AlphaFold API:
    https://alphafold.ebi.ac.uk/api/docs
    """
    
    BASE_URL = "https://alphafold.ebi.ac.uk"
    API_URL = f"{BASE_URL}/api"
    FILES_URL = f"{BASE_URL}/files"
    
    def __init__(self, cache_dir: Optional[Path] = None):
        """Initialize AlphaFold client.
        
        Args:
            cache_dir: Directory to cache downloaded structures
        """
        self.cache_dir = cache_dir or Path.home() / ".cryptic_ip" / "alphafold_cache"
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'CrypticIP/1.0 (https://github.com/Tommaso-R-Marena/cryptic-ip-binding-sites)'
        })
    
    def fetch_structure(self, uniprot_id: str, version: int = 4) -> Path:
        """Download AlphaFold structure for a UniProt ID.
        
        Args:
            uniprot_id: UniProt accession (e.g., 'P78563' for ADAR2)
            version: AlphaFold model version (default: 4)
            
        Returns:
            Path to downloaded PDB file
            
        Raises:
            ValueError: If UniProt ID not found in AlphaFold
            requests.HTTPError: If download fails
        """
        # Check cache first
        filename = f"AF-{uniprot_id}-F1-model_v{version}.pdb"
        cached_file = self.cache_dir / filename
        
        if cached_file.exists():
            logger.info(f"Using cached structure: {cached_file}")
            return cached_file
        
        # Download structure
        url = f"{self.FILES_URL}/{filename}"
        logger.info(f"Downloading {uniprot_id} from AlphaFold: {url}")
        
        try:
            response = self.session.get(url, timeout=30)
            response.raise_for_status()
            
            # Save to cache
            with open(cached_file, 'wb') as f:
                f.write(response.content)
            
            logger.info(f"Downloaded structure to {cached_file}")
            return cached_file
            
        except requests.HTTPError as e:
            if e.response.status_code == 404:
                raise ValueError(
                    f"UniProt ID {uniprot_id} not found in AlphaFold Database. "
                    f"Check ID or visit {self.BASE_URL}/entry/{uniprot_id}"
                )
            raise
    
    def get_metadata(self, uniprot_id: str) -> Dict:
        """Fetch AlphaFold prediction metadata.
        
        Args:
            uniprot_id: UniProt accession
            
        Returns:
            Dictionary with prediction metadata
        """
        url = f"{self.API_URL}/prediction/{uniprot_id}"
        
        try:
            response = self.session.get(url, timeout=10)
            response.raise_for_status()
            data = response.json()
            
            if not data:
                raise ValueError(f"No AlphaFold prediction for {uniprot_id}")
            
            # Extract first entry (usually only one)
            entry = data[0] if isinstance(data, list) else data
            
            return {
                'uniprot_id': entry['uniprotAccession'],
                'gene': entry.get('gene', ''),
                'organism': entry.get('organismScientificName', ''),
                'sequence_length': entry.get('uniprotSequenceLength', 0),
                'global_metric_value': entry.get('globalMetricValue', 0.0),
                'model_version': entry.get('latestVersion', 4),
                'model_date': entry.get('modelCreatedDate', ''),
            }
            
        except requests.HTTPError as e:
            if e.response.status_code == 404:
                raise ValueError(f"No metadata for {uniprot_id}")
            raise
    
    def fetch_batch(self, uniprot_ids: List[str], delay: float = 0.5) -> Dict[str, Path]:
        """Download multiple structures with rate limiting.
        
        Args:
            uniprot_ids: List of UniProt accessions
            delay: Delay between requests in seconds
            
        Returns:
            Dictionary mapping UniProt ID to file path
        """
        results = {}
        
        for i, uniprot_id in enumerate(uniprot_ids):
            try:
                path = self.fetch_structure(uniprot_id)
                results[uniprot_id] = path
                logger.info(f"Downloaded {i+1}/{len(uniprot_ids)}: {uniprot_id}")
            except Exception as e:
                logger.error(f"Failed to download {uniprot_id}: {e}")
                results[uniprot_id] = None
            
            # Rate limiting
            if i < len(uniprot_ids) - 1:
                time.sleep(delay)
        
        return results
    
    def fetch_proteome(self, proteome_id: str, output_dir: Path) -> int:
        """Download entire AlphaFold proteome.
        
        Note: Large proteomes (>10 GB). Download from FTP is recommended:
        https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/
        
        Args:
            proteome_id: UniProt proteome ID (e.g., 'UP000002311' for yeast)
            output_dir: Directory to save structures
            
        Returns:
            Number of structures downloaded
        """
        logger.warning(
            f"Downloading full proteome {proteome_id}. "
            f"This may take hours. Consider using FTP bulk download:"
            f"\n  wget https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/{proteome_id}_*.tar"
        )
        
        # For now, direct users to FTP
        raise NotImplementedError(
            "Bulk proteome download not implemented via API. "
            "Please use FTP download as documented in README.md"
        )


if __name__ == "__main__":
    # Example usage
    logging.basicConfig(level=logging.INFO)
    
    client = AlphaFoldClient()
    
    # Download ADAR2
    adar2_path = client.fetch_structure('P78563')
    print(f"ADAR2 structure: {adar2_path}")
    
    # Get metadata
    metadata = client.get_metadata('P78563')
    print(f"\nADAR2 metadata:")
    for key, value in metadata.items():
        print(f"  {key}: {value}")
