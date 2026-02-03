"""RCSB Protein Data Bank API client."""

import logging
import requests
from pathlib import Path
from typing import Dict, Optional, List
import json

logger = logging.getLogger(__name__)


class PDBClient:
    """Client for RCSB Protein Data Bank.
    
    Uses RCSB REST API:
    https://data.rcsb.org/
    """
    
    BASE_URL = "https://files.rcsb.org"
    API_URL = "https://data.rcsb.org/rest/v1"
    
    def __init__(self, cache_dir: Optional[Path] = None):
        """Initialize PDB client.
        
        Args:
            cache_dir: Directory to cache downloaded structures
        """
        self.cache_dir = cache_dir or Path.home() / ".cryptic_ip" / "pdb_cache"
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.session = requests.Session()
    
    def fetch_structure(self, pdb_id: str, format: str = 'pdb') -> Path:
        """Download PDB structure.
        
        Args:
            pdb_id: PDB identifier (e.g., '1ZY7')
            format: File format ('pdb', 'cif', 'xml')
            
        Returns:
            Path to downloaded structure file
        """
        pdb_id = pdb_id.upper()
        filename = f"{pdb_id}.{format}"
        cached_file = self.cache_dir / filename
        
        if cached_file.exists():
            logger.info(f"Using cached structure: {cached_file}")
            return cached_file
        
        # Download
        url = f"{self.BASE_URL}/download/{pdb_id}.{format}"
        logger.info(f"Downloading {pdb_id} from RCSB PDB: {url}")
        
        try:
            response = self.session.get(url, timeout=30)
            response.raise_for_status()
            
            with open(cached_file, 'wb') as f:
                f.write(response.content)
            
            logger.info(f"Downloaded structure to {cached_file}")
            return cached_file
            
        except requests.HTTPError as e:
            if e.response.status_code == 404:
                raise ValueError(f"PDB ID {pdb_id} not found")
            raise
    
    def get_entry_info(self, pdb_id: str) -> Dict:
        """Fetch PDB entry metadata.
        
        Args:
            pdb_id: PDB identifier
            
        Returns:
            Dictionary with entry information
        """
        pdb_id = pdb_id.upper()
        url = f"{self.API_URL}/core/entry/{pdb_id}"
        
        try:
            response = self.session.get(url, timeout=10)
            response.raise_for_status()
            data = response.json()
            
            # Extract key information
            return {
                'pdb_id': pdb_id,
                'title': data.get('struct', {}).get('title', ''),
                'resolution': data.get('rcsb_entry_info', {}).get('resolution_combined', [None])[0],
                'method': data.get('exptl', [{}])[0].get('method', ''),
                'deposition_date': data.get('rcsb_accession_info', {}).get('deposit_date', ''),
                'release_date': data.get('rcsb_accession_info', {}).get('initial_release_date', ''),
                'organism': data.get('rcsb_entry_container_identifiers', {}).get('source_organism_scientific_name', [''])[0],
            }
            
        except requests.HTTPError as e:
            if e.response.status_code == 404:
                raise ValueError(f"PDB ID {pdb_id} not found")
            raise
    
    def get_ligands(self, pdb_id: str) -> List[Dict]:
        """Get ligands in PDB structure.
        
        Args:
            pdb_id: PDB identifier
            
        Returns:
            List of ligand dictionaries
        """
        pdb_id = pdb_id.upper()
        url = f"{self.API_URL}/core/chemcomp/{pdb_id}"
        
        try:
            response = self.session.get(url, timeout=10)
            response.raise_for_status()
            data = response.json()
            
            ligands = []
            for comp in data.get('rcsb_chem_comp_container_identifiers', []):
                ligands.append({
                    'id': comp.get('comp_id', ''),
                    'name': comp.get('comp_name', ''),
                    'formula': comp.get('formula', ''),
                })
            
            return ligands
            
        except requests.HTTPError:
            return []
    
    def search_by_uniprot(self, uniprot_id: str) -> List[str]:
        """Find PDB entries for a UniProt ID.
        
        Args:
            uniprot_id: UniProt accession
            
        Returns:
            List of PDB IDs
        """
        url = f"{self.API_URL}/core/uniprot/{uniprot_id}"
        
        try:
            response = self.session.get(url, timeout=10)
            response.raise_for_status()
            data = response.json()
            
            pdb_ids = []
            for entry in data.get('rcsb_uniprot_container_identifiers', []):
                pdb_ids.extend(entry.get('rcsb_id', []))
            
            return list(set(pdb_ids))  # Remove duplicates
            
        except requests.HTTPError:
            return []


if __name__ == "__main__":
    # Example usage
    logging.basicConfig(level=logging.INFO)
    
    client = PDBClient()
    
    # Download ADAR2 crystal structure
    structure = client.fetch_structure('1ZY7')
    print(f"ADAR2 crystal structure: {structure}")
    
    # Get entry info
    info = client.get_entry_info('1ZY7')
    print(f"\nEntry information:")
    for key, value in info.items():
        print(f"  {key}: {value}")
    
    # Find all ADAR2 structures
    adar2_pdbs = client.search_by_uniprot('P78563')
    print(f"\nAll ADAR2 PDB entries: {adar2_pdbs}")
