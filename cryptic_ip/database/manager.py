"""
Manage and organize proteome structure databases.
"""

import pandas as pd
from pathlib import Path
from typing import List, Optional, Dict
from tqdm import tqdm


class ProteomeManager:
    """
    Manage proteome structure databases and metadata.
    """
    
    def __init__(self, proteome_dir: str):
        """
        Initialize manager for a proteome directory.
        
        Args:
            proteome_dir: Path to proteome structures
        """
        self.proteome_dir = Path(proteome_dir)
        if not self.proteome_dir.exists():
            raise ValueError(f"Proteome directory not found: {proteome_dir}")
        
        self.catalog = None
    
    def build_catalog(self, force: bool = False) -> pd.DataFrame:
        """
        Build catalog of all structures in proteome.
        
        Args:
            force: Rebuild even if catalog exists
            
        Returns:
            DataFrame with structure metadata
        """
        catalog_file = self.proteome_dir / "catalog.csv"
        
        if catalog_file.exists() and not force:
            print(f"Loading existing catalog from {catalog_file}")
            self.catalog = pd.read_csv(catalog_file)
            return self.catalog
        
        print(f"Building structure catalog for {self.proteome_dir}...")
        
        # Find all PDB files
        pdb_files = list(self.proteome_dir.glob("**/*.pdb"))
        print(f"Found {len(pdb_files)} structures")
        
        records = []
        for pdb_file in tqdm(pdb_files, desc="Cataloging structures"):
            # Extract UniProt ID from filename
            # Format: AF-{UNIPROT_ID}-F1-model_v4.pdb
            name = pdb_file.stem
            parts = name.split('-')
            
            if len(parts) >= 2:
                uniprot_id = parts[1]
            else:
                uniprot_id = name
            
            record = {
                'uniprot_id': uniprot_id,
                'filename': pdb_file.name,
                'filepath': str(pdb_file),
                'file_size': pdb_file.stat().st_size
            }
            
            # Try to extract pLDDT from file if available
            # (This would require parsing PDB file - placeholder)
            record['mean_plddt'] = None
            
            records.append(record)
        
        self.catalog = pd.DataFrame(records)
        
        # Save catalog
        self.catalog.to_csv(catalog_file, index=False)
        print(f"Catalog saved to {catalog_file}")
        
        return self.catalog
    
    def get_structure_path(self, uniprot_id: str) -> Optional[Path]:
        """
        Get path to structure for a given UniProt ID.
        
        Args:
            uniprot_id: UniProt identifier
            
        Returns:
            Path to PDB file or None if not found
        """
        if self.catalog is None:
            self.build_catalog()
        
        matches = self.catalog[self.catalog['uniprot_id'] == uniprot_id]
        if len(matches) == 0:
            return None
        
        return Path(matches.iloc[0]['filepath'])
    
    def filter_by_confidence(self, min_plddt: float = 70.0) -> pd.DataFrame:
        """
        Filter structures by average pLDDT confidence.
        
        Args:
            min_plddt: Minimum average pLDDT score
            
        Returns:
            Filtered catalog
        """
        if self.catalog is None:
            self.build_catalog()
        
        # Placeholder - requires parsing pLDDT from files
        print("Warning: pLDDT filtering not yet implemented")
        return self.catalog
    
    def get_statistics(self) -> Dict:
        """
        Get statistics about the proteome.
        
        Returns:
            Dictionary with proteome statistics
        """
        if self.catalog is None:
            self.build_catalog()
        
        total_size = self.catalog['file_size'].sum() / (1024**3)  # GB
        
        return {
            'total_structures': len(self.catalog),
            'total_size_gb': total_size,
            'unique_proteins': self.catalog['uniprot_id'].nunique()
        }
