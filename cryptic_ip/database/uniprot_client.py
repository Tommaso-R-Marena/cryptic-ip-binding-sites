"""UniProt API client for protein annotations."""

import logging
import requests
from typing import Dict, List, Optional
import time

logger = logging.getLogger(__name__)


class UniProtClient:
    """Client for UniProt REST API.
    
    Uses UniProt REST API:
    https://rest.uniprot.org
    """
    
    BASE_URL = "https://rest.uniprot.org"
    
    def __init__(self):
        """Initialize UniProt client."""
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'CrypticIP/1.0 (https://github.com/Tommaso-R-Marena/cryptic-ip-binding-sites)'
        })
    
    def get_protein_info(self, uniprot_id: str) -> Dict:
        """Fetch protein information.
        
        Args:
            uniprot_id: UniProt accession
            
        Returns:
            Dictionary with protein annotations
        """
        url = f"{self.BASE_URL}/uniprotkb/{uniprot_id}.json"
        
        try:
            response = self.session.get(url, timeout=10)
            response.raise_for_status()
            data = response.json()
            
            # Extract key information
            return {
                'accession': data['primaryAccession'],
                'gene_name': data.get('genes', [{}])[0].get('geneName', {}).get('value', ''),
                'protein_name': data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', ''),
                'organism': data.get('organism', {}).get('scientificName', ''),
                'taxonomy_id': data.get('organism', {}).get('taxonId', 0),
                'sequence_length': data.get('sequence', {}).get('length', 0),
                'mass': data.get('sequence', {}).get('molWeight', 0),
                'function': self._extract_function(data),
                'go_terms': self._extract_go_terms(data),
                'subcellular_location': self._extract_location(data),
            }
            
        except requests.HTTPError as e:
            if e.response.status_code == 404:
                raise ValueError(f"UniProt ID {uniprot_id} not found")
            raise
    
    def _extract_function(self, data: Dict) -> str:
        """Extract function annotation."""
        comments = data.get('comments', [])
        for comment in comments:
            if comment.get('commentType') == 'FUNCTION':
                texts = comment.get('texts', [])
                if texts:
                    return texts[0].get('value', '')
        return ''
    
    def _extract_go_terms(self, data: Dict) -> Dict[str, List[str]]:
        """Extract Gene Ontology terms."""
        go_terms = {'biological_process': [], 'molecular_function': [], 'cellular_component': []}
        
        for ref in data.get('uniProtKBCrossReferences', []):
            if ref.get('database') == 'GO':
                props = {p['key']: p['value'] for p in ref.get('properties', [])}
                term = props.get('GoTerm', '')
                category = props.get('GoEvidenceType', '')
                
                if 'P:' in term:
                    go_terms['biological_process'].append(term)
                elif 'F:' in term:
                    go_terms['molecular_function'].append(term)
                elif 'C:' in term:
                    go_terms['cellular_component'].append(term)
        
        return go_terms
    
    def _extract_location(self, data: Dict) -> List[str]:
        """Extract subcellular location."""
        locations = []
        comments = data.get('comments', [])
        for comment in comments:
            if comment.get('commentType') == 'SUBCELLULAR_LOCATION':
                for loc in comment.get('subcellularLocations', []):
                    location = loc.get('location', {}).get('value', '')
                    if location:
                        locations.append(location)
        return locations
    
    def get_sequence(self, uniprot_id: str) -> str:
        """Fetch protein sequence.
        
        Args:
            uniprot_id: UniProt accession
            
        Returns:
            Protein sequence (FASTA format)
        """
        url = f"{self.BASE_URL}/uniprotkb/{uniprot_id}.fasta"
        
        try:
            response = self.session.get(url, timeout=10)
            response.raise_for_status()
            return response.text
        except requests.HTTPError as e:
            if e.response.status_code == 404:
                raise ValueError(f"UniProt ID {uniprot_id} not found")
            raise
    
    def search_by_gene(self, gene_name: str, organism: Optional[str] = None) -> List[str]:
        """Search for proteins by gene name.
        
        Args:
            gene_name: Gene name to search
            organism: Optional organism filter
            
        Returns:
            List of UniProt accessions
        """
        query = f"gene:{gene_name}"
        if organism:
            query += f" AND organism_name:{organism}"
        
        url = f"{self.BASE_URL}/uniprotkb/search"
        params = {
            'query': query,
            'format': 'json',
            'size': 25
        }
        
        try:
            response = self.session.get(url, params=params, timeout=10)
            response.raise_for_status()
            data = response.json()
            
            return [entry['primaryAccession'] for entry in data.get('results', [])]
        except requests.HTTPError:
            return []


if __name__ == "__main__":
    # Example usage
    logging.basicConfig(level=logging.INFO)
    
    client = UniProtClient()
    
    # Get ADAR2 information
    info = client.get_protein_info('P78563')
    print(f"Protein: {info['protein_name']}")
    print(f"Gene: {info['gene_name']}")
    print(f"Organism: {info['organism']}")
    print(f"Length: {info['sequence_length']} aa")
    print(f"\nFunction: {info['function'][:200]}...")
    print(f"\nGO Terms (MF): {info['go_terms']['molecular_function'][:3]}")
