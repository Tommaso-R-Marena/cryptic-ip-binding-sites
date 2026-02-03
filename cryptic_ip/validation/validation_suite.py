"""
Comprehensive validation suite for pipeline testing.
"""

import pandas as pd
from typing import Dict, List
from pathlib import Path

from .adar2 import validate_adar2
from ..analysis import ProteinAnalyzer


class ValidationSuite:
    """
    Complete validation suite for testing pipeline on known examples.
    
    Positive controls: ADAR2, Pds5B, HDAC1/3 (should detect)
    Negative controls: PH domains (should reject)
    """
    
    POSITIVE_CONTROLS = {
        'ADAR2': {
            'alphafold_id': 'AF-P78563-F1',
            'pdb': '1ZY7',
            'ip_type': 'IP6',
            'description': 'Gold standard - completely buried IP6',
            'expected_score': 0.75
        },
        'Pds5B': {
            'pdb': '5HDT',
            'ip_type': 'IP6',
            'description': 'Cohesin regulator with buried IP6',
            'expected_score': 0.65
        },
        'HDAC1': {
            'pdb': '5ICN',
            'ip_type': 'IP4',
            'description': 'Histone deacetylase with IP4 at interface',
            'expected_score': 0.60
        }
    }
    
    NEGATIVE_CONTROLS = {
        'PLCd1_PH': {
            'pdb': '1MAI',
            'ip_type': 'IP3',
            'description': 'Classic surface-exposed PH domain',
            'expected_score': 0.30
        },
        'Btk_PH': {
            'pdb': '1BTK',
            'ip_type': 'IP4',
            'description': 'Kinase PH domain - membrane targeting',
            'expected_score': 0.35
        }
    }
    
    def __init__(self, data_dir: str = "data/validation"):
        """
        Initialize validation suite.
        
        Args:
            data_dir: Directory for validation structures
        """
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.results = []
    
    def run_positive_controls(self) -> pd.DataFrame:
        """
        Test pipeline on positive controls (known buried sites).
        
        Returns:
            DataFrame with results
        """
        print("\n" + "="*60)
        print("TESTING POSITIVE CONTROLS")
        print("="*60)
        
        results = []
        
        for name, info in self.POSITIVE_CONTROLS.items():
            print(f"\nTesting {name}...")
            try:
                if name == 'ADAR2':
                    result = validate_adar2()
                    results.append({
                        'protein': name,
                        'control_type': 'positive',
                        'ip_type': info['ip_type'],
                        'score': result['top_pocket_score'],
                        'expected': info['expected_score'],
                        'passed': result['validation_passed'],
                        'description': info['description']
                    })
                else:
                    print(f"  (PDB structure download and analysis to be implemented)")
            except Exception as e:
                print(f"  Error: {e}")
        
        return pd.DataFrame(results)
    
    def run_negative_controls(self) -> pd.DataFrame:
        """
        Test pipeline on negative controls (surface binding).
        
        Returns:
            DataFrame with results
        """
        print("\n" + "="*60)
        print("TESTING NEGATIVE CONTROLS")
        print("="*60)
        
        results = []
        
        for name, info in self.NEGATIVE_CONTROLS.items():
            print(f"\nTesting {name}...")
            print(f"  (Implementation pending - should score LOW)")
            # Placeholder for negative control testing
        
        return pd.DataFrame(results)
    
    def run_full_validation(self) -> Dict:
        """
        Run complete validation suite.
        
        Returns:
            Dictionary with validation summary
        """
        positive_results = self.run_positive_controls()
        negative_results = self.run_negative_controls()
        
        summary = {
            'positive_controls': positive_results,
            'negative_controls': negative_results,
            'separation_quality': self._calculate_separation(positive_results, negative_results)
        }
        
        self._print_summary(summary)
        return summary
    
    def _calculate_separation(self, positive: pd.DataFrame, negative: pd.DataFrame) -> Dict:
        """
        Calculate score separation between positive and negative controls.
        
        Args:
            positive: Positive control results
            negative: Negative control results
            
        Returns:
            Separation metrics
        """
        if len(positive) == 0 or len(negative) == 0:
            return {}
        
        pos_mean = positive['score'].mean()
        neg_mean = negative['score'].mean()
        
        return {
            'positive_mean': pos_mean,
            'negative_mean': neg_mean,
            'separation': pos_mean - neg_mean,
            'clear_separation': (pos_mean - neg_mean) > 0.3
        }
    
    def _print_summary(self, summary: Dict):
        """
        Print validation summary.
        
        Args:
            summary: Validation summary dictionary
        """
        print("\n" + "="*60)
        print("VALIDATION SUMMARY")
        print("="*60)
        
        if 'separation_quality' in summary and summary['separation_quality']:
            sep = summary['separation_quality']
            print(f"\nScore Separation:")
            print(f"  Positive controls: {sep.get('positive_mean', 0):.3f}")
            print(f"  Negative controls: {sep.get('negative_mean', 0):.3f}")
            print(f"  Separation: {sep.get('separation', 0):.3f}")
            print(f"  Clear separation: {sep.get('clear_separation', False)}")
        
        print("\n" + "="*60 + "\n")
