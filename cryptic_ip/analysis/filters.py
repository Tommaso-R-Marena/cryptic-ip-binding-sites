"""
Filtering and ranking of pocket candidates.
"""

import pandas as pd
from typing import Optional


class CandidateFilter:
    """
    Filter and rank pocket candidates based on confidence scores and criteria.
    """
    
    def __init__(self, min_score: float = 0.60, min_plddt: float = 70.0):
        """
        Initialize filter with thresholds.
        
        Args:
            min_score: Minimum composite score
            min_plddt: Minimum AlphaFold confidence (pLDDT)
        """
        self.min_score = min_score
        self.min_plddt = min_plddt
    
    def filter_by_score(self, results: pd.DataFrame) -> pd.DataFrame:
        """
        Filter pockets by composite score.
        
        Args:
            results: DataFrame with pocket analysis results
            
        Returns:
            Filtered DataFrame
        """
        return results[results['composite_score'] >= self.min_score].copy()
    
    def filter_by_confidence(self,
                           results: pd.DataFrame,
                           plddt_data: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        """
        Filter pockets by AlphaFold structure confidence.
        
        Args:
            results: DataFrame with pocket analysis results
            plddt_data: DataFrame with per-residue pLDDT scores
            
        Returns:
            Filtered DataFrame
        """
        if plddt_data is None:
            print("Warning: No pLDDT data provided, skipping confidence filter")
            return results
        
        # Filter pockets where pocket-lining residues have high confidence
        filtered = []
        for idx, row in results.iterrows():
            # This is a placeholder - in real implementation would
            # check pLDDT for specific pocket residues
            filtered.append(row)
        
        return pd.DataFrame(filtered)
    
    def filter_by_criteria(self,
                          results: pd.DataFrame,
                          min_basic: int = 4,
                          max_sasa: float = 10.0,
                          min_volume: float = 300,
                          max_volume: float = 800) -> pd.DataFrame:
        """
        Apply hard cutoffs for cryptic IP site criteria.
        
        Args:
            results: DataFrame with pocket analysis results
            min_basic: Minimum basic residue count
            max_sasa: Maximum average SASA
            min_volume: Minimum pocket volume
            max_volume: Maximum pocket volume
            
        Returns:
            Filtered DataFrame
        """
        mask = (
            (results['basic_residues'] >= min_basic) &
            (results['sasa'] <= max_sasa) &
            (results['volume'] >= min_volume) &
            (results['volume'] <= max_volume)
        )
        return results[mask].copy()
    
    def rank_candidates(self, results: pd.DataFrame) -> pd.DataFrame:
        """
        Rank candidates by composite score and add classification.
        
        Args:
            results: DataFrame with pocket analysis results
            
        Returns:
            Ranked DataFrame with classification
        """
        ranked = results.sort_values('composite_score', ascending=False).copy()
        ranked['rank'] = range(1, len(ranked) + 1)
        
        # Add classification
        from .scorer import PocketScorer
        scorer = PocketScorer()
        ranked['classification'] = ranked['composite_score'].apply(scorer.classify_site)
        
        return ranked
    
    def get_top_candidates(self, results: pd.DataFrame, n: int = 50) -> pd.DataFrame:
        """
        Get top N candidates for manual inspection.
        
        Args:
            results: DataFrame with pocket analysis results
            n: Number of top candidates to return
            
        Returns:
            Top N candidates
        """
        return self.rank_candidates(results).head(n)
