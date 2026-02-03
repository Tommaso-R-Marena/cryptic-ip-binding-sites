"""
Scoring system for cryptic IP binding site prediction.
"""

import numpy as np
from typing import Dict, Optional


class PocketScorer:
    """
    Score pockets based on cryptic IP binding site criteria.
    
    Criteria:
    1. Pocket depth >15 Å from surface
    2. Low solvent accessibility (SASA <5 ų)
    3. Strong positive electrostatic potential (>5 kT/e)
    4. ≥4 basic residues within 5 Å
    5. Volume 300-800 ų (IP3-IP6 range)
    """
    
    def __init__(self, weights: Optional[Dict[str, float]] = None):
        """
        Initialize scorer with optional custom weights.
        
        Args:
            weights: Dictionary of scoring weights
        """
        # Default weights (sum to 1.0)
        self.weights = weights or {
            'volume': 0.15,
            'depth': 0.25,
            'sasa': 0.30,
            'basic_residues': 0.20,
            'electrostatics': 0.10
        }
    
    def score_volume(self, volume: float) -> float:
        """
        Score pocket volume.
        
        Ideal range: 300-800 ų for IP3-IP6.
        
        Args:
            volume: Pocket volume in cubic Angstroms
            
        Returns:
            Score between 0 and 1
        """
        if 300 <= volume <= 800:
            return 1.0
        elif volume < 300:
            # Penalize small pockets
            return max(0, volume / 300)
        else:
            # Penalize large pockets (but less severely)
            return max(0, 1.0 - (volume - 800) / 1000)
    
    def score_depth(self, depth: float) -> float:
        """
        Score pocket depth.
        
        Deeper pockets are more likely to be buried structural sites.
        Target: >15 Å from surface.
        
        Args:
            depth: Pocket depth metric
            
        Returns:
            Score between 0 and 1
        """
        # Depth scoring: sigmoid function
        # Full credit at 15 Å, half credit at 10 Å
        if depth >= 15:
            return 1.0
        else:
            return 1 / (1 + np.exp(-(depth - 10) / 2))
    
    def score_sasa(self, sasa: float) -> float:
        """
        Score solvent accessibility.
        
        Lower SASA = more buried = better for cryptic sites.
        Target: <5 ų.
        
        Args:
            sasa: Average SASA of pocket residues
            
        Returns:
            Score between 0 and 1
        """
        if sasa <= 5:
            return 1.0
        elif sasa >= 50:
            return 0.0
        else:
            # Linear decay from 5 to 50 ų
            return (50 - sasa) / 45
    
    def score_basic_residues(self, count: int) -> float:
        """
        Score basic residue count.
        
        IP molecules are highly negatively charged and require
        multiple basic residues for coordination.
        Target: ≥4 basic residues.
        
        Args:
            count: Number of Arg/Lys/His residues near pocket
            
        Returns:
            Score between 0 and 1
        """
        if count >= 6:
            return 1.0
        elif count >= 4:
            return 0.8
        elif count >= 2:
            return 0.4
        else:
            return 0.0
    
    def score_electrostatics(self, potential: Optional[float]) -> float:
        """
        Score electrostatic potential.
        
        Positive potential required for IP binding.
        Target: >5 kT/e.
        
        Args:
            potential: Electrostatic potential in kT/e (None if unavailable)
            
        Returns:
            Score between 0 and 1
        """
        if potential is None:
            # Neutral score if electrostatics not calculated
            return 0.5
        
        if potential >= 5:
            return 1.0
        elif potential >= 3:
            return 0.7
        elif potential >= 1:
            return 0.4
        else:
            return 0.0
    
    def calculate_composite_score(self,
                                 volume: float,
                                 depth: float,
                                 sasa: float,
                                 basic_count: int,
                                 potential: Optional[float] = None) -> float:
        """
        Calculate weighted composite score.
        
        Args:
            volume: Pocket volume
            depth: Pocket depth
            sasa: Average SASA
            basic_count: Number of basic residues
            potential: Electrostatic potential (optional)
            
        Returns:
            Composite score between 0 and 1
        """
        scores = {
            'volume': self.score_volume(volume),
            'depth': self.score_depth(depth),
            'sasa': self.score_sasa(sasa),
            'basic_residues': self.score_basic_residues(basic_count),
            'electrostatics': self.score_electrostatics(potential)
        }
        
        composite = sum(scores[k] * self.weights[k] for k in scores)
        return composite
    
    def classify_site(self, score: float) -> str:
        """
        Classify pocket based on composite score.
        
        Args:
            score: Composite score
            
        Returns:
            Classification string
        """
        if score >= 0.75:
            return "High confidence cryptic IP site"
        elif score >= 0.60:
            return "Moderate confidence candidate"
        elif score >= 0.40:
            return "Low confidence - manual inspection recommended"
        else:
            return "Unlikely cryptic IP site"
