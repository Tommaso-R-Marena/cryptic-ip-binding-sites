"""Tests for scoring module."""

import pytest
import numpy as np
from cryptic_ip.analysis.scoring import CompositeScorer


class TestCompositeScorer:
    """Test suite for CompositeScorer class."""

    def test_initialization(self):
        """Test scorer initialization."""
        scorer = CompositeScorer()
        assert scorer is not None

    def test_score_calculation_basic(self):
        """Test basic score calculation."""
        scorer = CompositeScorer()

        pocket_data = {
            'volume': 500.0,
            'depth': 20.0,
            'sasa': 3.0,
            'potential': 6.0,
            'basic_residues': 6
        }

        score = scorer.calculate_composite_score(pocket_data)
        assert score > 0
        assert score <= 1.0

    def test_adar2_positive_control(self):
        """Test that ADAR2-like pocket scores highly."""
        scorer = CompositeScorer()

        # ADAR2 IP6 site characteristics
        adar2_pocket = {
            'volume': 650.0,  # Appropriate for IP6
            'depth': 18.0,    # Deep buried pocket
            'sasa': 2.0,      # Very low accessibility
            'potential': 7.5, # Strong positive potential
            'basic_residues': 6  # Multiple coordinating residues
        }

        score = scorer.calculate_composite_score(adar2_pocket)
        assert score > 0.7, f"ADAR2-like pocket scored too low: {score}"

    def test_ph_domain_negative_control(self):
        """Test that PH domain-like pocket scores lowly."""
        scorer = CompositeScorer()

        # Surface IP-binding PH domain characteristics
        ph_domain_pocket = {
            'volume': 400.0,
            'depth': 6.0,     # Shallow surface pocket
            'sasa': 65.0,     # High accessibility
            'potential': 3.0, # Moderate positive potential
            'basic_residues': 3
        }

        score = scorer.calculate_composite_score(ph_domain_pocket)
        assert score < 0.4, f"PH domain-like pocket scored too high: {score}"

    def test_score_separation(self):
        """Test that buried and surface sites are clearly separated."""
        scorer = CompositeScorer()

        buried_scores = []
        surface_scores = []

        # Generate multiple buried pocket examples
        for i in range(10):
            buried = {
                'volume': 500 + i * 20,
                'depth': 15 + i,
                'sasa': 1.0 + i * 0.3,
                'potential': 6.0 + i * 0.2,
                'basic_residues': 5 + i % 2
            }
            buried_scores.append(scorer.calculate_composite_score(buried))

        # Generate multiple surface pocket examples
        for i in range(10):
            surface = {
                'volume': 400 + i * 15,
                'depth': 5 + i * 0.5,
                'sasa': 50.0 + i * 3,
                'potential': 3.0 + i * 0.1,
                'basic_residues': 2 + i % 2
            }
            surface_scores.append(scorer.calculate_composite_score(surface))

        # Check for clear separation
        min_buried = min(buried_scores)
        max_surface = max(surface_scores)

        assert min_buried > max_surface, "No clear separation between buried and surface sites"

    def test_parameter_ranges(self):
        """Test score behavior with extreme parameter values."""
        scorer = CompositeScorer()

        # Test with minimum values
        min_pocket = {
            'volume': 100.0,
            'depth': 0.0,
            'sasa': 0.0,
            'potential': 0.0,
            'basic_residues': 0
        }
        min_score = scorer.calculate_composite_score(min_pocket)
        assert 0 <= min_score <= 1.0

        # Test with maximum reasonable values
        max_pocket = {
            'volume': 1000.0,
            'depth': 30.0,
            'sasa': 0.1,
            'potential': 10.0,
            'basic_residues': 10
        }
        max_score = scorer.calculate_composite_score(max_pocket)
        assert 0 <= max_score <= 1.0
        assert max_score > min_score

    def test_weight_customization(self):
        """Test custom weight settings."""
        custom_weights = {
            'sasa': 0.4,
            'depth': 0.2,
            'potential': 0.2,
            'basic_residues': 0.1,
            'volume': 0.1
        }

        scorer = CompositeScorer(weights=custom_weights)

        pocket = {
            'volume': 500.0,
            'depth': 15.0,
            'sasa': 3.0,
            'potential': 5.0,
            'basic_residues': 5
        }

        score = scorer.calculate_composite_score(pocket)
        assert 0 <= score <= 1.0
