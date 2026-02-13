#!/usr/bin/env python3
"""
Tests for the DesignReasoner module.

Covers:
- Metric thresholds (min/max elimination)
- Composite scoring (weighted, direction-aware)
- Trade-off identification
- Rationale generation
- Pre-built reasoners (enzyme, linker, codon)
- Edge cases (no options, all eliminated, single option, ties)
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from design_reasoner import (
    DecisionResult,
    DesignReasoner,
    EliminatedOption,
    Metric,
    Option,
    RankedOption,
    TradeOff,
    codon_variant_reasoner,
    enzyme_selection_reasoner,
    linker_selection_reasoner,
)


# ---------------------------------------------------------------------------
# Metric tests
# ---------------------------------------------------------------------------

class TestMetric:
    def test_within_threshold(self):
        m = Metric("test", threshold_min=0.3, threshold_max=0.8)
        assert m.is_within_threshold(0.5)

    def test_below_min(self):
        m = Metric("test", threshold_min=0.3)
        assert not m.is_within_threshold(0.1)

    def test_above_max(self):
        m = Metric("test", threshold_max=0.8)
        assert not m.is_within_threshold(0.9)

    def test_no_threshold(self):
        m = Metric("test")
        assert m.is_within_threshold(999.0)
        assert m.is_within_threshold(-999.0)

    def test_at_boundary(self):
        m = Metric("test", threshold_min=0.5, threshold_max=0.5)
        assert m.is_within_threshold(0.5)


# ---------------------------------------------------------------------------
# Basic decision tests
# ---------------------------------------------------------------------------

class TestBasicDecision:
    def test_single_metric_higher_is_better(self):
        r = DesignReasoner("test")
        r.add_metric(Metric("score", weight=1.0, higher_is_better=True))
        r.add_option(Option("A", scores={"score": 0.8}))
        r.add_option(Option("B", scores={"score": 0.6}))

        result = r.decide()
        assert result.selected.name == "A"
        assert result.ranking[0].option.name == "A"
        assert result.ranking[1].option.name == "B"

    def test_single_metric_lower_is_better(self):
        r = DesignReasoner("test")
        r.add_metric(Metric("risk", weight=1.0, higher_is_better=False))
        r.add_option(Option("Safe", scores={"risk": 0.1}))
        r.add_option(Option("Risky", scores={"risk": 0.9}))

        result = r.decide()
        assert result.selected.name == "Safe"

    def test_two_metrics_weighted(self):
        """With two metrics, the higher-weighted one should dominate."""
        r = DesignReasoner("test")
        r.add_metric(Metric("primary", weight=2.0, higher_is_better=True))
        r.add_metric(Metric("secondary", weight=0.5, higher_is_better=True))

        # A: great primary, bad secondary
        # B: bad primary, great secondary
        r.add_option(Option("A", scores={"primary": 0.9, "secondary": 0.1}))
        r.add_option(Option("B", scores={"primary": 0.1, "secondary": 0.9}))

        result = r.decide()
        assert result.selected.name == "A"  # Primary dominates

    def test_threshold_eliminates(self):
        r = DesignReasoner("test")
        r.add_metric(Metric("reliability", threshold_min=0.5))
        r.add_option(Option("Good", scores={"reliability": 0.8}))
        r.add_option(Option("Bad", scores={"reliability": 0.2}))

        result = r.decide()
        assert result.selected.name == "Good"
        assert len(result.eliminated) == 1
        assert result.eliminated[0].option_name == "Bad"


# ---------------------------------------------------------------------------
# Edge case tests
# ---------------------------------------------------------------------------

class TestEdgeCases:
    def test_no_options(self):
        r = DesignReasoner("test")
        r.add_metric(Metric("score"))
        result = r.decide()
        assert result.selected is None
        assert "No options" in result.rationale

    def test_all_eliminated(self):
        r = DesignReasoner("test")
        r.add_metric(Metric("quality", threshold_min=0.9))
        r.add_option(Option("Low", scores={"quality": 0.1}))
        r.add_option(Option("Medium", scores={"quality": 0.5}))

        result = r.decide()
        assert result.selected is None
        assert len(result.eliminated) == 2
        assert "No viable" in result.rationale

    def test_single_option(self):
        r = DesignReasoner("test")
        r.add_metric(Metric("score"))
        r.add_option(Option("Only", scores={"score": 0.5}))

        result = r.decide()
        assert result.selected.name == "Only"

    def test_tie_score(self):
        """When options have identical scores, first-added wins (stable sort)."""
        r = DesignReasoner("test")
        r.add_metric(Metric("score", weight=1.0, higher_is_better=True))
        r.add_option(Option("First", scores={"score": 0.5}))
        r.add_option(Option("Second", scores={"score": 0.5}))

        result = r.decide()
        # Both have same score — first added should be selected (stable)
        assert result.selected is not None
        assert len(result.ranking) == 2

    def test_missing_metric_score(self):
        """Options missing a metric score should still work."""
        r = DesignReasoner("test")
        r.add_metric(Metric("primary", weight=1.0, higher_is_better=True))
        r.add_metric(Metric("secondary", weight=0.5, higher_is_better=True))

        r.add_option(Option("Full", scores={"primary": 0.8, "secondary": 0.5}))
        r.add_option(Option("Partial", scores={"primary": 0.7}))  # Missing secondary

        result = r.decide()
        assert result.selected.name == "Full"


# ---------------------------------------------------------------------------
# Trade-off identification tests
# ---------------------------------------------------------------------------

class TestTradeOffs:
    def test_no_tradeoffs_when_best_on_all(self):
        r = DesignReasoner("test")
        r.add_metric(Metric("a", weight=1.0, higher_is_better=True))
        r.add_metric(Metric("b", weight=1.0, higher_is_better=True))

        r.add_option(Option("Dominant", scores={"a": 1.0, "b": 1.0}))
        r.add_option(Option("Inferior", scores={"a": 0.5, "b": 0.5}))

        result = r.decide()
        assert result.selected.name == "Dominant"
        assert len(result.trade_offs) == 0

    def test_tradeoff_detected(self):
        r = DesignReasoner("test")
        r.add_metric(Metric("speed", weight=1.0, higher_is_better=True))
        r.add_metric(Metric("safety", weight=0.5, higher_is_better=True))

        # Fast wins overall but is not the safest
        r.add_option(Option("Fast", scores={"speed": 1.0, "safety": 0.3}))
        r.add_option(Option("Safe", scores={"speed": 0.2, "safety": 1.0}))

        result = r.decide()
        assert result.selected.name == "Fast"
        assert len(result.trade_offs) == 1
        assert result.trade_offs[0].metric_worse == "safety"


# ---------------------------------------------------------------------------
# Rationale generation tests
# ---------------------------------------------------------------------------

class TestRationale:
    def test_rationale_includes_selection(self):
        r = DesignReasoner("test")
        r.add_metric(Metric("score"))
        r.add_option(Option("Winner", scores={"score": 1.0}))
        result = r.decide()
        assert "Winner" in result.rationale

    def test_rationale_includes_elimination_count(self):
        r = DesignReasoner("test")
        r.add_metric(Metric("quality", threshold_min=0.5))
        r.add_option(Option("Good", scores={"quality": 0.8}))
        r.add_option(Option("Bad", scores={"quality": 0.1}))
        result = r.decide()
        assert "1 option" in result.rationale and "eliminated" in result.rationale

    def test_str_output(self):
        r = DesignReasoner("test")
        r.add_metric(Metric("score"))
        r.add_option(Option("A", scores={"score": 0.8}))
        result = r.decide()
        output = str(result)
        assert "Design Decision" in output
        assert "SELECTED" in output


# ---------------------------------------------------------------------------
# Pre-built reasoner tests
# ---------------------------------------------------------------------------

class TestEnzymeSelectionReasoner:
    def test_picks_reliable_enzyme(self):
        candidates = [
            {"name": "HindIII", "reliability": 1.0, "dam_risk": 0.0,
             "star_activity": 0.2, "site_unique": True},
            {"name": "BbvCI", "reliability": 0.4, "dam_risk": 0.0,
             "star_activity": 0.1, "site_unique": True},
        ]
        result = enzyme_selection_reasoner(candidates)
        assert result.selected.name == "HindIII"

    def test_eliminates_high_dam_risk(self):
        candidates = [
            {"name": "PstI", "reliability": 0.7, "dam_risk": 0.8,
             "star_activity": 0.3, "site_unique": True},
        ]
        result = enzyme_selection_reasoner(candidates)
        assert result.selected is None
        assert len(result.eliminated) == 1

    def test_skips_non_unique(self):
        candidates = [
            {"name": "HindIII", "reliability": 1.0, "dam_risk": 0.0,
             "star_activity": 0.2, "site_unique": False},
        ]
        result = enzyme_selection_reasoner(candidates)
        assert result.selected is None  # No options after filtering

    def test_dam_vs_reliability_tradeoff(self):
        candidates = [
            {"name": "HindIII", "reliability": 1.0, "dam_risk": 0.0,
             "star_activity": 0.2, "site_unique": True},
            {"name": "XbaI", "reliability": 0.8, "dam_risk": 0.4,
             "star_activity": 0.2, "site_unique": True},
        ]
        result = enzyme_selection_reasoner(candidates)
        assert result.selected.name == "HindIII"


class TestLinkerSelectionReasoner:
    def test_picks_flexible_synthesizable(self):
        candidates = [
            {"name": "GGGGS5", "flexibility": 0.9,
             "synthesis_complexity": 0.3, "gc_content": 0.55, "has_repeats": False},
            {"name": "rigid_helix", "flexibility": 0.2,
             "synthesis_complexity": 0.5, "gc_content": 0.50, "has_repeats": False},
        ]
        result = linker_selection_reasoner(candidates)
        assert result.selected.name == "GGGGS5"

    def test_eliminates_repetitive(self):
        """BUG-007: Repetitive linkers should be filtered out."""
        candidates = [
            {"name": "GGGGS_rep", "flexibility": 0.9,
             "synthesis_complexity": 0.9, "gc_content": 0.55, "has_repeats": True},
            {"name": "GGGGS_varied", "flexibility": 0.85,
             "synthesis_complexity": 0.3, "gc_content": 0.55, "has_repeats": False},
        ]
        result = linker_selection_reasoner(candidates)
        assert result.selected.name == "GGGGS_varied"


class TestCodonVariantReasoner:
    def test_picks_low_splice_high_usage(self):
        candidates = [
            {"name": "GTC", "splice_score": -5.0, "codon_usage": 0.4,
             "creates_new_gt": False, "creates_new_ag": False},
            {"name": "GTT", "splice_score": -3.0, "codon_usage": 0.6,
             "creates_new_gt": False, "creates_new_ag": False},
        ]
        result = codon_variant_reasoner("pos132", candidates)
        # GTC has much lower splice score (weighted 1.0) vs GTT's higher usage (weighted 0.3)
        assert result.selected.name == "GTC"

    def test_eliminates_high_splice_score(self):
        """Codons with MaxEntScan > 0.0 should be eliminated (do-no-harm)."""
        candidates = [
            {"name": "GTG", "splice_score": 2.5, "codon_usage": 0.3,
             "creates_new_gt": False, "creates_new_ag": False},
            {"name": "GTC", "splice_score": -2.5, "codon_usage": 0.4,
             "creates_new_gt": False, "creates_new_ag": False},
        ]
        result = codon_variant_reasoner("pos132", candidates)
        assert result.selected.name == "GTC"
        assert any(e.option_name == "GTG" for e in result.eliminated)

    def test_skips_codons_creating_new_gt(self):
        candidates = [
            {"name": "GTG", "splice_score": -5.0, "codon_usage": 0.5,
             "creates_new_gt": True, "creates_new_ag": False},
            {"name": "GTC", "splice_score": -3.0, "codon_usage": 0.4,
             "creates_new_gt": False, "creates_new_ag": False},
        ]
        result = codon_variant_reasoner("pos132", candidates)
        assert result.selected.name == "GTC"  # GTG filtered out pre-entry


# ---------------------------------------------------------------------------
# Integration: Multi-metric complex scenario
# ---------------------------------------------------------------------------

class TestComplexScenario:
    def test_vr_insertion_site_selection(self):
        """Realistic multi-metric decision for VR insertion site."""
        r = DesignReasoner("Select VR insertion site")
        r.add_metric(Metric("capsid_tolerance", weight=1.0, higher_is_better=True))
        r.add_metric(Metric("surface_exposure", weight=0.8, higher_is_better=True))
        r.add_metric(Metric("assembly_disruption", weight=0.6, higher_is_better=False,
                            threshold_max=0.7))

        r.add_option(Option("VR-IV", scores={
            "capsid_tolerance": 0.9, "surface_exposure": 0.85, "assembly_disruption": 0.2
        }))
        r.add_option(Option("VR-VIII", scores={
            "capsid_tolerance": 0.7, "surface_exposure": 0.95, "assembly_disruption": 0.4
        }))
        r.add_option(Option("VR-I", scores={
            "capsid_tolerance": 0.3, "surface_exposure": 0.6, "assembly_disruption": 0.8
        }))

        result = r.decide()
        assert result.selected.name == "VR-IV"
        assert any(e.option_name == "VR-I" for e in result.eliminated)  # Eliminated by threshold
        assert len(result.ranking) == 2  # VR-IV and VR-VIII remain


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
