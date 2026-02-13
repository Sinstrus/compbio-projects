#!/usr/bin/env python3
"""
design_reasoner.py — General-purpose framework for DNA design decisions.

Addresses Plan item 4e (P0): The agent has strong verification but no systematic
framework for *design* decisions. This module provides the missing framework:

    ENUMERATE options → SCORE on metrics → IDENTIFY trade-offs →
    APPLY thresholds → SELECT best → DOCUMENT rationale

This is the generalized version of the reasoning already embedded in
optimize_intron_body_codons.py, made available for all design decisions.

Usage:
    from design_reasoner import DesignReasoner, Metric, Option

    reasoner = DesignReasoner(decision_name="Select cloning enzyme")

    # Define metrics (what matters for this decision)
    reasoner.add_metric(Metric("reliability", weight=1.0, higher_is_better=True,
                               threshold_min=0.5))
    reasoner.add_metric(Metric("dam_risk", weight=0.5, higher_is_better=False,
                               threshold_max=0.3))

    # Add options (what are the valid choices)
    reasoner.add_option(Option("HindIII", scores={"reliability": 1.0, "dam_risk": 0.0}))
    reasoner.add_option(Option("XbaI", scores={"reliability": 0.8, "dam_risk": 0.4}))

    # Decide
    result = reasoner.decide()
    print(result.selected)       # "HindIII"
    print(result.rationale)      # Human-readable explanation
    print(result.ranking)        # All options ranked
    print(result.trade_offs)     # Identified trade-offs
    print(result.eliminated)     # Options eliminated by thresholds
"""

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class Metric:
    """
    A scoring dimension for design decisions.

    Attributes:
        name: Human-readable metric name (e.g., "reliability", "gc_content").
        weight: Relative importance (0.0-1.0). Higher = more important.
        higher_is_better: True if higher scores are preferred.
        threshold_min: Options with scores below this are eliminated.
        threshold_max: Options with scores above this are eliminated.
        unit: Optional unit string for display (e.g., "%", "kcal/mol").
    """
    name: str
    weight: float = 1.0
    higher_is_better: bool = True
    threshold_min: Optional[float] = None
    threshold_max: Optional[float] = None
    unit: str = ""

    def is_within_threshold(self, value: float) -> bool:
        """Check if a value passes this metric's thresholds."""
        if self.threshold_min is not None and value < self.threshold_min:
            return False
        if self.threshold_max is not None and value > self.threshold_max:
            return False
        return True


@dataclass
class Option:
    """
    A candidate design choice.

    Attributes:
        name: Human-readable option name (e.g., "HindIII", "VR-IV insertion").
        scores: Dict of {metric_name: score_value}.
        metadata: Arbitrary additional data about this option.
    """
    name: str
    scores: Dict[str, float] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class TradeOff:
    """A trade-off between two metrics for a given option."""
    option_name: str
    metric_better: str
    metric_worse: str
    better_score: float
    worse_score: float
    description: str


@dataclass
class EliminatedOption:
    """An option eliminated by threshold violation."""
    option_name: str
    metric_name: str
    score: float
    threshold: str  # e.g., "min=0.5" or "max=0.3"
    reason: str


@dataclass
class RankedOption:
    """An option with its composite score and rank."""
    rank: int
    option: Option
    composite_score: float
    metric_contributions: Dict[str, float]  # Metric name → weighted contribution


@dataclass
class DecisionResult:
    """
    Complete result of a design decision.

    Attributes:
        decision_name: What decision was being made.
        selected: The selected option (or None if all eliminated).
        ranking: All viable options ranked by composite score.
        eliminated: Options eliminated by thresholds.
        trade_offs: Identified trade-offs for the selected option.
        rationale: Human-readable explanation of why selected option was chosen.
    """
    decision_name: str
    selected: Optional[Option]
    ranking: List[RankedOption]
    eliminated: List[EliminatedOption]
    trade_offs: List[TradeOff]
    rationale: str

    def __str__(self):
        lines = [f"=== Design Decision: {self.decision_name} ==="]

        if self.selected:
            lines.append(f"\n  SELECTED: {self.selected.name}")
        else:
            lines.append("\n  SELECTED: None (all options eliminated)")

        if self.ranking:
            lines.append("\n  Ranking:")
            for ro in self.ranking:
                contributions = ", ".join(
                    f"{k}={v:.2f}" for k, v in ro.metric_contributions.items()
                )
                lines.append(
                    f"    #{ro.rank}: {ro.option.name} "
                    f"(score={ro.composite_score:.3f}) [{contributions}]"
                )

        if self.eliminated:
            lines.append("\n  Eliminated:")
            for eo in self.eliminated:
                lines.append(f"    {eo.option_name}: {eo.reason}")

        if self.trade_offs:
            lines.append("\n  Trade-offs for selected:")
            for to in self.trade_offs:
                lines.append(f"    {to.description}")

        lines.append(f"\n  Rationale: {self.rationale}")
        return '\n'.join(lines)


# ---------------------------------------------------------------------------
# DesignReasoner — Main Class
# ---------------------------------------------------------------------------

class DesignReasoner:
    """
    General-purpose framework for making DNA design decisions.

    Implements the decision process:
        1. ENUMERATE options
        2. SCORE each option on relevant metrics
        3. IDENTIFY trade-offs
        4. APPLY thresholds (eliminate options outside acceptable ranges)
        5. SELECT the option with the best weighted composite score
        6. DOCUMENT the rationale

    This is the generalized version of the scoring logic in
    optimize_intron_body_codons.py's do_no_harm_check(), now available
    for any design decision.
    """

    def __init__(self, decision_name: str):
        self.decision_name = decision_name
        self.metrics: Dict[str, Metric] = {}
        self.options: List[Option] = []

    def add_metric(self, metric: Metric):
        """Register a scoring metric."""
        self.metrics[metric.name] = metric

    def add_option(self, option: Option):
        """Register a candidate option."""
        self.options.append(option)

    def _eliminate_by_thresholds(self) -> Tuple[List[Option], List[EliminatedOption]]:
        """Apply threshold filters to eliminate invalid options."""
        viable = []
        eliminated = []

        for option in self.options:
            is_viable = True
            for metric_name, metric in self.metrics.items():
                if metric_name not in option.scores:
                    continue
                score = option.scores[metric_name]
                if not metric.is_within_threshold(score):
                    threshold_str = ""
                    if metric.threshold_min is not None and score < metric.threshold_min:
                        threshold_str = f"min={metric.threshold_min}"
                    if metric.threshold_max is not None and score > metric.threshold_max:
                        threshold_str = f"max={metric.threshold_max}"
                    eliminated.append(EliminatedOption(
                        option_name=option.name,
                        metric_name=metric_name,
                        score=score,
                        threshold=threshold_str,
                        reason=f"{option.name} eliminated: {metric_name}={score:.3f} "
                               f"violates threshold ({threshold_str})",
                    ))
                    is_viable = False
                    break  # One violation is enough

            if is_viable:
                viable.append(option)

        return viable, eliminated

    def _compute_composite_score(self, option: Option) -> Tuple[float, Dict[str, float]]:
        """
        Compute weighted composite score for an option.

        Scores are normalized: if higher_is_better, score is used as-is.
        If lower_is_better, score is negated before weighting.
        """
        total_weight = sum(m.weight for m in self.metrics.values())
        if total_weight == 0:
            return 0.0, {}

        composite = 0.0
        contributions = {}

        for metric_name, metric in self.metrics.items():
            if metric_name not in option.scores:
                continue
            raw_score = option.scores[metric_name]
            # Normalize direction: always make "higher = better" for compositing
            directed_score = raw_score if metric.higher_is_better else -raw_score
            weighted = directed_score * (metric.weight / total_weight)
            composite += weighted
            contributions[metric_name] = weighted

        return composite, contributions

    def _identify_trade_offs(self, selected: Option) -> List[TradeOff]:
        """
        Identify metrics where the selected option is NOT the best.

        A trade-off exists when the selected option is the best overall
        but not the best on every individual metric.
        """
        trade_offs = []

        for metric_name, metric in self.metrics.items():
            if metric_name not in selected.scores:
                continue

            selected_score = selected.scores[metric_name]

            # Find the best score for this metric across all options
            best_score = selected_score
            best_option = selected.name
            for option in self.options:
                if metric_name not in option.scores:
                    continue
                score = option.scores[metric_name]
                if metric.higher_is_better:
                    if score > best_score:
                        best_score = score
                        best_option = option.name
                else:
                    if score < best_score:
                        best_score = score
                        best_option = option.name

            if best_option != selected.name:
                trade_offs.append(TradeOff(
                    option_name=selected.name,
                    metric_better="(composite)",
                    metric_worse=metric_name,
                    better_score=0.0,
                    worse_score=selected_score,
                    description=(
                        f"{selected.name} is not best on {metric_name} "
                        f"({selected_score:.3f}{metric.unit} vs "
                        f"{best_option}'s {best_score:.3f}{metric.unit})"
                    ),
                ))

        return trade_offs

    def _build_rationale(
        self,
        selected: Optional[Option],
        ranking: List[RankedOption],
        eliminated: List[EliminatedOption],
        trade_offs: List[TradeOff],
    ) -> str:
        """Generate a human-readable rationale for the decision."""
        if not selected:
            if eliminated:
                reasons = "; ".join(e.reason for e in eliminated)
                return f"No viable options. All eliminated: {reasons}"
            return "No options provided."

        parts = [f"Selected {selected.name}"]

        if ranking and len(ranking) > 1:
            margin = ranking[0].composite_score - ranking[1].composite_score
            parts.append(
                f"ranked #{ranking[0].rank} with composite score "
                f"{ranking[0].composite_score:.3f} "
                f"(margin={margin:.3f} over #{ranking[1].rank} {ranking[1].option.name})"
            )

        if trade_offs:
            tradeoff_names = [t.metric_worse for t in trade_offs]
            parts.append(f"accepting trade-offs on: {', '.join(tradeoff_names)}")

        if eliminated:
            parts.append(f"{len(eliminated)} option(s) eliminated by thresholds")

        return ". ".join(parts) + "."

    def decide(self) -> DecisionResult:
        """
        Execute the full decision process and return a DecisionResult.

        Steps:
        1. Apply thresholds to eliminate invalid options
        2. Compute composite scores for viable options
        3. Rank by composite score
        4. Select top-ranked option
        5. Identify trade-offs
        6. Build rationale
        """
        # Step 1: Eliminate
        viable, eliminated = self._eliminate_by_thresholds()

        # Step 2-3: Score and rank
        scored = []
        for option in viable:
            composite, contributions = self._compute_composite_score(option)
            scored.append((option, composite, contributions))

        scored.sort(key=lambda x: x[1], reverse=True)

        ranking = []
        for i, (option, composite, contributions) in enumerate(scored):
            ranking.append(RankedOption(
                rank=i + 1,
                option=option,
                composite_score=composite,
                metric_contributions=contributions,
            ))

        # Step 4: Select
        selected = ranking[0].option if ranking else None

        # Step 5: Trade-offs
        trade_offs = self._identify_trade_offs(selected) if selected else []

        # Step 6: Rationale
        rationale = self._build_rationale(selected, ranking, eliminated, trade_offs)

        return DecisionResult(
            decision_name=self.decision_name,
            selected=selected,
            ranking=ranking,
            eliminated=eliminated,
            trade_offs=trade_offs,
            rationale=rationale,
        )


# ---------------------------------------------------------------------------
# Convenience: Pre-built Reasoners for Common Decisions
# ---------------------------------------------------------------------------

def enzyme_selection_reasoner(
    candidates: List[Dict[str, Any]],
) -> DecisionResult:
    """
    Pre-built reasoner for restriction enzyme selection.

    Each candidate dict should have:
        name: str — Enzyme name
        reliability: float (0-1) — Reliability rating
        dam_risk: float (0-1) — Dam sensitivity risk (0=none, 1=fully blocked)
        site_unique: bool — Whether site is unique in the construct
        star_activity: float (0-1) — Star activity risk

    Returns DecisionResult.
    """
    reasoner = DesignReasoner(decision_name="Select restriction enzyme")

    reasoner.add_metric(Metric(
        name="reliability", weight=1.0, higher_is_better=True,
        threshold_min=0.3,
    ))
    reasoner.add_metric(Metric(
        name="dam_risk", weight=0.8, higher_is_better=False,
        threshold_max=0.7,
    ))
    reasoner.add_metric(Metric(
        name="star_activity", weight=0.3, higher_is_better=False,
    ))

    for c in candidates:
        if not c.get("site_unique", True):
            continue  # Pre-filter: site must be unique
        reasoner.add_option(Option(
            name=c["name"],
            scores={
                "reliability": c.get("reliability", 0.5),
                "dam_risk": c.get("dam_risk", 0.0),
                "star_activity": c.get("star_activity", 0.0),
            },
            metadata=c,
        ))

    return reasoner.decide()


def linker_selection_reasoner(
    candidates: List[Dict[str, Any]],
) -> DecisionResult:
    """
    Pre-built reasoner for linker design selection.

    Each candidate dict should have:
        name: str — Linker name/description
        flexibility: float (0-1) — Conformational flexibility
        synthesis_complexity: float (0-1) — Synthesis difficulty (0=easy, 1=hard)
        length_bp: int — Length in base pairs
        gc_content: float (0-1) — GC fraction
        has_repeats: bool — Whether linker contains problematic repeats

    Returns DecisionResult.
    """
    reasoner = DesignReasoner(decision_name="Select linker design")

    reasoner.add_metric(Metric(
        name="flexibility", weight=1.0, higher_is_better=True,
    ))
    reasoner.add_metric(Metric(
        name="synthesis_complexity", weight=0.8, higher_is_better=False,
        threshold_max=0.8,
    ))
    reasoner.add_metric(Metric(
        name="gc_content", weight=0.3, higher_is_better=False,
        threshold_min=0.25, threshold_max=0.65,
    ))

    for c in candidates:
        if c.get("has_repeats", False):
            continue  # BUG-007: skip repetitive linkers
        reasoner.add_option(Option(
            name=c["name"],
            scores={
                "flexibility": c.get("flexibility", 0.5),
                "synthesis_complexity": c.get("synthesis_complexity", 0.3),
                "gc_content": c.get("gc_content", 0.5),
            },
            metadata=c,
        ))

    return reasoner.decide()


def codon_variant_reasoner(
    position_label: str,
    candidates: List[Dict[str, Any]],
) -> DecisionResult:
    """
    Pre-built reasoner for synonymous codon selection (generalized
    from optimize_intron_body_codons.py logic).

    Each candidate dict should have:
        name: str — Codon sequence (e.g., "GTC", "GTG")
        splice_score: float — MaxEntScan score (lower = less splice risk)
        codon_usage: float (0-1) — Human codon usage frequency
        creates_new_gt: bool — Whether this codon creates a new GT dinucleotide
        creates_new_ag: bool — Whether this codon creates a new AG dinucleotide

    Returns DecisionResult.
    """
    reasoner = DesignReasoner(
        decision_name=f"Select codon at {position_label}"
    )

    reasoner.add_metric(Metric(
        name="splice_score", weight=1.0, higher_is_better=False,
        threshold_max=0.0,  # Do-no-harm: MaxEntScan < 0.0 = "Not Detected"
    ))
    reasoner.add_metric(Metric(
        name="codon_usage", weight=0.3, higher_is_better=True,
    ))

    for c in candidates:
        if c.get("creates_new_gt", False) or c.get("creates_new_ag", False):
            continue
        reasoner.add_option(Option(
            name=c["name"],
            scores={
                "splice_score": c.get("splice_score", -10.0),
                "codon_usage": c.get("codon_usage", 0.2),
            },
            metadata=c,
        ))

    return reasoner.decide()


# ---------------------------------------------------------------------------
# CLI Demo
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Demo: Enzyme selection
    print("=" * 60)
    candidates = [
        {"name": "HindIII", "reliability": 1.0, "dam_risk": 0.0,
         "star_activity": 0.2, "site_unique": True},
        {"name": "XbaI", "reliability": 0.8, "dam_risk": 0.4,
         "star_activity": 0.2, "site_unique": True},
        {"name": "BbvCI", "reliability": 0.4, "dam_risk": 0.0,
         "star_activity": 0.1, "site_unique": True},
        {"name": "PstI", "reliability": 0.7, "dam_risk": 0.8,
         "star_activity": 0.3, "site_unique": True},
    ]
    result = enzyme_selection_reasoner(candidates)
    print(result)

    # Demo: Custom multi-metric decision
    print("\n" + "=" * 60)
    reasoner = DesignReasoner(decision_name="Select VR insertion site")
    reasoner.add_metric(Metric("capsid_tolerance", weight=1.0, higher_is_better=True))
    reasoner.add_metric(Metric("surface_exposure", weight=0.8, higher_is_better=True))
    reasoner.add_metric(Metric("assembly_disruption", weight=0.6, higher_is_better=False,
                                threshold_max=0.7))

    reasoner.add_option(Option("VR-IV", scores={
        "capsid_tolerance": 0.9, "surface_exposure": 0.85, "assembly_disruption": 0.2
    }))
    reasoner.add_option(Option("VR-VIII", scores={
        "capsid_tolerance": 0.7, "surface_exposure": 0.95, "assembly_disruption": 0.4
    }))
    reasoner.add_option(Option("VR-I", scores={
        "capsid_tolerance": 0.3, "surface_exposure": 0.6, "assembly_disruption": 0.8
    }))

    result = reasoner.decide()
    print(result)
