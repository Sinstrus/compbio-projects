"""Tests for matplotlib_panel.size_axis_labels().

Runnable as a script (no pytest dependency) — prints PASS/FAIL per case.
Run: python figure_bank/test_matplotlib_panel.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path.home() / "projects"))

from figure_bank.matplotlib_panel import (
    PanelRenderer,
    PanelSizing,
    PlacementResult,
    _check_label_distinctiveness,
    add_headroom,
    auto_place_legend,
    clearance_px_from_subplots,
    headroom_ylim_top,
    max_headroom_label,
    place_corner_text,
    probe_legend_headroom,
    size_axis_labels,
)


# Body-box dimensions per scaffold (from SCAFFOLD_REGISTRY.md)
SCF024_BODY = (429, 298)
SCF016_BODY = (440, 142)
SCF014_BODY = (440, 350)
SCF018_BODY = (290, 170)
SCF019_BODY = (215, 170)
SCF023_BODY = (145, 110)

# figsize matched to body aspect at dpi=200
SCF024_FIGSIZE = (2.86, 1.99)   # 1.44:1
SCF016_FIGSIZE = (2.93, 0.95)   # 3.10:1
SCF014_FIGSIZE = (2.93, 2.33)   # 1.26:1
SCF018_FIGSIZE = (1.93, 1.13)   # 1.71:1
SCF019_FIGSIZE = (1.43, 1.13)   # 1.27:1


_results: list[tuple[str, bool, str]] = []


def _case(name: str, ok: bool, detail: str = "") -> None:
    _results.append((name, ok, detail))
    flag = "PASS" if ok else "FAIL"
    print(f"  [{flag}] {name}" + (f"  — {detail}" if detail else ""))


# ---------------------------------------------------------------------------
# T1 — default fits short labels in roomy body
def t1():
    s = size_axis_labels(
        body_w_px=SCF024_BODY[0], body_h_px=SCF024_BODY[1],
        figsize_inches=SCF024_FIGSIZE, dpi=200,
        labels=["1M", "310K", "98K", "31K"],
        axis="x", axis_priority="y",
    )
    ok = (s.strategy == "default" and s.fontsize_pt == 10.0
          and s.rotation_deg == 0 and s.achieved_padding_ratio >= 0.20)
    _case("T1 default fits short labels (FIG178 doses)", ok,
          f"strategy={s.strategy} pt={s.fontsize_pt} rot={s.rotation_deg} "
          f"pad={s.achieved_padding_ratio:.2f}")
    return s


# ---------------------------------------------------------------------------
# T2 — shrink for medium labels in tighter body
def t2():
    s = size_axis_labels(
        body_w_px=SCF019_BODY[0], body_h_px=SCF019_BODY[1],
        figsize_inches=SCF019_FIGSIZE, dpi=200,
        labels=["VHH3-Mos", "VHH3-Tc", "VHH3-VP2n", "VHH3-IR",
                "VHH27-A", "VHH27-B", "VHH27-C", "WT"],
        axis="x", axis_priority="y",
    )
    ok = s.strategy in ("shrunk", "rotated", "failed")
    _case("T2 medium labels in SCF019-sized body", ok,
          f"strategy={s.strategy} pt={s.fontsize_pt} rot={s.rotation_deg}")
    return s


# ---------------------------------------------------------------------------
# T3 — rotation needed for long labels
def t3():
    s = size_axis_labels(
        body_w_px=SCF018_BODY[0], body_h_px=SCF018_BODY[1],
        figsize_inches=SCF018_FIGSIZE, dpi=200,
        labels=[f"Condition_long_name_{i:02d}" for i in range(10)],
        axis="x", axis_priority="y",
    )
    ok = s.strategy in ("rotated", "shrunk", "failed")
    _case("T3 long labels in SCF018 body force rotation/shrink", ok,
          f"strategy={s.strategy} pt={s.fontsize_pt} rot={s.rotation_deg}")
    return s


# ---------------------------------------------------------------------------
# T4 — horizontal flip recommended when axis_priority allows
def t4():
    s = size_axis_labels(
        body_w_px=SCF018_BODY[0], body_h_px=SCF018_BODY[1],
        figsize_inches=SCF018_FIGSIZE, dpi=200,
        labels=[f"Long_condition_label_xx_{i:02d}" for i in range(10)],
        axis="x", axis_priority="x",
        allow_horizontal_bar_flip=True,
    )
    if s.strategy in ("default", "shrunk", "rotated"):
        ok = s.recommend_horizontal_bars is False
        detail = f"sized OK without flip ({s.strategy} pt={s.fontsize_pt} rot={s.rotation_deg})"
    else:
        ok = s.recommend_horizontal_bars is True
        detail = f"flip recommended (strategy={s.strategy})"
    _case("T4 axis_priority=x enables flip when needed", ok, detail)
    return s


# ---------------------------------------------------------------------------
# T5 — flip blocked when axis_priority=y
def t5():
    s = size_axis_labels(
        body_w_px=SCF018_BODY[0], body_h_px=SCF018_BODY[1],
        figsize_inches=SCF018_FIGSIZE, dpi=200,
        labels=[f"Long_condition_label_xx_{i:02d}" for i in range(10)],
        axis="x", axis_priority="y",
        allow_horizontal_bar_flip=True,
    )
    ok = s.recommend_horizontal_bars is False
    _case("T5 axis_priority=y blocks horizontal flip", ok,
          f"strategy={s.strategy} flip={s.recommend_horizontal_bars}")
    return s


# ---------------------------------------------------------------------------
# T6 — prefix-distinctiveness warning on shared-prefix labels
def t6():
    s = size_axis_labels(
        body_w_px=SCF014_BODY[0], body_h_px=SCF014_BODY[1],
        figsize_inches=SCF014_FIGSIZE, dpi=200,
        labels=["AGGTAAGTAGT", "AGGTAAGTAGC", "AGGTAAGTAAT"],
    )
    ok = any("prefix distinctiveness" in w for w in s.warnings)
    _case("T6 prefix-distinctiveness warning fires (FIG177-style labels)", ok,
          f"warnings={s.warnings}")
    return s


# ---------------------------------------------------------------------------
# T7 — single label
def t7():
    s = size_axis_labels(
        body_w_px=SCF024_BODY[0], body_h_px=SCF024_BODY[1],
        figsize_inches=SCF024_FIGSIZE, dpi=200,
        labels=["WT"],
    )
    ok = s.strategy == "default"
    _case("T7 single label is no-op default", ok,
          f"strategy={s.strategy} pt={s.fontsize_pt}")
    return s


# ---------------------------------------------------------------------------
# T8 — empty labels
def t8():
    s = size_axis_labels(
        body_w_px=SCF024_BODY[0], body_h_px=SCF024_BODY[1],
        figsize_inches=SCF024_FIGSIZE, dpi=200,
        labels=[],
    )
    ok = s.strategy == "default" and s.fontsize_pt == 10.0
    _case("T8 empty labels is no-op", ok,
          f"strategy={s.strategy} pt={s.fontsize_pt}")
    return s


# ---------------------------------------------------------------------------
# T9 — reproducibility
def t9():
    args = dict(
        body_w_px=SCF024_BODY[0], body_h_px=SCF024_BODY[1],
        figsize_inches=SCF024_FIGSIZE, dpi=200,
        labels=["1M", "310K", "98K", "31K"],
    )
    a = size_axis_labels(**args)
    b = size_axis_labels(**args)
    ok = a == b
    _case("T9 same inputs produce identical PanelSizing", ok,
          f"a==b: {ok}")
    return a


# ---------------------------------------------------------------------------
# T10 — reproduces FIG178 hand-tuned values (VOY doses)
def t10():
    s = size_axis_labels(
        body_w_px=SCF024_BODY[0], body_h_px=SCF024_BODY[1],
        figsize_inches=SCF024_FIGSIZE, dpi=200,
        labels=["1M", "310K", "98K", "31K"],
    )
    ok = s.fontsize_pt == 10.0 and s.rotation_deg == 0
    _case("T10 reproduces FIG178 hand-tune (VOY-1631 dose labels)", ok,
          f"hand-tune: pt=10 rot=0; sizer: pt={s.fontsize_pt} rot={s.rotation_deg}")
    return s


# ---------------------------------------------------------------------------
# T10b — reproduces FIG178 fraction labels
def t10b():
    s = size_axis_labels(
        body_w_px=SCF024_BODY[0], body_h_px=SCF024_BODY[1],
        figsize_inches=SCF024_FIGSIZE, dpi=200,
        labels=["F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8"],
    )
    ok = s.fontsize_pt == 10.0 and s.rotation_deg == 0
    _case("T10b reproduces FIG178 hand-tune (F1-F8 fraction labels)", ok,
          f"hand-tune: pt=10 rot=0; sizer: pt={s.fontsize_pt} rot={s.rotation_deg}")
    return s


# ---------------------------------------------------------------------------
# T11 — FIG177 donor labels (sizer's choice for tier+donor combo)
def t11():
    labels = [f"T{t}\nAGGTAAGT" for t in ("3", "4", "5", "6", "7", "7")]
    s = size_axis_labels(
        body_w_px=SCF014_BODY[0], body_h_px=SCF014_BODY[1],
        figsize_inches=SCF014_FIGSIZE, dpi=200,
        labels=labels,
    )
    ok = s.strategy in ("default", "shrunk", "rotated")
    _case("T11 FIG177 donor labels (sizer determines)", ok,
          f"strategy={s.strategy} pt={s.fontsize_pt} rot={s.rotation_deg} "
          f"warnings={list(s.warnings)}")
    return s


# ---------------------------------------------------------------------------
# T12 — min_padding=0 should always succeed at some pt
def t12():
    s = size_axis_labels(
        body_w_px=SCF024_BODY[0], body_h_px=SCF024_BODY[1],
        figsize_inches=SCF024_FIGSIZE, dpi=200,
        labels=["1M", "310K", "98K", "31K"],
        min_padding=0.0,
    )
    ok = s.strategy != "failed"
    _case("T12 min_padding=0 succeeds for fittable labels", ok,
          f"strategy={s.strategy} pt={s.fontsize_pt}")
    return s


# ---------------------------------------------------------------------------
# Validation: input errors
def t_val_negative_padding():
    try:
        size_axis_labels(
            body_w_px=400, body_h_px=300,
            figsize_inches=(2.0, 1.5), dpi=200,
            labels=["A", "B"], min_padding=-0.1,
        )
        _case("V1 negative min_padding raises ValueError", False,
              "no exception raised")
    except ValueError:
        _case("V1 negative min_padding raises ValueError", True)


def t_val_target_below_min():
    try:
        size_axis_labels(
            body_w_px=400, body_h_px=300,
            figsize_inches=(2.0, 1.5), dpi=200,
            labels=["A", "B"], min_pt=8.0, target_pt=7.0,
        )
        _case("V2 target_pt < min_pt raises ValueError", False,
              "no exception raised")
    except ValueError:
        _case("V2 target_pt < min_pt raises ValueError", True)


# ---------------------------------------------------------------------------
# Phase 2 — PanelRenderer
# ---------------------------------------------------------------------------

def t_panelrenderer_basic():
    """Construct PanelRenderer; fig_ax returns figure at right figsize/dpi."""
    r = PanelRenderer(
        body_w_px=429, body_h_px=298,
        figsize_inches=SCF024_FIGSIZE, dpi=200,
    )
    fig, ax = r.fig_ax()
    ok = (
        abs(fig.get_figwidth() - SCF024_FIGSIZE[0]) < 0.01
        and abs(fig.get_figheight() - SCF024_FIGSIZE[1]) < 0.01
        and fig.get_dpi() == 200
    )
    import matplotlib.pyplot as plt
    plt.close(fig)
    _case("P1 PanelRenderer.fig_ax returns configured figure", ok)


def t_panelrenderer_size_x_labels_default():
    """size_x_labels picks default sizing for FIG178 doses."""
    r = PanelRenderer(
        body_w_px=429, body_h_px=298,
        figsize_inches=SCF024_FIGSIZE, dpi=200,
        margins={"left": 0.16, "right": 0.97, "top": 0.95, "bottom": 0.13},
        axis_priority="y",
    )
    import numpy as np
    fig, ax = r.fig_ax()
    x = np.arange(4)
    ax.bar(x, [493314, 469416, 339036, 217838])
    ax.set_xticks(x)
    sizing = r.size_x_labels(ax, labels=["1M", "310K", "98K", "31K"])
    import matplotlib.pyplot as plt
    plt.close(fig)
    ok = (sizing.strategy == "default" and sizing.fontsize_pt == 10.0
          and sizing.rotation_deg == 0)
    _case("P2 PanelRenderer.size_x_labels picks default for short labels", ok,
          f"strategy={sizing.strategy} pt={sizing.fontsize_pt}")


def t_panelrenderer_collision_legend_vs_text():
    """Inject a synthetic overlap (text annotation on top of legend) and
    confirm check_collisions() flags it."""
    import numpy as np
    r = PanelRenderer(
        body_w_px=429, body_h_px=298,
        figsize_inches=SCF024_FIGSIZE, dpi=200,
    )
    fig, ax = r.fig_ax()
    x = np.arange(4)
    ax.bar(x, [1, 2, 3, 4], label="data")
    legend = ax.legend(loc='upper left')
    # Force a text annotation right where the legend sits
    fig.canvas.draw()
    leg_bb = legend.get_window_extent()
    inv = ax.transData.inverted()
    # Place text at legend's data-coord center
    cx_disp = (leg_bb.x0 + leg_bb.x1) / 2
    cy_disp = (leg_bb.y0 + leg_bb.y1) / 2
    cx_data, cy_data = inv.transform((cx_disp, cy_disp))
    ax.text(cx_data, cy_data, "OVERLAP", fontsize=20, fontweight='bold')
    collisions = r.check_collisions(ax)
    import matplotlib.pyplot as plt
    plt.close(fig)
    ok = any("OVERLAP" in c and "legend" in c for c in collisions)
    _case("P3 check_collisions detects synthetic legend/text overlap", ok,
          f"collisions={collisions}")


def t_panelrenderer_collision_clean():
    """Plain bar chart with sized labels should report no collisions."""
    import numpy as np
    r = PanelRenderer(
        body_w_px=429, body_h_px=298,
        figsize_inches=SCF024_FIGSIZE, dpi=200,
        margins={"left": 0.16, "right": 0.97, "top": 0.95, "bottom": 0.13},
    )
    fig, ax = r.fig_ax()
    x = np.arange(4)
    ax.bar(x, [493314, 469416, 339036, 217838])
    ax.set_xticks(x)
    r.size_x_labels(ax, labels=["1M", "310K", "98K", "31K"])
    fig.subplots_adjust(**r.margins)
    collisions = r.check_collisions(ax)
    import matplotlib.pyplot as plt
    plt.close(fig)
    ok = len(collisions) == 0
    _case("P4 check_collisions clean for sized bar chart", ok,
          f"collisions={collisions}")


def t_panelrenderer_save_data_uri():
    """save_as_data_uri returns a valid base64 PNG."""
    import numpy as np
    r = PanelRenderer(
        body_w_px=429, body_h_px=298,
        figsize_inches=SCF024_FIGSIZE, dpi=200,
    )
    fig, ax = r.fig_ax()
    x = np.arange(4)
    ax.bar(x, [1, 2, 3, 4])
    ax.set_xticks(x)
    r.size_x_labels(ax, labels=["A", "B", "C", "D"])
    uri = r.save_as_data_uri(fig)
    ok = (
        uri.startswith("data:image/png;base64,")
        and len(uri) > 100  # non-trivial PNG
    )
    _case("P5 save_as_data_uri returns valid base64 PNG", ok,
          f"uri_length={len(uri)}")


def t_dist_helper():
    ok = (
        _check_label_distinctiveness(["WT", "M1", "M2"]) is None
        and _check_label_distinctiveness(["AGGTAAGTAGT", "AGGTAAGTAGC", "AGGTAAGTAAT"]) is not None
        and _check_label_distinctiveness([]) is None
        and _check_label_distinctiveness(["A"]) is None
    )
    _case("V3 _check_label_distinctiveness helper", ok)


# ---------------------------------------------------------------------------
# T13 — tight clearance disqualifies rotation 90 (FIG177-style bug)
def t13():
    """11-mer donor labels with realistic SCF014 bottom margin: rot=90
    would clip into the bottom of the panel. Sizer should pick rot=30
    or shrink at rot=0 instead."""
    labels = ["AGGTAAGTAGT", "AGGTAAGTAGC", "AGGTAAGTAAT",
              "AGGTAAGTACT", "AGGTAAGTACC", "AGGTAAGTGAA"]
    figsize = SCF014_FIGSIZE
    # Realistic bottom margin (subplots_bottom=0.20 typical for rotated labels)
    clearance = clearance_px_from_subplots(
        figsize_inches=figsize, dpi=200,
        axis="x", subplots_bottom=0.20, axis_label_pt=11,
    )
    s = size_axis_labels(
        body_w_px=SCF014_BODY[0], body_h_px=SCF014_BODY[1],
        figsize_inches=figsize, dpi=200,
        labels=labels,
        axis="x", axis_priority="y",
        tick_label_clearance_px=clearance,
    )
    # 11-mer at pt=10 rot=90 needs ~90px vertical clearance; we have ~70 → must reject
    ok = not (s.rotation_deg == 90 and s.fontsize_pt >= 9)
    _case("T13 clearance gate prevents rot=90 clipping (FIG177 11-mers)", ok,
          f"clearance={clearance:.0f}px strategy={s.strategy} "
          f"pt={s.fontsize_pt} rot={s.rotation_deg}")
    return s


# ---------------------------------------------------------------------------
# T14 — generous clearance allows rotation 90
def t14():
    """Same labels as T13 but with a tall bottom margin → rot=90 OK."""
    labels = ["AGGTAAGTAGT", "AGGTAAGTAGC", "AGGTAAGTAAT",
              "AGGTAAGTACT", "AGGTAAGTACC", "AGGTAAGTGAA"]
    clearance = clearance_px_from_subplots(
        figsize_inches=SCF014_FIGSIZE, dpi=200,
        axis="x", subplots_bottom=0.40, axis_label_pt=11,
    )
    s = size_axis_labels(
        body_w_px=SCF014_BODY[0], body_h_px=SCF014_BODY[1],
        figsize_inches=SCF014_FIGSIZE, dpi=200,
        labels=labels, axis="x", axis_priority="y",
        tick_label_clearance_px=clearance,
    )
    ok = s.strategy != "failed"
    _case("T14 clearance generous → rotation 90 fits", ok,
          f"clearance={clearance:.0f}px strategy={s.strategy} "
          f"pt={s.fontsize_pt} rot={s.rotation_deg}")
    return s


# ---------------------------------------------------------------------------
# T15 — clearance helper math
def t15():
    """clearance_px_from_subplots default (no in-panel axis label):
    figsize_h=2 in, dpi=200, bottom=0.13 → 2 × 200 × 0.13 = 52 px clearance."""
    c = clearance_px_from_subplots(
        figsize_inches=(2.86, 2.0), dpi=200,
        axis="x", subplots_bottom=0.13,
    )
    expected = 2.0 * 200 * 0.13
    ok = abs(c - expected) < 0.5
    _case("T15 clearance_px_from_subplots default math", ok,
          f"got={c:.1f}px expected={expected:.1f}px (no in-panel axis label)")

    # With in-panel axis label, reserve space
    c2 = clearance_px_from_subplots(
        figsize_inches=(2.86, 2.0), dpi=200,
        axis="x", subplots_bottom=0.13, axis_label_pt=11,
    )
    expected2 = 2.0 * 200 * 0.13 - 11 * 200 / 72 * 1.5
    ok2 = abs(c2 - expected2) < 0.5
    _case("T15b clearance with in-panel axis label", ok2,
          f"got={c2:.1f}px expected={expected2:.1f}px")


# ---------------------------------------------------------------------------
# T16 — when budget OK at rot=0 but clearance also OK, default still wins
def t16():
    """Adding clearance constraint shouldn't break the easy default path."""
    s = size_axis_labels(
        body_w_px=SCF024_BODY[0], body_h_px=SCF024_BODY[1],
        figsize_inches=SCF024_FIGSIZE, dpi=200,
        labels=["1M", "310K", "98K", "31K"],
        tick_label_clearance_px=100.0,  # generous
    )
    ok = s.strategy == "default" and s.rotation_deg == 0
    _case("T16 clearance constraint doesn't disturb default path", ok,
          f"strategy={s.strategy} pt={s.fontsize_pt}")


# ---------------------------------------------------------------------------
# T17 — check_collisions(include_data=True) flags legend over Line2D
def t17():
    """Place a legend in the upper-left corner of an axes with a line that
    rises into that corner — the data should be flagged."""
    import matplotlib.pyplot as plt
    import numpy as np
    fig, ax = plt.subplots(figsize=(3.0, 2.0), dpi=200)
    x = np.array([1, 2, 3, 4])
    y = np.array([100, 80, 60, 40])  # descending — high values at upper left
    ax.plot(x, y, label="series A")
    ax.plot(x, y * 0.8, label="series B")
    ax.set_xlim(0.5, 4.5)
    ax.set_ylim(0, 110)
    ax.legend(loc="upper left")
    r = PanelRenderer(
        body_w_px=400, body_h_px=270,
        figsize_inches=(3.0, 2.0), dpi=200,
    )
    collisions = r.check_collisions(ax, include_data=True)
    plt.close(fig)
    ok = any("legend" in c and "Line2D" in c for c in collisions)
    _case("T17 check_collisions detects legend-over-line", ok,
          f"collisions={collisions}")


# ---------------------------------------------------------------------------
# T18 — auto_place_legend picks a clean compass loc
def t18():
    """Line rises to upper-left → upper left fails, but upper right or
    lower right should be clean."""
    import matplotlib.pyplot as plt
    import numpy as np
    fig, ax = plt.subplots(figsize=(3.0, 2.0), dpi=200)
    x = np.array([1, 2, 3, 4])
    ax.plot(x, [100, 80, 60, 40], label="A")
    ax.plot(x, [80, 64, 48, 32], label="B")
    ax.set_xlim(0.5, 4.5)
    ax.set_ylim(0, 110)
    result = auto_place_legend(ax, fontsize=8)
    plt.close(fig)
    ok = isinstance(result, PlacementResult) and not result.overlapped
    _case("T18 auto_place_legend finds clean loc on descending line", ok,
          f"loc={result.loc} headroom={result.headroom_applied!r} "
          f"overlapped={result.overlapped}")


# ---------------------------------------------------------------------------
# T19 — place_corner_text falls through corners
def t19():
    """Put a fat bar in the upper-right corner so prefer='upper right'
    fails. The helper should fall through to another anchor without
    needing headroom."""
    import matplotlib.pyplot as plt
    import numpy as np
    fig, ax = plt.subplots(figsize=(3.0, 2.0), dpi=200)
    x = np.arange(4)
    ax.bar(x, [10, 20, 30, 95], color="tab:blue")  # tall bar at right
    ax.set_xlim(-0.5, 3.5)
    ax.set_ylim(0, 100)
    result = place_corner_text(ax, "×49", prefer="upper right",
                                fontsize=10, fontweight="bold",
                                bbox=dict(boxstyle="round,pad=0.15",
                                          fc="white", ec="none", alpha=0.7))
    plt.close(fig)
    # Sanity: text must end up SOMEWHERE without overlap, NOT at upper-right
    chosen = result.loc
    ok = (
        not result.overlapped
        and isinstance(chosen, tuple)
        and chosen != (0.96, 0.96)  # must have moved off upper-right
    )
    _case("T19 place_corner_text falls through to clean corner", ok,
          f"loc={chosen} headroom={result.headroom_applied!r} "
          f"overlapped={result.overlapped}")


# ---------------------------------------------------------------------------
# T20 — escalation ladder: corners fail → half-decade → full-decade
def t20():
    """Construct an axes where data fills the entire panel area at every
    perimeter slot. Half-decade headroom should clear it."""
    import matplotlib.pyplot as plt
    import numpy as np
    fig, ax = plt.subplots(figsize=(3.0, 2.0), dpi=200)
    # Eight tall bars span the whole x range, height up to 95% of ylim →
    # every horizontal "upper" slot has a bar in it.
    x = np.arange(8)
    heights = [85, 90, 88, 92, 95, 89, 91, 87]
    ax.bar(x, heights, color="tab:blue", width=0.85)
    ax.set_xlim(-0.5, 7.5)
    ax.set_ylim(0, 100)  # linear so bars take 85-95% of vertical
    # Place a wide, short legend that won't fit anywhere along the top
    ax.plot([], [], label="series 1")
    ax.plot([], [], label="series 2")
    ax.plot([], [], label="series 3")
    result = auto_place_legend(ax, fontsize=8, ncol=3)
    plt.close(fig)
    # Either accepted with linear headroom applied, or has overlap with
    # 'full_decade' tag.  Critical check: SOMETHING happened past step 0.
    ok = (
        result.headroom_applied in ("half_decade", "full_decade")
        or not result.overlapped  # cleared by repositioning at step 0
    )
    _case("T20 escalation ladder reaches headroom step or clears earlier", ok,
          f"loc={result.loc} headroom={result.headroom_applied!r} "
          f"overlapped={result.overlapped} "
          f"artists={result.overlap_artists}")


# ---------------------------------------------------------------------------
# T21 — add_headroom is idempotent and monotonic across log_decades
def t21():
    """Calling add_headroom twice with the same step is a no-op; calling
    with a larger step grows the top relative to the original."""
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax.set_ylim(1, 100)  # original top = 100

    add_headroom(ax, log_decades=0.5)
    top_after_half = ax.get_ylim()[1]

    add_headroom(ax, log_decades=0.5)  # no-op — same step
    top_after_repeat = ax.get_ylim()[1]

    add_headroom(ax, log_decades=1.0)  # escalate
    top_after_full = ax.get_ylim()[1]

    plt.close(fig)
    half_ok = abs(top_after_half - 100 * (10 ** 0.5)) < 0.5
    idempotent_ok = abs(top_after_repeat - top_after_half) < 1e-6
    monotonic_ok = top_after_full > top_after_half + 1
    full_ok = abs(top_after_full - 100 * 10.0) < 0.5
    ok = half_ok and idempotent_ok and monotonic_ok and full_ok
    _case("T21 add_headroom is idempotent and escalates from original", ok,
          f"100→{top_after_half:.1f} (half) →{top_after_full:.1f} (full)")


# T22 — probe_legend_headroom is non-destructive and returns correct label
def t22():
    """probe_legend_headroom restores ylim after probing; returns correct level."""
    import matplotlib.pyplot as plt
    import numpy as np

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax.set_ylim(1e3, 1e6)
    # Fill the entire panel with data from top to bottom so every compass
    # position overlaps at the base ylim, forcing the probe to escalate.
    x = np.arange(8)
    for _ in range(5):
        ax.plot(x, np.logspace(3, 6, 8), marker='o')
    fig.canvas.draw()

    orig_ylim = ax.get_ylim()
    label = probe_legend_headroom(ax, fontsize=7, ncol=3)
    restored_ylim = ax.get_ylim()
    no_legend_after = ax.get_legend() is None

    plt.close(fig)

    # ylim must be restored to the original value
    restored_ok = abs(restored_ylim[0] - orig_ylim[0]) < 0.1 and \
                  abs(restored_ylim[1] - orig_ylim[1]) < 0.1
    # legend must be removed after probing
    # label must be a known headroom label
    label_ok = label in ("", "half_decade", "full_decade", "exhausted")
    # headroom_ylim_top must produce a value >= original top
    new_top = headroom_ylim_top(orig_ylim[1], label, yscale="log")
    top_ok = new_top >= orig_ylim[1]
    # max_headroom_label picks the more demanding of two
    max_ok = max_headroom_label(["", "half_decade"]) == "half_decade"

    ok = restored_ok and no_legend_after and label_ok and top_ok and max_ok
    _case("T22 probe_legend_headroom is non-destructive + two-pass helpers", ok,
          f"label={label!r} restored_ylim={restored_ylim[1]:.0f} new_top={new_top:.2g} "
          f"no_legend={no_legend_after}")


# ---------------------------------------------------------------------------
def main() -> int:
    print("Running matplotlib_panel sizer tests...\n")
    t1(); t2(); t3(); t4(); t5(); t6(); t7(); t8()
    t9(); t10(); t10b(); t11(); t12()
    t13(); t14(); t15(); t16()
    t17(); t18(); t19(); t20(); t21(); t22()
    t_panelrenderer_basic()
    t_panelrenderer_size_x_labels_default()
    t_panelrenderer_collision_legend_vs_text()
    t_panelrenderer_collision_clean()
    t_panelrenderer_save_data_uri()
    t_val_negative_padding()
    t_val_target_below_min()
    t_dist_helper()

    n_pass = sum(1 for _, ok, _ in _results if ok)
    n_total = len(_results)
    print(f"\n{n_pass}/{n_total} passed")
    return 0 if n_pass == n_total else 1


if __name__ == "__main__":
    sys.exit(main())
