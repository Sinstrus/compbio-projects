"""Preview generator — renders each scaffold's preview() to previews/SCF###-preview.svg.

Two variants per scaffold:
    previews/SCF###-preview.svg          (with slot outlines, for inspection)
    previews/SCF###-preview-clean.svg    (no outlines, what real figures look like)

Run after editing any scaffold:
    python -m figure_bank.scaffolds._previews regenerate
    python -m figure_bank.scaffolds._previews regenerate SCF002
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from typing import Iterable


_SCAFFOLDS = ("SCF001", "SCF002", "SCF003", "SCF005", "SCF006",
              "SCF007", "SCF008", "SCF009", "SCF010", "SCF011", "SCF012",
              "SCF013", "SCF014", "SCF015", "SCF016", "SCF017",
              "SCF018", "SCF019", "SCF020", "SCF021", "SCF022", "SCF023",
              "SCF024", "SCF025")
_MODULES = {
    "SCF001": "figure_bank.scaffolds.SCF001_single_panel",
    "SCF002": "figure_bank.scaffolds.SCF002_two_panel",
    "SCF003": "figure_bank.scaffolds.SCF003_three_panel",
    "SCF005": "figure_bank.scaffolds.SCF005_repcap_three_zone",
    "SCF006": "figure_bank.scaffolds.SCF006_architecture_matrix",
    "SCF007": "figure_bank.scaffolds.SCF007_pulldown_pipeline_5stage",
    "SCF008": "figure_bank.scaffolds.SCF008_card_workflow_3step",
    "SCF009": "figure_bank.scaffolds.SCF009_card_workflow_4step",
    "SCF010": "figure_bank.scaffolds.SCF010_card_workflow_5step",
    "SCF011": "figure_bank.scaffolds.SCF011_card_workflow_6step_wrap",
    "SCF012": "figure_bank.scaffolds.SCF012_card_workflow_7step_wrap",
    "SCF013": "figure_bank.scaffolds.SCF013_plot_grid_1panel",
    "SCF014": "figure_bank.scaffolds.SCF014_plot_grid_2panel",
    "SCF015": "figure_bank.scaffolds.SCF015_plot_grid_3panel",
    "SCF016": "figure_bank.scaffolds.SCF016_plot_grid_4panel",
    "SCF017": "figure_bank.scaffolds.SCF017_plot_grid_5panel_ribbon",
    "SCF018": "figure_bank.scaffolds.SCF018_plot_grid_6panel",
    "SCF019": "figure_bank.scaffolds.SCF019_plot_grid_8panel",
    "SCF020": "figure_bank.scaffolds.SCF020_plot_grid_9panel",
    "SCF021": "figure_bank.scaffolds.SCF021_plot_grid_10panel",
    "SCF022": "figure_bank.scaffolds.SCF022_plot_grid_11panel",
    "SCF023": "figure_bank.scaffolds.SCF023_plot_grid_12panel",
    "SCF024": "figure_bank.scaffolds.SCF024_plot_grid_4panel_tall",
    "SCF025": "figure_bank.scaffolds.SCF025_design_panel_overview",
}

_PREVIEW_DIR = Path(__file__).parent / "previews"


def regenerate(targets: Iterable[str] | None = None) -> list[Path]:
    targets = list(targets) if targets else list(_SCAFFOLDS)
    written = []
    _PREVIEW_DIR.mkdir(exist_ok=True)
    for scf in targets:
        if scf not in _MODULES:
            print(f"unknown scaffold: {scf}", file=sys.stderr)
            continue
        mod = importlib.import_module(_MODULES[scf])
        outlined = mod.preview(show_outlines=True)
        clean = mod.preview(show_outlines=False)
        out_path = _PREVIEW_DIR / f"{scf}-preview.svg"
        clean_path = _PREVIEW_DIR / f"{scf}-preview-clean.svg"
        out_path.write_text(outlined, encoding="utf-8")
        clean_path.write_text(clean, encoding="utf-8")
        written.extend([out_path, clean_path])
        print(f"wrote {out_path.relative_to(Path.cwd()) if Path.cwd() in out_path.parents else out_path}")
        print(f"wrote {clean_path.relative_to(Path.cwd()) if Path.cwd() in clean_path.parents else clean_path}")
    return written


def main(argv: list[str]) -> int:
    if not argv or argv[0] != "regenerate":
        print("usage: python -m figure_bank.scaffolds._previews regenerate [SCF### ...]",
              file=sys.stderr)
        return 2
    regenerate(argv[1:])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
