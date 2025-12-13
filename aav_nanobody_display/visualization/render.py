"""
Visualization module using py3Dmol.

Generates interactive HTML visualizations of AAV capsids with nanobody insertions.
"""

import logging
from pathlib import Path
from typing import Optional, Dict, List, Any, Literal
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)

try:
    import py3Dmol
    HAS_PY3DMOL = True
except ImportError:
    HAS_PY3DMOL = False
    logger.warning("py3Dmol not installed. Install with: pip install py3dmol")


# Color schemes
COLOR_SCHEMES = {
    "capsid_default": {
        "vp3_core": "lightgray",
        "vr_loops": "steelblue",
        "nanobody": "firebrick",
        "linker": "gold",
    },
    "publication": {
        "vp3_core": "#E8E8E8",
        "vr_loops": "#4682B4",
        "nanobody": "#DC143C",
        "linker": "#FFD700",
    },
    "colorblind_safe": {
        "vp3_core": "#DDDDDD",
        "vr_loops": "#0072B2",
        "nanobody": "#D55E00",
        "linker": "#F0E442",
    },
}


@dataclass
class VisualizationConfig:
    """Configuration for capsid visualization."""
    width: int = 1000
    height: int = 800
    background_color: str = "white"
    color_scheme: str = "capsid_default"
    show_surface: bool = False
    surface_opacity: float = 0.5
    cartoon_style: bool = True
    show_labels: bool = True
    label_font_size: int = 14


class CapsidVisualizer:
    """
    Create interactive 3D visualizations of AAV capsids.
    
    Generates standalone HTML files that can be opened in any browser.
    """
    
    def __init__(self, config: Optional[VisualizationConfig] = None):
        """
        Initialize the visualizer.
        
        Args:
            config: Visualization configuration
        """
        if not HAS_PY3DMOL:
            raise ImportError(
                "py3Dmol is required for visualization. "
                "Install with: pip install py3dmol"
            )
        
        self.config = config or VisualizationConfig()
        self.colors = COLOR_SCHEMES.get(
            self.config.color_scheme, 
            COLOR_SCHEMES["capsid_default"]
        )
    
    def create_view(self) -> "py3Dmol.view":
        """Create a new py3Dmol view with configured dimensions."""
        view = py3Dmol.view(
            width=self.config.width,
            height=self.config.height
        )
        view.setBackgroundColor(self.config.background_color)
        return view
    
    def visualize_monomer(
        self,
        pdb_path: Path,
        vr_loops: Optional[Dict[str, tuple]] = None,
        highlight_residues: Optional[List[int]] = None,
        output_path: Optional[Path] = None
    ) -> "py3Dmol.view":
        """
        Visualize a single VP monomer with optional VR loop highlighting.
        
        Args:
            pdb_path: Path to monomer PDB
            vr_loops: Dict of {loop_name: (start, end)} to highlight
            highlight_residues: Additional residues to highlight
            output_path: If provided, save HTML to this path
            
        Returns:
            py3Dmol view object
        """
        view = self.create_view()
        
        # Load structure
        with open(pdb_path) as f:
            pdb_data = f.read()
        view.addModel(pdb_data, "pdb")
        
        # Base style - cartoon
        if self.config.cartoon_style:
            view.setStyle({}, {"cartoon": {"color": self.colors["vp3_core"]}})
        else:
            view.setStyle({}, {"stick": {"color": self.colors["vp3_core"]}})
        
        # Highlight VR loops
        if vr_loops:
            for loop_name, (start, end) in vr_loops.items():
                view.setStyle(
                    {"resi": f"{start}-{end}"},
                    {"cartoon": {"color": self.colors["vr_loops"]}}
                )
                
                if self.config.show_labels:
                    # Add label at midpoint
                    midpoint = (start + end) // 2
                    view.addLabel(
                        loop_name,
                        {
                            "position": {"resi": midpoint},
                            "backgroundColor": self.colors["vr_loops"],
                            "fontColor": "white",
                            "fontSize": self.config.label_font_size,
                        }
                    )
        
        # Highlight specific residues
        if highlight_residues:
            for res in highlight_residues:
                view.setStyle(
                    {"resi": res},
                    {"cartoon": {"color": "yellow"}, "stick": {}}
                )
        
        # Add surface if requested
        if self.config.show_surface:
            view.addSurface(
                py3Dmol.VDW,
                {"opacity": self.config.surface_opacity, "color": "white"}
            )
        
        view.zoomTo()
        
        # Save HTML if output path provided
        if output_path:
            self._save_html(view, output_path)
        
        return view
    
    def visualize_capsid_with_nanobody(
        self,
        capsid_pdb: Path,
        nanobody_chain: str = "B",
        capsid_chain: str = "A",
        insertion_site: Optional[tuple] = None,
        output_path: Optional[Path] = None
    ) -> "py3Dmol.view":
        """
        Visualize capsid monomer with attached nanobody.
        
        Args:
            capsid_pdb: Path to combined PDB (capsid + nanobody)
            nanobody_chain: Chain ID of nanobody
            capsid_chain: Chain ID of capsid
            insertion_site: (start, end) residues of insertion site
            output_path: If provided, save HTML
            
        Returns:
            py3Dmol view object
        """
        view = self.create_view()
        
        with open(capsid_pdb) as f:
            pdb_data = f.read()
        view.addModel(pdb_data, "pdb")
        
        # Style capsid
        view.setStyle(
            {"chain": capsid_chain},
            {"cartoon": {"color": self.colors["vp3_core"]}}
        )
        
        # Style nanobody
        view.setStyle(
            {"chain": nanobody_chain},
            {"cartoon": {"color": self.colors["nanobody"]}}
        )
        
        # Highlight insertion site
        if insertion_site:
            start, end = insertion_site
            view.setStyle(
                {"chain": capsid_chain, "resi": f"{start}-{end}"},
                {"cartoon": {"color": self.colors["linker"]}}
            )
        
        # Labels
        if self.config.show_labels:
            view.addLabel(
                "Nanobody",
                {
                    "position": {"chain": nanobody_chain, "resi": 60},
                    "backgroundColor": self.colors["nanobody"],
                    "fontColor": "white",
                    "fontSize": self.config.label_font_size,
                }
            )
            view.addLabel(
                "VP",
                {
                    "position": {"chain": capsid_chain, "resi": 400},
                    "backgroundColor": self.colors["vp3_core"],
                    "fontColor": "black",
                    "fontSize": self.config.label_font_size,
                }
            )
        
        view.zoomTo()
        
        if output_path:
            self._save_html(view, output_path)
        
        return view
    
    def visualize_full_capsid(
        self,
        capsid_pdb: Path,
        highlight_chains: Optional[List[str]] = None,
        nanobody_chains: Optional[List[str]] = None,
        output_path: Optional[Path] = None,
        surface_only: bool = False
    ) -> "py3Dmol.view":
        """
        Visualize full 60-mer capsid.
        
        Args:
            capsid_pdb: Path to full capsid PDB
            highlight_chains: Chains to highlight (e.g., show ASU)
            nanobody_chains: Chains that are nanobodies
            output_path: If provided, save HTML
            surface_only: If True, show only surface representation
            
        Returns:
            py3Dmol view object
        """
        view = self.create_view()
        
        with open(capsid_pdb) as f:
            pdb_data = f.read()
        view.addModel(pdb_data, "pdb")
        
        if surface_only:
            view.setStyle({}, {"sphere": {"radius": 0.5, "color": self.colors["vp3_core"]}})
        else:
            # Base cartoon style
            view.setStyle({}, {"cartoon": {"color": self.colors["vp3_core"]}})
        
        # Highlight specific chains
        if highlight_chains:
            for chain in highlight_chains:
                view.setStyle(
                    {"chain": chain},
                    {"cartoon": {"color": self.colors["vr_loops"]}}
                )
        
        # Style nanobody chains differently
        if nanobody_chains:
            for chain in nanobody_chains:
                view.setStyle(
                    {"chain": chain},
                    {"cartoon": {"color": self.colors["nanobody"]}}
                )
        
        view.zoomTo()
        
        if output_path:
            self._save_html(view, output_path)
        
        return view
    
    def compare_with_without_nanobody(
        self,
        capsid_only_pdb: Path,
        capsid_with_nb_pdb: Path,
        output_path: Optional[Path] = None
    ) -> str:
        """
        Create a side-by-side comparison HTML.
        
        Args:
            capsid_only_pdb: Path to capsid without nanobody
            capsid_with_nb_pdb: Path to capsid with nanobody
            output_path: Where to save HTML
            
        Returns:
            HTML string
        """
        # Read PDB data
        with open(capsid_only_pdb) as f:
            pdb1 = f.read().replace("'", "\\'").replace("\n", "\\n")
        with open(capsid_with_nb_pdb) as f:
            pdb2 = f.read().replace("'", "\\'").replace("\n", "\\n")
        
        html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>AAV Capsid Comparison</title>
    <script src="https://3dmol.org/build/3Dmol-min.js"></script>
    <style>
        body {{ margin: 0; font-family: Arial, sans-serif; }}
        .container {{ display: flex; width: 100%; height: 100vh; }}
        .viewer {{ flex: 1; position: relative; border: 1px solid #ccc; }}
        .label {{ position: absolute; top: 10px; left: 10px; 
                  background: rgba(0,0,0,0.7); color: white; 
                  padding: 5px 10px; border-radius: 4px; z-index: 100; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="viewer" id="viewer1">
            <div class="label">Without Nanobody</div>
        </div>
        <div class="viewer" id="viewer2">
            <div class="label">With Nanobody</div>
        </div>
    </div>
    <script>
        // Viewer 1 - Capsid only
        let viewer1 = $3Dmol.createViewer('viewer1', {{backgroundColor: 'white'}});
        viewer1.addModel('{pdb1}', 'pdb');
        viewer1.setStyle({{}}, {{cartoon: {{color: '{self.colors["vp3_core"]}'}}}});
        viewer1.zoomTo();
        viewer1.render();
        
        // Viewer 2 - Capsid with nanobody
        let viewer2 = $3Dmol.createViewer('viewer2', {{backgroundColor: 'white'}});
        viewer2.addModel('{pdb2}', 'pdb');
        viewer2.setStyle({{chain: 'A'}}, {{cartoon: {{color: '{self.colors["vp3_core"]}'}}}});
        viewer2.setStyle({{chain: 'B'}}, {{cartoon: {{color: '{self.colors["nanobody"]}'}}}});
        viewer2.zoomTo();
        viewer2.render();
        
        // Sync rotation between viewers
        viewer1.linkViewer(viewer2);
        viewer2.linkViewer(viewer1);
    </script>
</body>
</html>
"""
        
        if output_path:
            output_path = Path(output_path)
            with open(output_path, "w") as f:
                f.write(html)
            logger.info(f"Saved comparison to: {output_path}")
        
        return html
    
    def _save_html(self, view: "py3Dmol.view", output_path: Path):
        """Save a py3Dmol view as standalone HTML."""
        output_path = Path(output_path)
        
        # Get the HTML from py3Dmol
        html = view._make_html()
        
        # Wrap in full HTML document
        full_html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>AAV Capsid Visualization</title>
    <script src="https://3dmol.org/build/3Dmol-min.js"></script>
</head>
<body style="margin: 0;">
    {html}
</body>
</html>
"""
        
        with open(output_path, "w") as f:
            f.write(full_html)
        
        logger.info(f"Saved visualization to: {output_path}")


def visualize_capsid(
    pdb_path: Path,
    output_html: Optional[Path] = None,
    vr_loops: Optional[Dict[str, tuple]] = None,
    **config_kwargs
) -> "py3Dmol.view":
    """
    Convenience function to quickly visualize a capsid structure.
    
    Args:
        pdb_path: Path to PDB file
        output_html: Where to save HTML (optional)
        vr_loops: VR loops to highlight {name: (start, end)}
        **config_kwargs: Additional configuration options
        
    Returns:
        py3Dmol view object
    """
    config = VisualizationConfig(**config_kwargs)
    viz = CapsidVisualizer(config)
    return viz.visualize_monomer(pdb_path, vr_loops, output_path=output_html)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    print("Visualization module loaded.")
    print(f"py3Dmol available: {HAS_PY3DMOL}")
    print("\nUsage:")
    print("  viz = CapsidVisualizer()")
    print("  viz.visualize_monomer('capsid.pdb', output_path='view.html')")
