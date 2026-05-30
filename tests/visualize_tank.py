"""
Quick visualization script for tank geometry.

Run with:
    python tests/visualize_tank.py
"""

import sys
from pathlib import Path

# Allow running directly without installing the package.
sys.path.insert(0, str(Path(__file__).parent.parent))

from openfoam_tank_mesh.Profile import CylinderCapsTankProfile, SphereProfile

FIG_DIR = Path(__file__).parent / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
DPI = 600

# --- Sphere (revolved / axisymmetric) ---
sphere = SphereProfile(
    radius=1.0,
    fill_level=0.4,
    outlet_radius=0.05,
    bulk_cell_size=0.05,
    wall_tan_cell_size=0.02,
    wall_cell_size=0.005,
)
sphere.plot_3d(
    wedge_angle=20,
    save_path=str(FIG_DIR / "sphere_revolved.png"),
    dpi=DPI,
    show=True,
    sim_domain_plane_only=True,
)

# --- Cylinder with ellipsoidal caps (revolved) ---
cylinder = CylinderCapsTankProfile(
    cylinder_radius=0.5,
    cylinder_height=1.0,
    cap_height=0.3,
    fill_level=0.4,
    outlet_radius=0.05,
    bulk_cell_size=0.05,
    wall_tan_cell_size=0.02,
    wall_cell_size=0.005,
)
cylinder.plot_3d(
    wedge_angle=20,
    save_path=str(FIG_DIR / "cylinder_caps_revolved.png"),
    dpi=DPI,
    show=True,
    sim_domain_plane_only=True,
)

# --- Sphere (horizontal / empty_2d) ---
sphere_h = SphereProfile(
    radius=1.0,
    fill_level=0.4,
    outlet_radius=0.05,
    bulk_cell_size=0.05,
    wall_tan_cell_size=0.02,
    wall_cell_size=0.005,
)
sphere_h.plot_3d_horizontal(
    save_path=str(FIG_DIR / "sphere_horizontal.png"),
    dpi=DPI,
    show=True,
    sim_domain_plane_only=True,
)

print(f"Saved figures to {FIG_DIR.resolve()}")
