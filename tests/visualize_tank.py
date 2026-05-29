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

# --- Sphere (revolved / axisymmetric) ---
sphere = SphereProfile(
    radius=1.0,
    fill_level=0.5,
    outlet_radius=0.05,
    bulk_cell_size=0.05,
    wall_tan_cell_size=0.02,
    wall_cell_size=0.005,
)
sphere.plot_3d(wedge_angle=30)

# --- Cylinder with ellipsoidal caps (revolved) ---
cylinder = CylinderCapsTankProfile(
    cylinder_radius=0.5,
    cylinder_height=1.0,
    cap_height=0.3,
    fill_level=0.5,
    outlet_radius=0.05,
    bulk_cell_size=0.05,
    wall_tan_cell_size=0.02,
    wall_cell_size=0.005,
)
cylinder.plot_3d(wedge_angle=30)

# --- Sphere (horizontal / empty_2d) ---
sphere_h = SphereProfile(
    radius=1.0,
    fill_level=0.2,
    outlet_radius=0.05,
    bulk_cell_size=0.05,
    wall_tan_cell_size=0.02,
    wall_cell_size=0.005,
)
sphere_h.plot_3d_horizontal()
