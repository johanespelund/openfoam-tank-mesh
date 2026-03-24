from openfoam_tank_mesh.CylinderTank import CylinderTank

INCH = 0.0254


class KSiteTank(CylinderTank):
    """
    A class to represent the K-Site tank.
    See Stochl & Knoll (1991) for details:
        https://ntrs.nasa.gov/citations/19910015845
    """

    def __init__(self, fill_level: float, outlet_radius: float) -> None:
        super().__init__(
            name="K-Site",
            fill_level=fill_level,
            outlet_radius=outlet_radius,
            cylinder_radius=0.5 * 87.6 * INCH,
            cylinder_height=1.5 * INCH,
            cap_height=0.5 * 73 * INCH,
        )
        return None
