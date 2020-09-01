


class GlobalPhysicalProperties(object):
    """GlobalPhysicalProperties provides physical properties of the system."""

    # assign default values
    def __init__(self):
        self.g = 9.81                   # m/s**2
        self.density_rock = 2650.       # kg/m**3
        self.C_rock = 1000.             # J/kg-K
        self.k_rock = 2.1               # W/m-K
