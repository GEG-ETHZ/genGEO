
class Results(object):
    """Results is an empty object to be filled with results."""
    def __init__(self, fluid):
        self.end_T_C = None
        self.fluid =  fluid

    def finalState(self):
        return FluidState.getStateFromTQ(self.end_T_C, 0, self.fluid)
