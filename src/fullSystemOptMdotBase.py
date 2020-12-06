

class FullSystemOptMdotBase(object):
    """
    FullSystemSolverBase provides a solver to determine the optimum flow
    rate for a minimum or maximum of a given output variable.
    """

    def __init__(self, system):
        self.full_system = system

    def getTargetVar(self, system):
        raise Exception('GenGeo::no target variable provided to find Minimum '
                        ' or Maximum of!')

    def getDirection(self):
        raise Exception('GenGeo::no direction provided to find Minimum '
                        ' or Maximum getDirection must be '
                        '-1 for minimum or'
                        ' 1 for maximum!')

    def solve(self):
        raise Exception('GenGeo::no algorithm provided to solve opt m_dot '
                        'for minimum or maximum of output variable!')
