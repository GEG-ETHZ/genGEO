
# compute relative error between two values and compare to threshold
def check_rel_error(eps, val, ref):
    return abs(((val - ref) / ref)) < eps
