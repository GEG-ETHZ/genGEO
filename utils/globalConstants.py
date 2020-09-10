from pathlib import Path

def getProjectRoot() -> Path:
    return Path(__file__).parent.parent


class globalConstants(object):
    """globalConstants carries global constants to convert units."""

    secPerYear = 3600. * 24. * 365.25
    kelvin2celsius = 273.15
