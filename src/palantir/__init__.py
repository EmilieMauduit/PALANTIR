# src/__init__.py

from .prediction_tools import *
from .output_analysis import *
from .scripts import *
from .logger import *

__all__ = [
    "prediction_tools",
    "output_analysis",
    "scripts",
    "setup_logging"
]

__author__ = "Emilie Mauduit"
__copyright__ = "Copyright 2023, palantir"
__credits__ = ["Emilie Mauduit", "Jean-Mathias Griessmeier"]
__license__ = "MIT"
__version__ = "0.1.0"
__maintainer__ = "Emilie Mauduit"
__email__ = "emilie.mauduit@obspm.fr"
