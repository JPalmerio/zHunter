from pathlib import Path
import logging
from .MainGraphicsWidget import MainGraphicsWidget

__version__ = "0.10.3"

log = logging.getLogger(__name__)

__all__ = ["MainGraphicsWidget"]

ROOT_DIR = Path(__file__).parents[0]
DIRS = {
    "ROOT": ROOT_DIR,
    "UI": ROOT_DIR / "ui",
    "DATA": ROOT_DIR / "data",
    "CONFIG": ROOT_DIR / "config",
}
