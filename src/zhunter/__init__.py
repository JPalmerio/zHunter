from pathlib import Path
from zhunter.MainGraphicsWidget import MainGraphicsWidget

__version__ = "1.0.0"

__all__ = ["MainGraphicsWidget"]

__ROOT_DIR__ = Path(__file__).parent
ZHUNTER_DIR = Path(__file__).parent
DIRS = {"ROOT": __ROOT_DIR__, "UI": __ROOT_DIR__ / "ui", "DATA": __ROOT_DIR__ / "data"}
