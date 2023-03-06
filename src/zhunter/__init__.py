from pathlib import Path
import logging

__version__ = "0.10.2"

log = logging.getLogger(__name__)

ROOT_DIR = Path(__file__).parents[0]
DIRS = {"ROOT": ROOT_DIR, "UI": ROOT_DIR / "ui", "DATA": ROOT_DIR / "data"}
log.debug("DIRS: {}".format(DIRS))
