import logging
import sys
from pathlib import Path
from PyQt6 import QtWidgets

from zhunter import __version__
from zhunter.colors import ZHUNTER_LOGO
from zhunter.gui import MainGUI

logging.getLogger("PyQt6").setLevel(logging.INFO)
logging.getLogger("matplotlib").setLevel(logging.INFO)
log = logging.getLogger(__name__)
logging.basicConfig(
    stream=sys.stdout,
    level=logging.DEBUG,
    # format="%(asctime)s.%(msecs)03d | %(levelname)s | [%(name)s] - %(funcName)s : %(message)s",
    format="%(asctime)s.%(msecs)03d | %(levelname)-8s | %(funcName)s - %(filename)s:%(lineno)d : %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def main():
    log.info(
        "\n"
        + 36 * "-"
        + ZHUNTER_LOGO
        + "\n"
        + 14 * " "
        + f"v{__version__}\n"
        + 36 * "-"
    )
    app = QtWidgets.QApplication(sys.argv)
    gui = MainGUI()
    gui.show()
    sys.exit(app.exec())


def testmode1D():
    log.info(
        "\n"
        + 36 * "-"
        + ZHUNTER_LOGO
        + "\n"
        + """
 _ ____    _____ _____ ____ _____   __  __  ___  ____  _____ 
/ |  _ \  |_   _| ____/ ___|_   _| |  \/  |/ _ \|  _ \| ____|
| | | | |   | | |  _| \___ \ | |   | |\/| | | | | | | |  _|  
| | |_| |   | | | |___ ___) || |   | |  | | |_| | |_| | |___ 
|_|____/    |_| |_____|____/ |_|   |_|  |_|\___/|____/|_____|
"""
        + "\n"
        + 14 * " "
        + f"v{__version__}\n"
        + 36 * "-"
    )
    app = QtWidgets.QApplication(sys.argv)
    gui = MainGUI()
    gui.mode = "1D"
    fname = Path(
        str(Path(__file__).resolve().parents[2])
        + "/dev/data/test_input_files/XSHOOTER_bintable_1D.fits"
    )
    gui.config["fnames"]["data"] = fname
    gui.plot()
    gui.add_specsys(z=6.317, sys_type="abs")
    gui.show()
    sys.exit(app.exec())


def testmode2D():
    log.info(
        "\n"
        + 36 * "-"
        + ZHUNTER_LOGO
        + "\n"
        + """
 ____  ____    _____ _____ ____ _____   __  __  ___  ____  _____ 
|___ \|  _ \  |_   _| ____/ ___|_   _| |  \/  |/ _ \|  _ \| ____|
  __) | | | |   | | |  _| \___ \ | |   | |\/| | | | | | | |  _|  
 / __/| |_| |   | | | |___ ___) || |   | |  | | |_| | |_| | |___ 
|_____|____/    |_| |_____|____/ |_|   |_|  |_|\___/|____/|_____|
"""
        + "\n"
        + 14 * " "
        + f"v{__version__}\n"
        + 36 * "-"
    )
    app = QtWidgets.QApplication(sys.argv)
    gui = MainGUI()
    gui.mode = "2D"
    fname = Path(
        str(Path(__file__).resolve().parents[2]) + "/dev/data/test_input_files/2D.fits"
    )
    gui.config["fnames"]["data"] = fname
    gui.plot()
    gui.add_specsys(z=0.15135, sys_type="em")
    gui.show()
    sys.exit(app.exec())
