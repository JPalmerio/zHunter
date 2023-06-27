from PyQt6 import uic
from PyQt6 import QtWidgets
from zhunter.initialize import DIRS


class WavelengthCorrectionWindow(QtWidgets.QDialog):
    def __init__(self, parent):
        super(WavelengthCorrectionWindow, self).__init__(parent)
        self.parent = parent
        uic.loadUi(DIRS["UI"] / "wavelength_correction.ui", self)
