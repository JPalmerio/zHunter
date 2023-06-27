from PyQt6 import uic
from PyQt6 import QtWidgets
from zhunter.initialize import DIRS


class SmoothingWindow(QtWidgets.QDialog):
    def __init__(self, parent):
        super(SmoothingWindow, self).__init__(parent)
        self.parent = parent
        uic.loadUi(DIRS["UI"] / "smoothing.ui", self)
