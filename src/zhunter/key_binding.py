from PyQt6 import uic
from PyQt6 import QtWidgets
from zhunter.initialize import DIRS


class KeyBindingHelpDialog(QtWidgets.QDialog):
    def __init__(self, parent):
        super(KeyBindingHelpDialog, self).__init__(parent)
        self.parent = parent
        uic.loadUi(DIRS["UI"] / "key_bindings.ui", self)
