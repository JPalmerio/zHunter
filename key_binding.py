from PyQt5 import uic
from PyQt5 import QtWidgets
from pathlib import Path

ROOT_DIR = Path(__file__).parent.resolve()


class KeyBindingHelpDialog(QtWidgets.QDialog):
    def __init__(self,parent):
        super(KeyBindingHelpDialog, self).__init__(parent)
        self.parent = parent
        uic.loadUi(ROOT_DIR/'key_bindings.ui', self)
