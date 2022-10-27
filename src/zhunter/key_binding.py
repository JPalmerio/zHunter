from PyQt5 import uic
from PyQt5 import QtWidgets
from zhunter import DIRS


class KeyBindingHelpDialog(QtWidgets.QDialog):
    def __init__(self,parent):
        super(KeyBindingHelpDialog, self).__init__(parent)
        self.parent = parent
        uic.loadUi(DIRS['UI']/'key_bindings.ui', self)
