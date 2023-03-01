"""
    Adapted for PyQt6 from https://github.com/yjg30737/pyqt-checkbox-list-widget
"""

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import QListWidgetItem
from PyQt6 import QtWidgets


class ToolTipListWidget(QtWidgets.QListWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Note that *args, **kwargs is necessary to be able to
        # be promoted from Qt Designer
        self.__initUi()

    def __initUi(self):
        self.setMouseTracking(True)
        self.itemEntered.connect(self.__showToolTip)

    def addItem(self, aitem) -> None:
        item = ""
        text = ""
        if isinstance(aitem, str):
            item = QListWidgetItem()
            text = aitem
            item.setText(text)
        elif isinstance(aitem, QListWidgetItem):
            item = aitem
            text = item.text()
        self.setItemWidget(item, QtWidgets.QWidget())
        super().addItem(item)

    def __showToolTip(self, item: QListWidgetItem):
        text = item.text()
        text_width = self.fontMetrics().boundingRect(text).width()
        width = self.width()
        if text_width > width:
            item.setToolTip(text)
        else:
            item.setToolTip("")


class CheckBoxListWidget(ToolTipListWidget):
    checkedSignal = pyqtSignal(int, Qt.CheckState)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.itemChanged.connect(self.__sendCheckedSignal)

    def __sendCheckedSignal(self, item):
        r_idx = self.row(item)
        state = item.checkState()
        self.checkedSignal.emit(r_idx, state)

    def addItems(self, items) -> None:
        for item in items:
            self.addItem(item)

    def addItem(self, item) -> None:
        if isinstance(item, str):
            item = QListWidgetItem(item)
        item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
        item.setCheckState(Qt.CheckState.Unchecked)
        super().addItem(item)

    def toggleState(self, state):
        for i in range(self.count()):
            item = self.item(i)
            if item.checkState() != state:
                item.setCheckState(state)

    def getCheckedRows(self):
        return self.__getFlagRows(Qt.CheckState.Checked)

    def getUncheckedRows(self):
        return self.__getFlagRows(Qt.CheckState.Unchecked)

    def __getFlagRows(self, flag: Qt.CheckState):
        flag_lst = []
        for i in range(self.count()):
            item = self.item(i)
            if item.checkState() == flag:
                flag_lst.append(i)

        return flag_lst

    def removeCheckedRows(self):
        self.__removeFlagRows(Qt.CheckState.Checked)

    def removeUncheckedRows(self):
        self.__removeFlagRows(Qt.CheckState.Unchecked)

    def __removeFlagRows(self, flag):
        flag_lst = self.__getFlagRows(flag)
        flag_lst = reversed(flag_lst)
        for i in flag_lst:
            self.takeItem(i)
