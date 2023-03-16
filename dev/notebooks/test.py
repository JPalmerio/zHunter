import pyqtgraph as pg
from PyQt6 import QtWidgets
from PyQt6 import QtCore
from PyQt6 import QtGui
import sys
import numpy as np
from zhunter.misc import set_up_linked_vb
from itertools import cycle
import logging
from collections import defaultdict
from zhunter.colors import COLORS

color_style = "kraken9"
colors = COLORS[color_style]


log = logging.getLogger(__name__)
logging.basicConfig(
    stream=sys.stdout,
    level=logging.DEBUG,
    format="%(levelname)s [%(name)s] %(message)s",
)


ALLOWED_MODES = ["1D", "2D"]

qt_keys = (
    (getattr(QtCore.Qt.Key, attr), attr[4:])
    for attr in dir(QtCore.Qt.Key)
    if attr.startswith("Key_")
)
keys_mapping = defaultdict(lambda: "unknown", qt_keys)


class KeyPressViewBox(pg.ViewBox):
    """
    A custom ViewBox that can catch key presses
    and return the mouse position at the time
    of the keyPressEvent.
    """

    sigKeyPress = QtCore.pyqtSignal(str, float, float)

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.mousePoint = QtCore.QPointF()
        # StrongFocus necessary otherwise doesn't catch events
        self.setFocusPolicy(QtCore.Qt.FocusPolicy.StrongFocus)

    def itemChange(self, change, value):
        if change == self.GraphicsItemChange.ItemSceneChange and value:
            # automatically connect the signal when added to a scene
            value.sigMouseMoved.connect(self.mouse_moved)
            self.setFocus()
        return super().itemChange(change, value)

    def mouse_moved(self, pos):
        self.mousePoint = self.mapSceneToView(pos)

    def keyPressEvent(self, ev):
        key = keys_mapping[ev.key()]
        log.debug(
            "KeyPressViewBox caught keyPressEvent! "
            f"Key: {key}, Mouse position: [{self.mousePoint.x()},{self.mousePoint.y()}]"
        )
        self.sigKeyPress.emit(key, self.mousePoint.x(), self.mousePoint.y())


def __set_up_fit_plot(self, name):
        # Reuse most of the 2D setup but add a plot for
        # residuals
        self.__set_up_2D_plot(name=name)
        self.ax_res = self.addPlot(row=3, col=0, name="residuals")
        self.vb_resh = self.addViewBox(row=3, col=1, name="residuals histogram")

        # No need to adjust plot proportions because
        # it is called in __set_up_2D_plot, but
        # in fit mode, there is an additional residual plot
        # below 1D plot, so adjust it accordingly
        self.__adjust_res_axis()

        # Link plot views together
        self.ax_res.vb.setXLink(self.ax1D.vb)
        self.vb_resh.setYLink(self.ax_res.vb)

        # Remove mouse interactions on side histogram
        self.vb_resh.setMouseEnabled(x=False, y=True)

    # Placeholders
    # self.flux_1D_res = pg.PlotCurveItem(
    #     np.zeros(2),
    #     np.zeros(1),
    #     stepMode="center",
    #     pen=pg.mkPen(color=self.colors["spec"]),
    # )
    # self.flux_1D_resh = pg.PlotCurveItem(
    #     np.zeros(1),
    #     np.zeros(1),
    #     pen=pg.mkPen(color=self.colors["spec"]),
    # )
    # self.ax_res.addItem(self.flux_1D_res)
    # self.vb_resh.addItem(self.flux_1D_resh)
    #
    # Adjust plots
    # # Change the ratios of sizes of PlotItems
    # self.ci.layout.setRowStretchFactor(1, 2)
    # self.ci.layout.setRowStretchFactor(2, 5)
    # self.ci.layout.setRowStretchFactor(3, 2)
    # # Strech column 0 (where 1D and 2D plots are) to make it bigger in x than the side histograms
    # self.ci.layout.setColumnStretchFactor(0, 5)

    # def __adjust_res_axis(self):
    #     Hide the axis to make the plot prettier and more compact
    #     self.ax1D.hideAxis("bottom")
    #     # Change all the widths to be equal so things are aligned
    #     self.ax_res.getAxis("left").setWidth(60)
    #     self.ax_res.getAxis("bottom").setHeight(30)

    #     # Remove padding so that panning with keyboard preserves x and y range
    #     self.ax_res.vb.setDefaultPadding(padding=0.00)
    #     self.ax_res.showGrid(x=True, y=True)


def main():
    app = QtWidgets.QApplication(sys.argv)
    main_window = QtWidgets.QMainWindow()
    main_window.setBaseSize(1600, 900)
    win = MainPlotWidget(main_window)
    main_window.setCentralWidget(win)
    main_window.show()
    win.set_up_plot(name="test", mode="1D", colors=colors)
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
