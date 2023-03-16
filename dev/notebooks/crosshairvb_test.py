from PyQt6.QtWidgets import QApplication
from PyQt6 import QtWidgets
from PyQt6 import QtGui
from PyQt6 import QtCore
from PyQt6.QtCore import pyqtSignal

from collections import defaultdict
import pyqtgraph as pg
import logging
import sys

log = logging.getLogger(__name__)
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s [%(name)s] %(message)s')




# Get key mappings from Qt namespace
qt_keys = (
    (getattr(QtCore.Qt.Key, attr), attr[4:])
    for attr in dir(QtCore.Qt.Key)
    if attr.startswith("Key_")
)
keys_mapping = defaultdict(lambda: "unknown", qt_keys)


def keyPressHandler(event):
    key = keys_mapping[event.key()]
    log.info(f"You pressed the key: {key}")


class CrosshairViewBox(pg.ViewBox):
    # sigKeyPress = QtCore.pyqtSignal(object)

    def __init__(self, *args, **kwds):
        pg.ViewBox.__init__(self, *args, **kwds)
    
        self.crosshair_x = pg.InfiniteLine(
            angle=90, movable=False
        )
        self.crosshair_y = pg.InfiniteLine(
            angle=0, movable=False
        )
        self.addItem(self.crosshair_x, ignoreBounds=True)
        self.addItem(self.crosshair_y, ignoreBounds=True)

    def activate_crosshairs(self):
        self.scene().sigMouseMoved.connect(self.move_crosshair)

    def move_crosshair(self, ev):
        """
            This is triggered on a sigMouseMoved which sends a position as an event
        """
        pos = ev
        if self.sceneBoundingRect().contains(pos):
            mousePoint = self.mapSceneToView(ev)
            self.crosshair_x.setPos(mousePoint.x())
            self.crosshair_y.setPos(mousePoint.y())

    def keyPressEvent(self, ev):
        """
            In this case, ev is a QKeyEvent, for which I cannot get
            the position of the mouse. So I have to use the crosshairs.
        """
        self.scene().keyPressEvent(ev)
        posx = self.crosshair_x.getPos()[0]
        posy = self.crosshair_y.getPos()[1]
        print(f"Mouse position from crosshairs: [{posx}, {posy}]")
        # self.sigKeyPress.emit(ev)


def main():
    app = QtWidgets.QApplication(sys.argv)
    win = pg.GraphicsLayoutWidget(show=True)
    win.ci.setBorder((50, 50, 100)) # this is to see where the Items' bounds are
    chvb = CrosshairViewBox()
    # chvb.sigKeyPress.connect(keyPressHandler)
    win.addItem(chvb)
    chvb.activate_crosshairs()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()