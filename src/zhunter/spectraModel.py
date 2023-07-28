from PyQt6 import QtGui
from PyQt6 import QtCore
import logging

log = logging.getLogger(__name__)


class SpectraModel(QtCore.QAbstractListModel):
    def __init__(self, *args, spectra=None, **kwargs):
        super(SpectraModel, self).__init__(*args, **kwargs)
        self.spectra = spectra or []

    def data(self, index, role):
        if role == QtCore.Qt.ItemDataRole.DisplayRole:
            # See below for the data structure.
            status, spectrum = self.spectra[index.row()]
            return "{:s}".format(spectrum.name)

        if role == QtCore.Qt.ItemDataRole.ForegroundRole:
            status, spectrum = self.spectra[index.row()]
            return QtGui.QBrush(QtGui.QColor(spectrum.plot_properties['color']))

        if role == QtCore.Qt.ItemDataRole.DecorationRole:
            status, spectrum = self.spectra[index.row()]
            return QtGui.QColor(spectrum.plot_properties['color'])

    # def sort(self):
    #     self.spectra.sort(reverse=True)

    def rowCount(self, index):
        return len(self.spectra)

    def delete(self, index):
        # Remove the item and refresh.
        _, _spec = self.spectra[index.row()]
        log.debug(f"Received request to delete spectrum {_spec.name}")
        _spec.undraw()
        del self.spectra[index.row()]
        self.layoutChanged.emit()
        self.sort()

    def clear(self):
        for _, spectrum in self.spectra:
            spectrum.clear()
        self.spectra = []

    def get_color(self, index):
        _, _spec = self.spectra[index.row()]
        return _spec.plot_properties['color']

    def get_spectrum(self, index):
        _, _spec = self.spectra[index.row()]
        return _spec
