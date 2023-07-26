from PyQt6 import uic
from PyQt6 import QtWidgets
from zhunter.initialize import DIRS
import astropy.units as u

from astropy.units.core import UnitConversionError
import logging

from zhunter.data_handler import DataHandler

log = logging.getLogger(__name__)


class UnitsWindow(QtWidgets.QDialog):

    def __init__(self, parent):
        super(UnitsWindow, self).__init__(parent)
        self.parent = parent
        uic.loadUi(DIRS["UI"] / "units.ui", self)

        self.data = None
        self.connect_signals_and_slots()

    def load_data(self, data):
        """Load a DataHandler instance

        Parameters
        ----------
        data : DataHandler
            DataHandler instance containing the data
        """
        if not isinstance(data, DataHandler):
            raise TypeError("Data provided must be a DataHandler instance")

        self.data = data
        self.data.sigUnitsUpdated.connect(self.refresh)
        self.refresh()

    def refresh(self):
        """Refresh combobox"""
        if self.data is None:
            log.warning("Refreshed units window without loading data")
            return

        self.cbb_axis.clear()
        self.cbb_axis.addItems(list(self.data.units.keys()))

    def clear(self):
        self.cbb_axis.clear()
        self.data = None

    def connect_signals_and_slots(self):
        self.cbb_axis.currentIndexChanged.connect(self.update_current_units_disp)
        self.set_units_btn.clicked.connect(self.set_units)
        self.convert_units_btn.clicked.connect(self.convert_units)

    def get_valid_units(self, textbox):
        """Read provided textbox, strip string and
        check that it is a valid unit.
        If not, return None

        Parameters
        ----------
        textbox : PyQt6.QtWidgets.QLineEdit
            Where to read the string to check if valid.

        Returns
        -------
        Unit
            Validated units or None
        """
        try:
            unit_str = str(textbox.text()).strip()
            if unit_str:
                valid_units = u.Unit(unit_str)
                return valid_units
            else:
                return None

        except ValueError as e:
            log.warning(f"Invalid units: {e}")
            QtWidgets.QMessageBox.warning(
                self, "Invalid units", f"{e}"
            )
            return None

    def set_units(self):
        """Set the units of the currently selected
        axis in the combobox to what is specified in
        the appropriate QLineEdit
        """
        axis = str(self.cbb_axis.currentText())

        # Avoid case where combobox is empty
        if not axis:
            return

        new_units = self.get_valid_units(textbox=self.txb_set_units)

        # If user provided units are valid
        if new_units:
            log.info(f"Redefining {axis} axis units as {new_units}")
            self.data.set_units({axis: new_units})
            self.update_current_units_disp()

    def convert_units(self):
        """Convert units of the currently selected
        axis in the combobox to what is specified in
        the appropriate QLineEdit
        """
        axis = str(self.cbb_axis.currentText())

        # Avoid case where combobox is empty
        if not axis:
            return

        current_units = u.Unit(self.data.units[axis])
        new_units = self.get_valid_units(textbox=self.txb_convert_units)

        if new_units:
            log.info(f"Converting '{axis}' axis from '{current_units}' to '{new_units}'")
            try:
                self.data.convert_units(axis, new_units)
                self.update_current_units_disp()
            except UnitConversionError as e:
                log.warning(f"Invalid unit conversion: {e}")
                QtWidgets.QMessageBox.warning(
                    self, "Invalid unit conversion", f"{e}"
                )

    def update_current_units_disp(self):
        axis = str(self.cbb_axis.currentText())
        # Avoid case where combobox is empty
        if not axis:
            self.current_units_label.clear()
        else:
            self.current_units_label.setText(str(self.data.units[axis]))
