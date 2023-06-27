from PyQt6 import uic
from PyQt6 import QtWidgets
from zhunter.initialize import DIRS
import astropy.units as u
from astropy.units.core import UnitConversionError
import logging

log = logging.getLogger(__name__)

# To allow converting e.g. from Hz to nm
u.set_enabled_equivalencies(u.spectral())


class UnitsWindow(QtWidgets.QDialog):
    def __init__(self, parent):
        super(UnitsWindow, self).__init__(parent)
        self.parent = parent
        uic.loadUi(DIRS["UI"] / "units.ui", self)

    def activate(self, data):
        self.data = data

        self.cbb_axis.currentIndexChanged.connect(self.update_current_units_disp)
        self.cbb_axis.addItems(list(self.data["units"].keys()))

        self.connect_signals_and_slots()

    def connect_signals_and_slots(self):
        self.set_units_btn.clicked.connect(self.set_units)
        self.convert_units_btn.clicked.connect(self.convert_units)

    def get_valid_units(self, textbox):
        try:
            valid_units = u.Unit(str(textbox.text()))
            return valid_units

        except ValueError as e:
            log.warning(f"Invalid units: {e}")
            QtWidgets.QMessageBox.warning(
                self, "Invalid units", f"{e}"
            )
            return None

    def set_units(self):
        axis = str(self.cbb_axis.currentText())

        new_units = self.get_valid_units(textbox=self.txb_set_units)

        # If user provided units are valid
        if new_units:
            log.info(f"Redefining {axis} axis units as {new_units}")
            self.data.set_unit({axis: new_units})

            self.update_current_units_disp()

    def convert_units(self):

        axis = str(self.cbb_axis.currentText())

        current_units = u.Unit(self.data["units"][axis])
        new_units = self.get_valid_units(textbox=self.txb_convert_units)

        if new_units:
            log.debug(f"Converting {axis} units from {current_units} to {new_units}")
            raise NotImplementedError
            # try:
            #     self.data[quant] = self.data[quant].to(
            #         new_units,
            #         # Add equivalency in case converting Fnu/Flam
            #         equivalencies=u.spectral_density(self.data['wvlg'])
            #     )
            #     self.data["units"][quant] = new_units
            #     self.update_current_units_disp()

            # except UnitConversionError as e:
            #     log.warning(f"Invalid unit conversion: {e}")
            #     QtWidgets.QMessageBox.warning(
            #         self, "Invalid unit conversion", f"{e}"
            #     )

            # self.data.set_1D_displayed()

    def update_current_units_disp(self):
        axis = str(self.cbb_axis.currentText())

        self.current_units_label.setText(str(self.data["units"][axis]))
