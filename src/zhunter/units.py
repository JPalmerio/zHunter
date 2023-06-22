from PyQt6 import uic
from PyQt6 import QtWidgets
from zhunter import DIRS
import astropy.units as u
from astropy.units.core import UnitConversionError
import logging

log = logging.getLogger(__name__)

# To allow converting e.g. from Hz to nm
u.set_enabled_equivalencies(u.spectral())


class UnitsWindow(QtWidgets.QWidget):
    def __init__(self, data):
        super().__init__()
        uic.loadUi(DIRS["UI"] / "units.ui", self)

        self.data = data
        self.units = {
            'wvlg': data['wvlg'].unit,
            'flux_1D': data['flux_1D'].unit,
            'flux_2D': data['flux_2D'].unit,
            'spat': data['spat'].unit,
        }

        self.cbb_quantity.currentIndexChanged.connect(self.update_current_units_disp)
        self.cbb_quantity.addItems(list(self.units.keys()))

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
        quant = str(self.cbb_quantity.currentText())

        new_units = self.get_valid_units(textbox=self.txb_set_units)

        # If user provided units are valid
        if new_units:
            self.data[quant] = self.data[quant].value * new_units
            self.units[quant] = new_units

            self.update_current_units_disp()

    def convert_units(self):

        quant = str(self.cbb_quantity.currentText())

        current_units = u.Unit(self.units[quant])
        new_units = self.get_valid_units(textbox=self.txb_convert_units)

        if new_units:
            log.debug(f"Converting {current_units} to {new_units}")
            try:
                self.data[quant] = self.data[quant].to(
                    new_units,
                    # Add equivalency in case converting Fnu/Flam
                    equivalencies=u.spectral_density(self.data['wvlg'])
                )
                self.units[quant] = new_units
                self.update_current_units_disp()

            except UnitConversionError as e:

                log.warning(f"Invalid unit conversion: {e}")
                QtWidgets.QMessageBox.warning(
                    self, "Invalid unit conversion", f"{e}"
                )

            # self.data.set_1D_displayed()

    def update_current_units_disp(self):
        ax = str(self.cbb_quantity.currentText())

        # self.txb_current_unit.setText(str(self.units[ax]))
        self.current_units_label.setText(str(self.units[ax]))
