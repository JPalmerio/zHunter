import pandas as pd
from pathlib import Path
import logging

log = logging.getLogger(__name__)


def create_line_ratios(input_fname, sep=',', output_fname='line_ratio.csv', save=True):
    """
        Takes a line list and calculate all possible ratios between the
        line wavelengths. Then keep only 1 < ratio <= 2 and write them
        to a file.
        Input file must contain 2 columns for the name and the wavelength
    """
    lines = pd.read_csv(input_fname, sep=sep, comment='#')
    ratio = []
    ratio_name = []

    for n1, w1 in zip(lines['name'], lines['wvlg']):
        for n2, w2 in zip(lines['name'], lines['wvlg']):
            ratio_name.append('/'.join([n1.strip(),n2.strip()]))
            ratio.append(w1/w2)

    df_ratios = pd.DataFrame({'ratio':ratio, 'name':ratio_name})
    usable_ratios = (df_ratios['ratio'] > 1) & (df_ratios['ratio'] <= 2)
    df_ratios = df_ratios[usable_ratios]
    df_ratios = df_ratios.sort_values('ratio')

    line_dir = Path(str(input_fname)).resolve().parent
    output_fname = line_dir/output_fname
    if save:
        df_ratios.to_csv(output_fname, index=False)
        log.info("Saved line ratios in {}".format(output_fname))

    return df_ratios
