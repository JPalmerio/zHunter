import shutil
import logging
import yaml

from pathlib import Path
import zhunter.io as io
from .misc import create_line_ratios
from astropy.io.ascii import read as ascii_read
from zhunter.colors import COLORS

log = logging.getLogger(__name__)

ROOT_DIR = Path(__file__).parents[0]

DIRS = {
    "ROOT": ROOT_DIR,
    "UI": ROOT_DIR / "ui",
    "DATA": ROOT_DIR / "data",
    "CONFIG": ROOT_DIR / "config",
}

LINE_TYPES = ['intervening', 'emission', 'GRB']

line_dir = DIRS["DATA"] / "lines/"


def get_config_fname():
    """Find the configuration file for zhunter.
    First tries to load a user-defined config file. If it doesn't find one,
    it looks for a default config file. If it fails again it means this is
    the first time zhunter is being used and it creates the default config
    file and saves it in '~/.config/zhunter/default_config.yaml'

    Returns
    -------
    Path
        `Path` object of the configuration file name.
    """

    log.debug("Looking for config file")
    # Try loading a user-defined config file
    config_fname = Path("~/.config/zhunter/user_config.yaml").expanduser()
    # If user-defined config file doesn't exist, try loading default config file
    if not config_fname.exists():
        config_fname = Path("~/.config/zhunter/default_config.yaml").expanduser()
        # if the default config file doesn't exist, it means it is the first
        # time zhunter is installed on the computer, so copy the file
        if not config_fname.exists():
            log.debug(
                "No configuration under '~/.config/zhunter/default_config.yaml'"
                " on this computer, creating one"
            )
            config_fname.parent.mkdir(parents=True, exist_ok=True)
            default_config_fname = DIRS["CONFIG"] / "default_config.yaml"
            shutil.copyfile(default_config_fname, config_fname)

    log.debug(f"Config file found:\n{config_fname}")

    return config_fname


def define_paths(input_fnames, default=True):
    """Define the paths to the various files used throughout the code.
    The fnames is used to get the name of certain paths.
    If default is True, will look in './data/lines/' for the line lists.
    Otherwise the full path to the file is needed.

    Parameters
    ----------
    fnames : dict
        Dictionary containing the files names with the keys
        ['intervening_lines', 'emission_lines', 'GRB_lines']
    default : bool, optional
        Where to search for the files. If `True`, will search in './data/lines/',
        if `False`, the full path must be provided for the files in the config.

    Returns
    -------
    dict
        Dictionary containing the file names.
    """
    fnames = {}
    for line_type in LINE_TYPES:
        if default:
            fnames[line_type+'_lines'] = line_dir / input_fnames[line_type+'_lines']
        else:
            fnames[line_type+'_lines'] = Path(input_fnames[line_type+'_lines'])

    fnames["line_ratio"] = DIRS["DATA"] / "lines/line_ratio.csv"

    fnames["tellurics"] = (
        DIRS["DATA"] / "tellurics/sky_transimission_opt_to_nir.ecsv.gz"
    )

    fnames["sky_bkg"] = (
        DIRS["DATA"] / "sky_background/sky_background_norm_opt_to_nir.ecsv.gz"
    )

    # Make sure that all files are well defined
    for f in fnames.values():
        if not f.exists():
            raise FileNotFoundError(f"File '{f}' does not exist.")

    fnames['config'] = input_fnames['config']

    # Get last opened file
    if Path(input_fnames['last_opened']).expanduser().exists():
        fnames["last_opened"] = input_fnames['last_opened']
    else:
        fnames["last_opened"] = Path('~').expanduser()

    return fnames


def load_colors(style="kraken17", default="kraken17"):
    """Load colors into dictionary.

    Parameters
    ----------
    style : str
        Name of the style to load. Available styles are :
        ["kraken9", "kraken17", "cvd", "old", "cyberpunk"]
    Returns
    -------
    colors : dict
        Dictionary containing the name of various items
        to be plotted and their associated color.
    """

    # Define color style
    try:
        colors = COLORS[style]
    except KeyError:
        log.error(
            "This color palette doest not exist. Please use one of "
            f"{list(COLORS.keys())}. Falling back to default colors: '{default}'"
        )
        colors = COLORS[default]
    return colors


def load_config(fname):
    """Load a yaml configuration file

    Returns
    -------
    dict
        Configuration loaded from the yaml file.
    """

    # Load config file
    with open(fname, "r") as f:
        config = yaml.safe_load(f)
        log.info(f"Loaded configuration from:\n{fname}")

    config['fnames']['config'] = fname
    return config


def load_line_lists(fnames, calc_ratio=True):
    """
    Load the input line lists and check the format is ok.

    Parameters
    ----------
    fnames : dict
        Dictionary containing the names of the files for
        each type of line list. Expected keys are ['intervening_lines',
        'emission_lines', 'GRB_lines']
    calc_ratio : bool, optional, default=True
        If True, will recalculate the line ratios using the list of lines
        of `intervening`. If False, will try to load a file name
        from the `fnames` dictionary using the 'line_ratio' key.

    Returns
    -------
    lines : dict
        Dictionary containing the line lists loaded in memory.
        The keys are ['intervening', 'emission', 'GRB',
        'line_ratio']

    """
    lines = {}

    for line_type in LINE_TYPES:
        log.debug(
            f"Reading {line_type} lines from:\n{fnames[line_type+'_lines']}"
        )
        lines[line_type] = io.read_line_list(fnames[line_type+'_lines'])

    if calc_ratio:
        ratios = create_line_ratios(fnames["intervening_lines"])
    else:
        log.debug(f"Read line ratios fron:\n{fnames['line_ratio']}")
        ratios = ascii_read(fnames["line_ratio"])

    lines['ratios'] = ratios

    return lines


def update_last_opened(fnames):
    # Get the path to the configuration file
    config_fname = fnames['config']
    last_opened = fnames['last_opened']

    with open(str(config_fname), "r") as f:
        config = yaml.safe_load(f)

    config['fnames']['last_opened'] = str(last_opened)

    with open(str(config_fname), "w") as f:
        yaml.dump(config, stream=f)
        log.debug(f"Updated last opened file with:\n{last_opened}")
