import shutil
import logging
import yaml

from pathlib import Path
from zhunter import DIRS

log = logging.getLogger(__name__)


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
            log.debug("No configuration on this computer, creating one")
            config_fname.parent.mkdir(parents=True, exist_ok=True)
            default_config_fname = DIRS["CONFIG"] / "default_config.yaml"
            shutil.copyfile(default_config_fname, config_fname)

    log.debug(f"Config file found:\n{config_fname}")

    return config_fname


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

    return config
