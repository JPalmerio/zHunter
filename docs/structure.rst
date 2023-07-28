Code organization
=================

zHunter is divided in modules, each responsible for its own part.
The main modules are:

    - `gui.py`: This module takes care of the interaction with the user.
    This is where the actions are coded, what button leads to what effect, etc.
    - `spectrum.py`: This module has two main classes : OneDSpectrum and TwoDSpectrum which
    contain all the relevant methods. For example each spectrum instance can be plotted,
    smoothed, using dedicated methods.
    - GraphicsWidget (composed of `OneDGraphicsWidget.py`, `MainGraphicsWidget.py`,
    `LineFitGraphicsWidget.py`, `VelocityGraphicsWidget.py`): These modules contain the
    various widget that will display and take care of plotting.


Configuration
=============
zHunter can be configure via its configuration file located in `~/.config/zhunter/`