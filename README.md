# *z*Hunter
*z*Hunter is a Graphical User Interface (GUI) tool to visualize and perform basic manipulation of 1D and 2D astronomical spectra.
It is originally developped to help find (hunt for) the redshift *z* of transient sources observed spectroscopically, hence its name.
It uses [`Python 3.8`](https://www.python.org/downloads/release/python-383/) and is based on the [`pyqtgraph`](https://pyqtgraph.readthedocs.io/en/latest/introduction.html) library for speed (as opposed to the more commonly used [`matplotlib`](https://matplotlib.org/)).


# Installation

If you use a virtual environment manager for python (which we recommend), you can create an environment specific for *z*Hunter with:

```
$ conda create -n zHunter python=3.8.3
$ pip install -r requirements.txt
```


There is no installation yet, just clone the repository and launch the code by running:

```
$ python zhunter.py
```

You can make sure the code works by loading the example file `example_2D.fits` (*hint*: GRB redshift is around 6.3).
