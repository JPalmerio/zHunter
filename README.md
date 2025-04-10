[![DOI](https://zenodo.org/badge/465226687.svg)](https://doi.org/10.5281/zenodo.15189495)

# *z*Hunter

*z*Hunter is a Graphical User Interface (GUI) tool to visualize and perform basic manipulation of 1D and 2D astronomical spectra.
It is originally developed to help find (hunt for) the redshift *z* of transient sources observed spectroscopically, hence its name.


# Installation

If you use a virtual environment manager for python (recommended!), you can create an environment specific for *z*Hunter with:

```
conda create -n zhunter python=3.10
```
Then don't forget to activate it when using zhunter with:
```
conda activate zhunter
```

## Using pip

```
$ pip install zhunter
```

If you want the latest development you can clone the project, move to the root of the project, switch to the `dev` branch and do:
```
$ pip install -e .
```

## Launching the GUI
If the installation went smoothly, you can launch the GUI by simply typing in your terminal:
```
$ zhunter
```
