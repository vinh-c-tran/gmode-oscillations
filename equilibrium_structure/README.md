# Equilibrium Structure 
This first folder/set of code solves for the equilibrium structure of the star. In this case, this means determining first the various particle fractions and meson field values as a function of baryon number density. All of the functions are stored in the various `.py` files and the usage is shown in the `Solver Example.ipynb` file. It is meant to be called in a Jupyter notebook for easy visualization. 

# Various Files
## `nlw_classes.py`
This file defines the `eos` class which stores all of the parameters for a given model or EOS such as the meson-meson and meson-baryon coupling constants. Then a given `eos` is initialized as an `eos` object or an instance of this `eos` class. See the `nlw_eos.py` file which gives examples of this. Also defined in this file are the `particle` class and children classes `baryon`, `lepton`, and `meson` which store attributes of a given particle such as mass and charge and functions which calculate things like their energy or momentum. At the bottom various particles are instantiated such as `Neutron` for example. 

## `nlw_eos.py` 
This file has the various models instantiated as `eos` objects for use in this work. 

## `solver_func.py`
This file sets up the system of nonlinear equations that we have to solve to determine the equilibrium structure. See the paper or my thesis for more details on what these equations exactly are. 
