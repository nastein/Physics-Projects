# DARKCAST

Darkcast is the companion software package to the paper [*Serendipity in dark photon searches*](https://arxiv.org/abs/1801.04847) and is a framework for recasting constraints from dark photon searches into other models. The Darkcast package is written as a module in Python and has no external dependencies, sans Python itself. To begin recasting, download the source and try running some of the examples:
```bash
wget https://gitlab.com/philten/darkcast/-/archive/master/darkcast-master.tar.gz
tar -xzvf darkcast-master.tar.gz
mv darkcast-master darkcast
cd darkcast/examples
```

There are a number of examples provided corresponding to the figures of [*Serendipity in dark photon searches*](https://arxiv.org/abs/1801.04847). Each can be run as:
```bash
python <example>.py
```
where `<example>.py` is the relevant example. If the Python plotting module `matplotlib` is available, all the examples will produce plots in the PDF format. Otherwise, only text output of the recasting will be produced, which can be read in by the user's plotting utility of choice.

1. `bfrac.py`: this example calculates the $`\mu\mu`$, $`ee`$, $`\nu\nu`$ and hadronic branching fractions for all the available models in Darkcast. This will produce the output `bfrac_<model>_<channel>.txt` where each file gives the branching fraction for the given model and channel. The first column is the $`X`$ boson mass in GeV and the second column is the branching fraction. The plots `bfrac_<model>.pdf` will also be produced.
2. `visible.py`: loads all the available visible limits and recasts them to all the available models. Depending on the machine, this example can take some time. To avoid clutter a `recast` directory is created in the directory the example is run from and text file are produced in the format `recast/<model>/<limit>.lmt`. For each limit the first column is the $`X`$ boson mass in GeV, the second column is the lower global coupling limit, and if the limit is two-sided an optional third column gives the upper global coupling limit.
3. `invisible.py`: is the same as `visible.py` but now recasts all the invisible limits. The text file limits are written to the same directory structure as `visible.py`
4. `user.py`: demonstrates how new user defined models and limits can be written, and outputs the limits `recast/<model>/user_limit.lmt` for the user limit defined by `user_limit.py` using the models `dark_photon` and `user_model`. Here, `user_model` is defined by `user_model.py`. Optionally, the plot `user.pdf` is produced. The following files are used by `user_limit.py`, where details on how each of these files can be used are given in `user_limit.py`.
  * `user_limit_single.lmt`: defines an example lower bound limit.
  * `user_limit_double.lmt`: defines an example double-sided limit.
  * `user_limit_rvalue.lmt`: defines an example full limit using r-values.
  * `user_limit.prd`: defines the production mechanisms for the limit.
5. `logo.py`: draws the Darkcast logo.

The following is a simple usage example which recasts the prompt LHCb dark photon limits to the $`B`$ boson model.
```python
# Load the module.
import darkcast

# Change any global parameters, here the speed of light (m/s).
darkcast.pars.c = 3e8

# Load a limit, here the LHCb prompt limit.
limit = darkcast.Limit('LHCb_Aaij2017rft_prompt')

# Print the notes and BibTex for the limit.
print limit.notes
print limit.bibtex

# Load a model for recasting.
model = darkcast.Model('B_boson')

# Recast from the limit model to the new model.
recast = limit.recast(model)

# Write out the recast limit.
recast.write('LHCb_B_boson.lmt)
```

## References

When using Darkcast, please cite [*Serendipity in dark photon searches*](https://arxiv.org/abs/1801.04847) as published in JHEP. Individual citations are also provided for each limit via `limit.bibtex` and a comprehensive list of references is provided in the [refs](refs) directory.

## Licensing

Darkcast is licensed under the GNU GPL version 2, or later and is copyrighted (C) 2020 by Philip Ilten, Yotam Soreq, Mike Williams, and Wei Xue.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the [GNU General Public License](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html) or [LICENSE](LICENSE) for more details.

## Structure

The structure of Darkcast is as follows, beginning at the top level.

* `README.md`: is this file.
* `__init__.py`: initializes the Darkcast package.
* `efficiency.py`: defines the `Efficiency` class used to calculate efficiency ratios.
* `limit.py`: defines the classes used to create a limit.
* `model.py`: defines the classes needed to create model.
* `pars.py`: contains all the parameters used by Darkcast.
* `production.py`: defines the needed classes for calculating production ratios.
* `utils.py`: contains auxiliary utilities which are not physics related.

### [limits](limits)

All available limits distributed with Darkcast are provided here. Each limit is defined by the file `<experiment>_<INSPIRE tag>_<optional information>.py`. For each limit, the corresponding bounds are provided by the `*.lmt` files. Some limits require complex production mechanism fractions to be defined; this is done with the `*.prd` files.

When limits are loaded in Darkcast via:
```python
import darkcast
limit = darkcast.Limit('name')
```
the following search paths are used.
1. The current directory within the Python interpreter.
2. The paths defined by the environment variable `DARKCAST_LIMIT_PATH`.
3. The `limits` directory of the Darkcast package.

Any data required by the limits is searched along the following paths:
  1. The absolute path, if the absolute path is given.
  2. The current directory within the Python interpreter.
  3. The paths defined by the environment variable `DARKCAST_DATA_PATH`.
  4. The Darkcast package directory.

### [models](models)

All available models distributed with Darkcast are provided here. Each model is defined by the file `<model>.py`. Currently, dark photon, $`B`$ boson, $`B-L`$ boson, and protophobic models are defined. When models are loaded in Darkcast via:
```python
import darkcast
model = darkcast.Model('name')
```
the following search paths are used.
1. The current directory within the Python interpreter.
2. The paths defined by the environment variable `DARKCAST_MODEL_PATH`.
3. The `models` directory of the Darkcast package.

### [reach](reach)

A relatively comprehensive set of projected future limits are provided. While an attempt has been made to define the relevant efficiencies and production mechanisms for these projections, no validation has been performed, and recasting these projections may give nonsensical results.

### [recast](recast)

The recasted bounds for all valid limits in `limits` and `reach` are provided here in `recast/limits` and `recast/reach` respectively. A sub-directory is provided for each model defined in `models`. These recast limits are generated by running `python examples/visible.py`, `python examples/invisible.py`, and `python examples/reach.py`. For reference details on each recast limit, please consult the `bibtex` member of the corresponding limit:
```python
import darkcast
limit = darkcast.Limit('BaBar_Lees2014xha') 
print limit.bibtex
```
Note that the recasted reaches have not been fully validated.

### [refs](refs)

This directory contains a `.tex` and `.bib` file that provide easy access to
dark-sector references.

### [vmd](vmd)

This directory contains all the VMD data needed to calculate limits. Twelve curves are provided, giving the $`\mathcal{R}_\mu^\mathcal{F}(m)`$ curves of figure 1 from [*Serendipity in dark photon searches*](https://arxiv.org/abs/1801.04847). Here, $`\mathcal{F}`$ specifies the final state. Each final state is divided into its $`\omega`$, $`\phi`$, and $`\rho^0`$ components, including the $`\omega-\phi`$ interference term. All curves are given as a function of $`m`$ in GeV.

## LXPLUS (ATLAS specific) Directions

Depending on the default environmnent for a user on LXPLUS, it is possible that Darkcast will not work out of the box. Specifically, LHCb and CMS users appear to have no issue, while ATLAS users do. This typically is because the Python version is too old (less than `2.7`). Additionally, `matplotlib` is not necessarily available by default and so plotting function may not work out of the box.

During a recent Darkcast tutorial, Caterina Doglioni suggested the following prescription for setting up a working Darkcast environment for ATLAS users via a `conda` environment.

```bash
# Requires a graphical session (conda will look for DISPLAY), e.g. ssh -Y.
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
chmod 755 Miniconda2-latest-Linux-x86_64.sh
# When asked to put the miniconda directory in your PATH say no.
export PATH=$PWD/miniconda2/bin:$PATH
conda create --name=darkcastenv python=2.7
source activate darkcastenv
conda install -c conda-forge matplotlib
# Every time you log in then run the following.
export PATH=/afs/cern.ch/user/d/doglioni/Work/miniconda2/bin:$PATH
source activate darkcastenv
```

## History

* **v0.2** - 09/04/2019
  * Second beta release of the Darkcast code. A number of limits have been added, including reach for future experiments.
* **v0.1** - 31/05/2018
  * First beta release of the full Darkcast code.