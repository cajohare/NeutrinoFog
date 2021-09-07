[![arXiv](https://img.shields.io/badge/arXiv-2105.04565-B31B1B.svg)](https://arxiv.org/abs/2109.XXXX)
[![MIT Licence](https://badges.frapsoft.com/os/mit/mit.svg?v=103)](https://opensource.org/licenses/mit-license.php)

# NeutrinoFog
Code to reproduce the results from my paper arXiv:[2109.XXXX]

[<img align="right" src="plots/plots_png/NuFloorExplanation.png" height="350">](https://github.com/cajohare/NeutrinoFog/raw/master/plots/plots_png/NuFloorExplanation.png)

---

# Requirements
* [`CMasher`](https://cmasher.readthedocs.io/)
* gfortran (or equivalent f95 compiler, but make sure to edit src/like/Makefile)

# Running the code
The heavy lifting of the likelihood analysis is done in a short fortran program (src/like/like.f95), which can be executed via python. To start, compile the fortran likelihood code by doing
```
make -C src/like
```
Then to run one of the python scripts do e.g.
```
python src/runNuFloor_SI.py Xe
```
using one of the following scrips:
* [`src/runNuFloor_SI.py`](https://github.com/cajohare/NeutrinoFog/tree/main/src/runNuFloor_SI.py) - SI limits
* [`src/runNuFloor_SDp.py`](https://github.com/cajohare/NeutrinoFog/tree/main/src/runNuFloor_SDp.py) - SD (proton) limits
* [`src/runNuFloor_SDn.py`](https://github.com/cajohare/NeutrinoFog/tree/main/src/runNuFloor_SDn.py) - SD (neutron) limits
* [`src/runNuFloor_NuFluxUncertainties.py`](https://github.com/cajohare/NeutrinoFog/tree/main/src/runNuFloor_NuFluxUncertainties.py) - To test the effect of neutrino flux uncertainties
Each takes about 30 minutes per target on my machine, but the resolutions can be lowered to generate them faster.

# Contents
* [`src/`](https://github.com/cajohare/NeutrinoFog/tree/main/src) - Contains python and fortran code to produce the main results
* [`plots/`](https://github.com/cajohare/NeutrinoFog/tree/main/plots) - Contains all the plots in pdf and png formats
* [`notebooks/`](https://github.com/cajohare/NeutrinoFog/tree/main/notebooks) - Jupyter notebooks for plotting results
* [`data/`](https://github.com/cajohare/NeutrinoFog/tree/main/notebooks) - data files, including neutrino fluxes, experimental exclusion limits, and the data produced by this code

# Python modules
* [`WIMPFuncs.py`](https://github.com/cajohare/NeutrinoFog/blob/master/code/WIMPFuncs.py) - Functions needed to calculate WIMP rates
* [`NeutrinoFuncs.py`](https://github.com/cajohare/NeutrinoFog/blob/master/code/NeutrinoFuncs.py) - Functions needed to calculate neutrino rates
* [`Like.py`](https://github.com/cajohare/NeutrinoFog/blob/master/code/Like.py) - Functions for wrapping and running the likelihood code
* [`LabFuncs.py`](https://github.com/cajohare/NeutrinoFog/blob/master/code/LabFuncs.py) - Various utilities
* [`PlotFuncs.py`](https://github.com/cajohare/NeutrinoFog/blob/master/code/PlotFuncs.py) - Plotting functions

# Notebooks
* [`Main.ipynb`](https://github.com/cajohare/NeutrinoFog/blob/master/code/Main.ipynb) - Figures 1 and 3b of the paper
* [`Explanation.ipynb`](https://github.com/cajohare/NeutrinoFog/blob/master/code/Explanation.ipynb) - Figure 2
* [`SI.ipynb`](https://github.com/cajohare/NeutrinoFog/blob/master/code/SI.ipynb) - Figure 3a and S1
* [`SDp.ipynb`](https://github.com/cajohare/NeutrinoFog/blob/master/code/SDp.ipynb) - Figure S2a
* [`SDn.ipynb`](https://github.com/cajohare/NeutrinoFog/blob/master/code/SDn.ipynb) - Figure S2b

---

# If you just want .txt files for the neutrino floors themselves here they are:
## Spin independent
* Xenon: [NeutrinoFloor_Xe_SI.txt](https://github.com/cajohare/NeutrinoFog/raw/master/data/floors/NeutrinoFloor_Xe_SI.txt)
* Germanium: [NeutrinoFloor_Ge_SI.txt](https://github.com/cajohare/NeutrinoFog/raw/master/data/floors/NeutrinoFloor_Ge_SI.txt)
* Argon: [NeutrinoFloor_Ar_SI.txt](https://github.com/cajohare/NeutrinoFog/raw/master/data/floors/NeutrinoFloor_Ar_SI.txt)
* Helium: [NeutrinoFloor_He_SI.txt](https://github.com/cajohare/NeutrinoFog/raw/master/data/floors/NeutrinoFloor_He_SI.txt)
* Fluorine: [NeutrinoFloor_F_SI.txt](https://github.com/cajohare/NeutrinoFog/raw/master/data/floors/NeutrinoFloor_F_SI.txt)
* CaWO4: [NeutrinoFloor_CaWO4_SI.txt](https://github.com/cajohare/NeutrinoFog/raw/master/data/floors/NeutrinoFloor_CaWO4_SI.txt)
* NaI: [NeutrinoFloor_NaI_SI.txt](https://github.com/cajohare/NeutrinoFog/raw/master/data/floors/NeutrinoFloor_NaI_SI.txt)

## Spin dependent (proton)
* Xenon: [NeutrinoFloor_Xe_SDp.txt](https://github.com/cajohare/NeutrinoFog/raw/master/data/floors/NeutrinoFloor_Xe_SDp.txt)
* Germanium: [NeutrinoFloor_Ge_SDp.txt](https://github.com/cajohare/NeutrinoFog/raw/master/data/floors/NeutrinoFloor_Ge_SDp.txt)
* Silicon: [NeutrinoFloor_Si_SDp.txt](https://github.com/cajohare/NeutrinoFog/raw/master/data/floors/NeutrinoFloor_Si_SDp.txt)
* Fluorine: [NeutrinoFloor_F_SDp.txt](https://github.com/cajohare/NeutrinoFog/raw/master/data/floors/NeutrinoFloor_F_SDp.txt)
* NaI: [NeutrinoFloor_NaI_SDp.txt](https://github.com/cajohare/NeutrinoFog/raw/master/data/floors/NeutrinoFloor_NaI_SDp.txt)

## Spin dependent (neutron)
* Xenon: [NeutrinoFloor_Xe_SDn.txt](https://github.com/cajohare/NeutrinoFog/raw/master/data/floors/NeutrinoFloor_Xe_SDn.txt)
* Germanium: [NeutrinoFloor_Ge_SDn.txt](https://github.com/cajohare/NeutrinoFog/raw/master/data/floors/NeutrinoFloor_Ge_SDn.txt)
* Silicon: [NeutrinoFloor_Si_SDn.txt](https://github.com/cajohare/NeutrinoFog/raw/master/data/floors/NeutrinoFloor_Si_SDn.txt)
* Fluorine: [NeutrinoFloor_F_SDn.txt](https://github.com/cajohare/NeutrinoFog/raw/master/data/floors/NeutrinoFloor_F_SDn.txt)
* NaI: [NeutrinoFloor_NaI_SDn.txt](https://github.com/cajohare/NeutrinoFog/raw/master/data/floors/NeutrinoFloor_NaI_SDn.txt)

---

If you need any further assistance or have any questions, contact me at ciaran.aj.ohare@gmail.com. And if you do use anything here please cite the paper, [O'Hare 2021](https://arxiv.org/abs/2109.?????)
```
@article{OHare:2021xxx,
    author = "O'Hare, Ciaran A. J.",
    title = "{Fog on the horizon: a new definition of the neutrino floor for dark matter searches}",
    eprint = "2109.XXXXX",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    month = "9",
    year = "2021"
}
```
