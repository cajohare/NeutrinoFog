[![arXiv](https://img.shields.io/badge/arXiv-2105.04565-B31B1B.svg)](https://arxiv.org/abs/2109.XXXX)
[![MIT Licence](https://badges.frapsoft.com/os/mit/mit.svg?v=103)](https://opensource.org/licenses/mit-license.php)


# NeutrinoFog
Python-3 Code to reproduce the results from my paper arXiv:[2109.XXXX]

The heavy lifting of the likelihood analysis is done in a simple fortran program (src/like/like.f95), which is wrapped in python and can be run by executing one of the python scripts:
* src/runNuFloor_SI
* src/runNuFloor_SDp
* src/runNuFloor_SDn
* src/runNuFloor_NuFluxUncertainties

If you need any assistance or have any questions, contact me at ciaran.aj.ohare@gmail.com

If you use anything here please cite the paper, [O'Hare 2021](https://arxiv.org/abs/2109.?????)
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

To start, compile the fortran likelihood code by doing
```
cd src/like/
make all
```

Then to run one of the python scripts do e.g.
```
cd ..
python runNuFloor_SI.py Xe
```

# Requirements
* [`CMASHER`](...)
* [`numba`](http://numba.pydata.org/)
* [`tqdm`](https://pypi.org/project/tqdm/)

# Contents
* [`src/`](https://github.com/cajohare/DarkPhotonCookbook/tree/main/code) contains python and fortran code to produce the main results
* [`plots/`](https://github.com/cajohare/DarkPhotonCookbook/tree/main/plots) contains all the plots in pdf and png formats
* [` notebooks/`](https://github.com/cajohare/DarkPhotonCookbook/tree/main/notebooks) jupyter notebooks for plotting results

# Python modules
* [`LabFuncs.py`](https://github.com/cajohare/DarkPhotonCookbook/blob/master/code/LabFuncs.py)
* [`PlotFuncs.py`](https://github.com/cajohare/DarkPhotonCookbook/blob/master/code/PlotFuncs.py)
* [`WIMPFuncs.py`](https://github.com/cajohare/DarkPhotonCookbook/blob/master/code/WIMPFuncs.py)
* [`NeutrinoFuncs.py`](https://github.com/cajohare/DarkPhotonCookbook/blob/master/code/NeutrinoFuncs.py)
* [`Like.py`](https://github.com/cajohare/DarkPhotonCookbook/blob/master/code/Like.py)

# Notebooks
* [`Main.ipynb`](https://github.com/cajohare/DarkPhotonCookbook/blob/master/code/Main.ipynb)
* [`Explanation.ipynb`](https://github.com/cajohare/DarkPhotonCookbook/blob/master/code/Explanation.ipynb)
* [`SI.ipynb`](https://github.com/cajohare/DarkPhotonCookbook/blob/master/code/SI.ipynb)
* [`SDp.ipynb`](https://github.com/cajohare/DarkPhotonCookbook/blob/master/code/SDp.ipynb)
* [`SDn.ipynb`](https://github.com/cajohare/DarkPhotonCookbook/blob/master/code/SDn.ipynb)
* [`Write_NuFloors.ipynb`](https://github.com/cajohare/DarkPhotonCookbook/blob/master/code/Write_NuFloors.ipynb)

# Example result
[<img align="right" src="plots/plots_png/NuFloorExplanation.png" height="250">](https://github.com/cajohare/DarkPhotonCookbook/raw/master/plots/plots_png/NuFloorExplanation.png)
### [Latitude dependence:](https://github.com/cajohare/DarkPhotonCookbook/blob/master/code/Explanation.ipynb)
This plot is to show
### &nbsp;
### &nbsp;
### &nbsp;
### &nbsp;
---

# If you just want txt files for the neutrino floors themselves here they are:

## Spin independent
* Xenon: [.txt](https://github.com/cajohare/AxionLimits/raw/master/data/floors/NeutrinoFloor_Xe_SI.txt)
* Germanium: [.txt](https://github.com/cajohare/AxionLimits/raw/master/data/floors/NeutrinoFloor_Ge_SI.txt)
* Argon: [.txt](https://github.com/cajohare/AxionLimits/raw/master/data/floors/NeutrinoFloor_Ar_SI.txt)
* Helium: [.txt](https://github.com/cajohare/AxionLimits/raw/master/data/floors/NeutrinoFloor_He_SI.txt)
* Fluorine: [.txt](https://github.com/cajohare/AxionLimits/raw/master/data/floors/NeutrinoFloor_F_SI.txt)
* CaWO4: [.txt](https://github.com/cajohare/AxionLimits/raw/master/data/floors/NeutrinoFloor_CaWO4_SI.txt)
* NaI: [.txt](https://github.com/cajohare/AxionLimits/raw/master/data/floors/NeutrinoFloor_NaI_SI.txt)

## Spin dependent (proton)
* Xenon: [.txt](https://github.com/cajohare/AxionLimits/raw/master/data/floors/NeutrinoFloor_Xe_SDp.txt)
* Germanium: [.txt](https://github.com/cajohare/AxionLimits/raw/master/data/floors/NeutrinoFloor_Ge_SDp.txt)
* Silicon: [.txt](https://github.com/cajohare/AxionLimits/raw/master/data/floors/NeutrinoFloor_Si_SDp.txt)
* Fluorine: [.txt](https://github.com/cajohare/AxionLimits/raw/master/data/floors/NeutrinoFloor_F_SDp.txt)
* NaI: [.txt](https://github.com/cajohare/AxionLimits/raw/master/data/floors/NeutrinoFloor_NaI_SDp.txt)

## Spin dependent (neutron)
* Xenon: [.txt](https://github.com/cajohare/AxionLimits/raw/master/data/floors/NeutrinoFloor_Xe_SDn.txt)
* Germanium: [.txt](https://github.com/cajohare/AxionLimits/raw/master/data/floors/NeutrinoFloor_Ge_SDn.txt)
* Silicon: [.txt](https://github.com/cajohare/AxionLimits/raw/master/data/floors/NeutrinoFloor_Si_SDn.txt)
* Fluorine: [.txt](https://github.com/cajohare/AxionLimits/raw/master/data/floors/NeutrinoFloor_F_SDn.txt)
* NaI: [.txt](https://github.com/cajohare/AxionLimits/raw/master/data/floors/NeutrinoFloor_NaI_SDn.txt)
