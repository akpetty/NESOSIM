# NASA's Eulerian Snow on Sea Ice Model (NESOSIM)
**Alek Petty**

The NASA Eulerian Snow On Sea Ice Model (NESOSIM) is a new open source model that produces daily estimates of the depth and density of snow on sea ice across the polar oceans.

![NESOSIM schematic](schamtic.jpg?raw=true "NESOSIM V1 schematic")

NESOSIM has been developed in a three-dimensional Eulerian framework and includes several parameterizations that represent key mechanisms of snow variability through the snow accumulation/growth season. NESOSIM currently has two vertical snow layers to broadly represent the evolution of both old/compacted snow and new/fresh snow. 

NESOSIM is being made available as an open source project to encourage continued model development and active engagement with the snow on sea ice community. The model code is written in Python, an open source programming language (Python Software Foundation, https://www.python.org/), to better enable future community development efforts. Our hope is that the model will continue to evolve as additional snow processes are incorporated, especially as new field and remote sensing snow observations are collected and made available for calibration/validation. Obvious examples of planned future improvements include the incorporation of snow melt and rain on snow processes, which are not currently included in this initial model version, enabling the model to be run year-round.

For more details of the model physics and preliminary results/calibration efforts see the following discussion paper:

Petty, A. A., M. Webster, L. N. Boisvert, T. Markus, The NASA Eulerian Snow on Sea Ice Model (NESOSIM): Initial model development and analysis, Geosci. Mod. Dev.


## Getting Started

### Conda installation

We used conda to install and maintain our Python (2.7) environment. The relevant dependancies of this conda Python environment can be found in nesosim.yml, providing a direct means of cloning this Python environment directly (including the same module versions). 

The code was also tested on a remote server using a new default conda 2.7 Python enviornment and the following additions:

```
conda create -n <ENVNAME> python=2.7
conda install scipy
conda install matplotlib
conda install basemap
conda install h5py
conda install netCDF4
```
This conda Python environment can then be activated with 

```
source activate <ENVNAME>
```

Alternatively you can try running NESOSIM in your own Python environment and installing the libraries you may be missing, although there will likely be differences in the module versions used. I would appreciate this feedback, especially from those trying to run NESOSIM in Python 3.

Further information about installing Conda/Python, and a brief introduction to using Python can be found on my NASA Cryospheric Sciences meetup repo: https://github.com/akpetty/cryoscripts.

### Model code

The primary model code can be found in Main/

Descriptions should be included at the top of each Python script. 

The model can be run by calling 

```
python run.py
```
which calls the main model function in NESOSIM.py, including functions given in commonFuncs.py


### Forcing Data

NESOSIM is forced with daily inputs of snowfall and near-surface winds (from reanalyses), sea ice concentration (from satellite passive microwave data) and sea ice drift (from satellite feature tracking), during the accumulation season (August through April).  
The various forcing data used to run NESOSIM are described in Petty et al., (2018, GMD).

The TestData folder contains gridded (100 km) test forcing data for 2010-2011: ERA-I (snowfall and winds) and MEDIAN-SF (snowfall), Bootstrap (ice concentration), NSIDCv3 (ice drift). 

All the gridded data used in Petty et al., (2018, GMD) will be made available to download shortly, which should be copied to a local Data folder (unless you change the DataPath variable in NESOSIM.py)

### Output

Add info here.

### Plotting Scripts

Add info here once these scripts are uploaded.

Contact me if you any any questions (alek.a.petty@nasa.gov)!

Alek



