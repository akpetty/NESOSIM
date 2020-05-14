# The NASA Eulerian Snow on Sea Ice Model (NESOSIM) v1.0
**Contact: Alek Petty / alek.a.petty@nasa.gov / www.alekpetty.com**

The NASA Eulerian Snow On Sea Ice Model (NESOSIM) is a three-dimensional, two-layer (vertical), Eulerian snow on sea ice budget model developed with the primary aim of producing daily estimates of the sea ice snow depth and snow density across the polar oceans.  

NESOSIM v1.0 includes several parameterizations that represent key mechanisms of snow variability through the snow accumulation/growth season, and two snow layers to broadly represent the evolution of both old/compacted snow and new/fresh snow. 


![NESOSIM schematic](schematic.jpg?raw=true "NESOSIM v1.0 schematic")


NESOSIM is being made available as an open source project to encourage continued model development and active engagement with the snow on sea ice community. The model code is written in Python, an open source programming language (Python Software Foundation, https://www.python.org/), to better enable future community development efforts. Our hope is that the model will continue to evolve as additional snow processes are incorporated, especially as new field and remote sensing snow observations are collected and made available for calibration/validation. Obvious examples of planned future improvements include the incorporation of snow melt and rain on snow processes, which are not currently included in this initial model version, enabling the model to be run year-round.

For more details of the model physics and preliminary results/calibration efforts see the following discussion paper:

Petty, A. A., M. Webster, L. N. Boisvert, T. Markus, The NASA Eulerian Snow on Sea Ice Model (NESOSIM): Initial model development and analysis, Geosci. Mod. Dev.

Versions:
 - v1.0: This initial NESOSIM model version is configured to run only for the Arctic Ocean through the accumulation season (August 15th to May 1st).
 - v1.1: This new version of NESOSIM includes a few small updates including: upgrade to Python 3, a bigger domain, introduction of CloudSat scaling parameters (Cabaj et al., 2020)

## Getting Started

### Conda installation

I recommend using the included conda environment file - nesosim.yml - to ensure consistency in the Python 2.7 configuration when running this code. Note that you might struggle loading the pickled data in Python 3, hence the delay in updating (need to then regrid all the forcing data).

```
conda env create -f nesosim27.yml
```

Alternatively you can try generating your own conda environment using the following packages

```
conda create -n nesosim27 python=2.7 scipy matplotlib basemap h5py netCDF4 xarray proj4

```
The conda Python environment can be activated with 

```
source activate nesosim27
```

Or if you really want you can try simply installing the libraries it flags may be missing when you try and run the scripts, although there will likely be differences in the module versions used. 

Further information about installing Conda/Python, and a brief introduction to using Python can be found on my NASA Cryospheric Sciences meetup repo: https://github.com/akpetty/cryoscripts.

### Model code

The NESOSIM model source code can be found in 

```
Scripts/NESOSIM.py
```
which also needs various functions included in commonFuncs.py. This file is best run using the seperate configuration/run script.

```
python run.py
```

Also included in this repo:
```
Scripts/Analysis/
```
- Analysis scripts used in Petty et al., (2018, GMD).

```
Scripts/GetData/
```
 - Scripts to download the raw forcing data files.

```
Scripts/Gridding/
```
- Scripts to grid these raw forcing data to the model domain.

```
Scripts/Plotting
```
- Plotting scripts used in Petty et al., (2018, GMD).

```
saveBudgetsNCDF.py
```
- Convert the xarray model output to netCDF files (both stored in Output).

Descriptions should be included at the top of each Python script. 


### Forcing Data


NESOSIM (v1.0) is forced with daily inputs of snowfall and near-surface winds (from reanalyses), sea ice concentration (from satellite passive microwave data) and sea ice drift (from satellite feature tracking), during the accumulation season (August through April).  

The various forcing data used to run NESOSIM are described in Petty et al., (2018, GMD).

```
Forcings/
```
 - Currently empty because of size constraints, but all the gridded forcing data used in Petty et al., (2018, GMD) have been made available on the NASA Cryospheric Sciences Lab website: https://neptune.gsfc.nasa.gov/csb/index.php?section=516, which can be downloaded, unzipped, and copied to a the Forcings folder (unless you change the DataPath variable in NESOSIM.py). Note that the files in Scripts/gridding shows how these forcings were generated from the raw data.


```
TestForcings/
```

 - Includes sample gridded (100 km) test forcing data for 2010-2011: ERA-I (snowfall and winds) and MEDIAN-SF (snowfall), Bootstrap (ice concentration), NSIDCv3 (ice drift). 


### Ancillary Data

```
AncData/
```
- Contains temperature-scaled intitial snow depths and the drifting soviet station data that were used for model calibration.

### Model Output

 NetCDF files used by the various analysis/plotting scripts included in this repo, which include all the model variables, can be found out:
```
Output/budgets
```

NetCDF files only including the primary model variables are also generated and stored in:
```
Output/final
```
A single NetCDF file is provided for each annual accumulation season model run (August 15th to May 1st of the following year), providing daily data on the 100 km polar stereographic model domain. 
In this final output, each NetCDF file includes the following variables: 
 - Effective snow depth (snowDepth)
 - Grid-cell snow volume (snowVol)
 - Bulk snow density (density)
 - Reanalysis-derived daily cumulative snowfall (Precip)
 - Passive microwave ice concentration data (iceConc)
 - Day of the year (day). 

 The final NetCDF data files for the ERA-I and MEDIAN forced simulations are also hosted on the NASA Cryospheric Sciences Lab website: https://neptune.gsfc.nasa.gov/csb/index.php?section=516


Contact me if you any any questions (alek.a.petty@nasa.gov)!

Alek



