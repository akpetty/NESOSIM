# The NASA Eulerian Snow on Sea Ice Model (NESOSIM) v1.0
**Contact: Alek Petty / alek.a.petty@nasa.gov / www.alekpetty.com**

The NASA Eulerian Snow On Sea Ice Model (NESOSIM) is a three-dimensional, two-layer (vertical), Eulerian snow on sea ice budget model developed with the primary aim of producing daily estimates of the depth and density fo snow on sea ice across the polar oceans.  

NESOSIM includes several parameterizations that represent key mechanisms of snow variability through the snow accumulation/growth season, and two snow layers to broadly represent the evolution of both old/compacted snow and new/fresh snow. 

![NESOSIM schematic](schematic.jpg?raw=true "NESOSIM v1 schematic")

NESOSIM is being made available as an open source project to encourage continued model development and active engagement with the snow on sea ice community. The model code is written in Python, an open source programming language (Python Software Foundation, https://www.python.org/), to better enable future community development efforts. Our hope is that the model will continue to evolve as additional snow processes are incorporated, especially as new field and remote sensing snow observations are collected and made available for calibration/validation. Obvious examples of planned future improvements include the incorporation of snow melt and rain on snow processes, which are not currently included in this initial model version, enabling the model to be run year-round.

For more details of the model physics and preliminary results/calibration efforts see the following discussion paper:

Petty, A. A., M. Webster, L. N. Boisvert, T. Markus, The NASA Eulerian Snow on Sea Ice Model (NESOSIM): Initial model development and analysis, Geosci. Mod. Dev.

Versions:
 v1.0: This initial NESOSIM model version is configured to run only for the Arctic Ocean through the accumulation season (August 15th to May 1st). This was the version described in Petty et al., (2018) so please refer to that source code (click on the releases tab aobve) for that specific code version.
 v1.1: This latest version of NESOSIM includes a few minor updates. Thanks to Alex Cabaj for the help with some of this. Changes include: 
  - Upgrade of the code to Python 3.
  - 50 km (increased from 100 km) grid resolution.
  - An extended Arctic domain to cover all the peripheral seas.
  - Switched from Basemap to pyproj/cartopy (surprisingly painful).
  - Introduction of CloudSat scaling parameters by Alex Cabaj (Cabaj et al., 2020).

## Getting Started

### Conda installation

I recommend using the included conda environment file - nesosim3.yml - to ensure consistency in the Python 3.7 configuration when running this code. Note that in version 1.1 we have upgraded the code to run in Python 3 (3.7) so your earlier environment will likely not work anymore. The code changes were pretty small. I have found some issues with conda giving access to all the necessary libraries, so you might need to use pip to install things like xarray and cartopy

```
conda env create -f environment.yml
```

Alternatively you can try generating your own conda environment using the following packages

```
conda create -n nesosim3 python=3.7 scipy matplotlib basemap h5py netCDF4 xarray proj4 cdsapi

```
The conda Python environment can be activated with 

```
source activate nesosim3
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


NESOSIM (v1.1) is forced with daily (50 km) gridded inputs of snowfall and near-surface winds (from reanalyses), sea ice concentration (from satellite passive microwave data) and sea ice drift (from satellite feature tracking), during the accumulation season (August through April).  

The various forcing data used to run NESOSIM are described in Petty et al., (2018, GMD) but have ben updated to an extended 50 km (along with a few other tweaks) in this v1.1 configuration.

```
forcings/
```
 - Currently empty because of size constraints, but all the gridded forcing data used in Petty et al., (2018, GMD) have been made available on the NASA Cryospheric Sciences Lab website: https://neptune.gsfc.nasa.gov/csb/index.php?section=516, which can be downloaded, unzipped, and copied to a the Forcings folder (unless you change the DataPath variable in NESOSIM.py). Note that the files in Scripts/gridding shows how these forcings were generated from the raw data.


```
test_forcings/
```

 - In v1.1 this now includes gridded (50 km) test forcing data for 2018-2020: ERA-5 snowfall and winds, CDR ice concentration, OSI SAF sea ice drift vectors. 


### Ancillary Data

```
anc_data/
```
- Contains temperature-scaled intitial snow depths and the drifting soviet station data that were used for model calibration.

### Model Output

 NetCDF files used by the various analysis/plotting scripts included in this repo, which include all the model variables, can be found out:
```
output/budgets
```

NetCDF files only including the primary model variables are also generated and stored in:
```
output/final
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


Do contact me if you any any questions or thoughts on anything included here (alek.a.petty@nasa.gov). Cheers!

Alek



