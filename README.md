# argo_clustering_cokriging

This repository houses code for https://arxiv.org/abs/2211.04012, "A multivariate functional-data mixture model for spatio-temporal data: inference and cokriging"; Moritz Korte-Stapff, Drew Yarger, Stilian Stoev, and Tailen Hsing.


There are two main components of our code. The first is the code that implements our methodology in a package, based on the folder "fstmr." Each file contains the relevant code for different parts of the methodology. There are three classes: one_var (for univariate data), mult_var_ind (for multivariate data with spatial independence), and mult_var_spat (our main methodology). The second part of the code is in paper/code, which contains the code for the data analysis, simulations, images, and processing results that are described in the paper. 


The batch files in paper/ may give more insight about the organization of the code and its order. code/paper/ contains 5 folders: 1. src: utility code that does not fit in the package 2. data_preparation: code used to process the data 3. data_analysis: code that runs each aspect of the data analysis, 4. data_analysis_results: code that processes results from data_analysis, 5. simulation: the simulation code.  Most data_analysis files have a cousin file in data_analysis_results that processes the results from its cousin (get it? because the folders are their parents). For the simulation, both code that runs and plots results are in the file paper/code/simulation.

We describe the workflow in the analysis in a series of shell scripts in paper/. The data analysis consists of number 00, 01, 02, 03, and 04. 00 and 01 get and process the data, 02 runs model estimation and prediction for oxygen, 03 runs the leave out experiments, and 04 runs model estimation and prediction for nitrate. Also, 05 runs the simulation code. 

We describe what needs to be done to run 00_get_data.sh, since there are a number of moving parts:

1. The SOCCOM data should be downloaded with location specified in 00_get_data.sh as 'soccom_location'
2. The Argo snapshot should be downloaded with location specified in 00_get_data.sh as 'argo_location'
3. The grid.nc file for SOSE should be downloaded and placed as paper/data/grid.nc
4. The Roemmich and Gilson product should be downloaded with location specified in 'RG_location' of 00_get_data.sh.

Major sources of data:

1. SOCCOM Data: http://doi.org/10.6075/J0T43SZG

2. Argo Data: http://doi.org/10.17882/42182#85023

Minor sources of data:

1. RG climatology: https://sio-argo.ucsd.edu/RG_Climatology.html, ``2004-2018 1/6 x 1/6 degree mean field [t=0]'' available for download at https://sio-argo.ucsd.edu/gilson/argo_climatology/RG_ArgoClim_33pfit_2019_mean.nc.gz

2. SOSE grid: http://sose.ucsd.edu/BSOSE6_iter133_solution.html, with a direct download at http://sose.ucsd.edu/SO6/SETUP/grid.nc

3. Sea ice extent: National Snow and Ice Data Center, https://masie_web.apps.nsidc.org/pub/DATASETS/NOAA/G02135/south/monthly/geotiff/09_Sep/S_201409_concentration_v3.0.tif


