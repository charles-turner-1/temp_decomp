# Kappa_r Temperature Decomposition.

This contains the MATLAB code to compute an estimate of kappa_r, along with the uncertainty on the estimate.
It also contains a number of scripts which will generate sample fields as a supplement to the paper (DOI here).

.mat files contain sample data: these have been restricted to a single point in time for the gamma calculation and an array of grid cells approximating the location of the real world RAPID array for the full 240 year time series due to file size contraints. Please contact me if you wish to run the sample analyses on the full dataset.
	* KappaRandPolyfit.mat contains a global Kr field, along with ordinary least squares fits of dT/dCnat at each point. These fields are generated in the script `computeGlobalappaRValues.m`.
	* RAPID_TMP_DIC_Fields.mat contains temperature and DIC (along with decomposed components) at the location of the RAPID array. The full global fields are generated in the script `GenerateFullRedistTempTimeseries.m`. The location of the RAPID array is then selected by the function globalToRapid.
	* gammaValues.mat contains the gamma values presented in Figure A1, which are used to adjust Cnat to correct for outgassing. A sample gamma value, for the 2090-2099 decadal mean gamma value, is produced in `GenerateFullRedistTempTimeseries.m`. Generating the full timeseries will require the full, 4D fields.
	* grid_data_areas_masks.mat contains data pertaining to the ORCA1 grid: longitudes, latitudes,grid cell volumes, masks for ocean basins, etc.
	* section_location_76pt.mat contains data identifying the location of the RAPID array.

.m files contain scripts and functions which allow the user to generate estimates of kappa_r, gamma (the correction factor for outgassing), estimates of excess and redistributed temperature, and to peruse some of the calculations made in the paper. To fully reproduce the results in the paper, the full 4D fields are necessary: please get in touch if you would like these fields.

.m scripts are set up such that they should currently run on the sample data provided in the repository. 
Additionally, the code used to generate the full decomposed fields is retained in these scripts but is commented out: the procedure is essentially the same.


