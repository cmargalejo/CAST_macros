# Quick access to the data
The final ROOT dataframes containing background data after X-ray cuts in the 0 to 12 keV energy range are kept here. Fiducial, veto or tracking cuts not applied yet, only X-ray cuts. They are purely ROOT dataframes to try to ensure their readability and use in the future.

It is important to take into account that these DF still have to be rotated by 45º if one wants to compare them with the existing simulations, as the readout in the physical setup is rotated itself. Also, the positions are defined from the point of view of the incoming photon, whereas the simulations that "see" the readout from behind, a specular image is needed.
These files are:
* background_df_all_runs_incArDataset2HitChannels_v2.4.0_final_withRateAna_withAllSC.root
* background_df_all_runs_dataset2_subsets2and3_v2.4.0_extrapolatedCuts_withRateAna.root
* MISSING

# Using the macros

In this directory there are macros that specific for each of the CAST datasets. There is one macro to generate a ROOT dataframe of each dataset, because each of them has different cuts that we need to apply.

The macros do compute the cuts for each of the energy ranges and apply them. However, for some observables (sigmas) the cut must be different for each run, so it is not computed in the macro but it is given as external data. In that cases, the cuts were computed based on the efficiency and the value for eahc run is listed in .csv files that are read my the relevant macro. The calibration factors for each run are also included using .csv files. All these .csv files are stored in the `data` directory.



## Dataset 1: Ar data in which only hit channels were saved
### applyBackgroundCuts_v2.4.0.C
Macro to compute the cuts for each energy range, apply them to all the runs, and generate a ROOT dataframe keeping only the events that pass the X-ray cuts. Fiducial or veto cut are still not applied here. It is important to know that some of the observables have a different cut for each run. In this case, the values are input using .csv files which are also available in this repository. The calibration factors are also included using .csv files.

Usage:
```
restRoot
.L applyBackgroundCuts_v2.4.0.C
cut()
```

### XRayTube_FOM_scalingCalculator_v2.4.0.C
Macro to compute the scaling factor to be applied to the energy ranges different that the one for 5.9 keV. The results are used by the former macro `applyBackgroundCuts_v2.4.0.C`. This means that if you just want to re-run the analysis as it was originally done, you don't need this macro as its results are incorporated into the former macro `applyBackgroundCuts_v2.4.0.C`.
### applyFiducialAndVetoCuts.C
Macro that receives the ROOT dataframe generated with the former macro `applyBackgroundCuts_v2.4.0.C`, and allows to apply veto and fiducial cuts.

Usage:
```
restRoot
.L applyFiducialAndVetoCuts.C
cut("cutstring")
```
where the cutstring can be, e.g., `energy_keV>0`.

## Dataset 2: Ar data in which all channels were saved
The cuts are different to those in dataset 1, so different macros are used. However, the procedure is exactly the same. The macros are:
* `XRayTube_FOM_scalingCalculator_v2.4.0_dataset2.C`
* `applyBackgroundCuts_v2.4.0_dataset2.C`
* `applyFiducialAndVetoCuts_dataset2.C`

## Dataset 3 : Xe data
This analysis ended up being a bit more complex. There are also a lot of data, so the macro has to be compiled and executed as a program to avoid problems of it stopping or running out of memory. Use the following command to compile it:
``` 
g++ -o optimizeCuts_Xe_v2.4.0 optimizeCuts_Xe_v2.4.0.C -O3 `root-config --cflags --libs` -L/programas/rest/v2.4.0/lib -lRestFramework -lRestConnectors -lRestDetector -lRestGeant4 -lRestLegacy -lRestRaw -lRestTrack
```
### optimizeCuts_Xe_v2.4.0.C
Macro that computes and applies the cuts. Like tha macros for the other datasets, it generates a ROOT dataframe containing only the X-ray-like events, but no veto or fiducial cut. It has to be run in two steps, one to generate a DF of each run, and also 2 txt files containting the duration of each run and the run number. Then these 2 .txt files are used as arguments for the `analyze` function, which gives the background rate.

Usage:
`./optimizeCuts_Xe_v2.4.0 /R10850_00001*` to run it on a single run. However, to run it over all the runs one can use the bash script `launchMacroXe_v2.4.0.sh`:
`./launchMacroXe_v2.4.0.sh`

### XRayTube_FOM_scalingCalculator_Xe_v2.4.0.C
Same as the cases with Ar. You don't need to run it again as the results are already incorporated in `optimizeCuts_Xe_v2.4.0.C`, but it is still kept here for reference. 




