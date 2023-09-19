# Using the macros

In this directory there are macros that specific for each of CAST datasets. There is one macro to generate a ROOT dataframe of each dataset, because each of them has different cuts that we need to apply. The calibration factors for each run are also included.

## Dataset 1: Ar data in which only hit channels were saved
* **applyBackgroundCuts_v2.4.0.C**  : macro to compute the cuts for each energy range, apply them to all the runs, and generate a ROOT dataframe keeping only the events that pass the X-ray cuts. Fiducial or veto cut are still not applied here. It is important to know that some of the observables have a different cut for each run. In this case, the values are input using .csv files which are also available in this repository. The calibration factors are also included using .csv files.

Usage:
```
restRoot
.L applyBackgroundCuts_v2.4.0.C
cut()
```

### **XRayTube_FOM_scalingCalculator_v2.4.0.C** : macro to compute the scaling factor to be applied to the energy ranges different that the one for 5.9 keV. The results are used by the former macro `applyBackgroundCuts_v2.4.0.C`. This means that if you just want to re-run the analysis as it was originally done, you don't need this macro as its results are incorporated into the former macro `applyBackgroundCuts_v2.4.0.C`.
### **applyFiducialAndVetoCuts.C** : macro that receives the ROOT dataframe generated with the former macro `applyBackgroundCuts_v2.4.0.C`, and allows to apply veto and fiducial cuts.

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
`g++ -o optimizeCuts_Xe_v2.4.0 optimizeCuts_Xe_v2.4.0.C -O3 `root-config --cflags --libs` -L/programas/rest/v2.4.0/lib -lRestFramework -lRestConnectors -lRestDetector -lRestGeant4 -lRestLegacy -lRestRaw -lRestTrack`.
### `optimizeCuts_Xe_v2.4.0.C` : macro that computes and applies the cuts. Like tha macros for the other datasets, it generates a ROOT dataframe containing only the X-ray-like events, but no veto or fiducial cut. It has to be run in two steps, one to generate a DF of each run, and also 2 txt files containting the duration of each run and the run number. Then these 2 .txt files are used as arguments for the `analyze` function, which gives the background rate.

Usage:
`./optimizeCuts_Xe_v2.4.0 /R10850_00001*` to run it on a single run. However, to run it over all the runs one can use the bash script `launchMacroXe_v2.4.0.sh`:
`./launchMacroXe_v2.4.0.sh`




