# Using the macros

In this directory there are macros that specific for each of CAST datasets. There is one macro to generate a ROOT dataframe of each dataset, because each of them has different cuts that we need to apply. The calibration factors for each run are also included.

 *Dataset 1: Ar data in which only hit channels were saved
    *applyBackgroundCuts_v2.4.0.C'  : macro to compute the cuts for each energy range, apply them to all the runs, and generate a ROOT dataframe keeping only the events that pass the X-ray cuts. Fiducial or veto cut are still not applied here.
    Usage:
    ```
   restRoot
    .L applyBackgroundCuts_v2.4.0.C
    cut()
    ```
    It is important to know that some of the observables have a different cut for each run. In this case, the values are input using .csv files which are also available in this repository. The calibration factors are also included using .csv files.
    *XRayTube_FOM_scalingCalculator_v2.4.0.C : macro to compute the scaling factor to be applied to the energy ranges different that the one for 5.9 keV. The results are incorporated into the former macro `applyBackgroundCuts_v2.4.0.C`
    *applyFiducialAndVetoCuts.C : macro that receives the ROOT dataframe generated with the former macro 'applyBackgroundCuts_v2.4.0.C', and allows to apply veto and fiducial cuts.
    Usage:
    ```
    restRoot
    .L applyFiducialAndVetoCuts.C
    cut("cutstring")
    ```
    where the cutstring can be, e.g., 'energy_keV>0'.

*Dataset 2: Ar data in which all channels were saved
*Dataset 3 : Xe data


