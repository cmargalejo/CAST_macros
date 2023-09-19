# CAST_macros
ROOT and Python scripts used in the analysis of CAST Micromegas data analysis.

## Cut definition
Under the directory `cutDefinition` one can find sample macros that are relevant to optimize and apply
the X-ray cuts to the background data.
The specific macros used for each of the 2019-2021 datasets still have to be uploaded.
* `sample_optimizeCuts_FOMcalculator_CAST.C` is a ROOT macro to optimize the X-ray cuts.
*  `SampleApplyBackgroundCuts.C`is a ROOT macro to apply the cuts, including the X-ray tube extrapolation strategy.

## Limit calculation
Under the directory `limitCalculation` one can find a python script to compute the axion-photon coupling as well
as some data and plots. It is still in a preliminary state.

## Spot simulations data
It contains the data of simulations for photons of different energies coming from a source placed at infinity. There is also a simulation of the X-ray finger data done for 3 keV.
Data provided by Jaime Ruz (I think the original author was Mike from LLNL).

**TODOs**:
* Include all the datasets, each with their own software efficiency, gas efficiency, background level and candidates.
* Write the unbinned approach.
