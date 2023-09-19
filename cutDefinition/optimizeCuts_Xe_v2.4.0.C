// Author: Cristina Margalejo and Sebastian Schmidt, November 2022
// Updated for the CAST Xe 2020-2021 data taking campaign, February 2023
// Based on former macro by Konrad Altenmüller, June 2022
//
//
//
#include <TCanvas.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TMath.h>
#include <string.h>
#include "Math/Minimizer.h"
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <string>
#include "ROOT/RDataFrame.hxx"
#include <fstream>
// #include "ROOT/RVec.hxx"


#include "prettyprint.hpp"

// trying to solve the issue of REST libraries not loading in September 2023.
#include "TSystem.h"


using namespace ROOT::VecOps;

// custom type alias for a bit easier reading of the map mapping
// each run number to the correct cuts for X / Y gauss sigma.
typedef std::map<std::pair<int, int>, std::tuple<std::pair<double, double>, std::pair<double, double>, std::pair<double, double>>> RunCutMap;
//               ^--- run number
//                                        ^---  & ^-- X cuts         ^---  ^ ^--- Y cuts
typedef std::map<int, int> EnergyCalibMap;
typedef std::vector<std::vector<std::pair<double, double>>> EnergyCutsVec;
//           ^--- energy index
//                       ^--- observable in specific order
//                                        ^---  & ^--- cut low & high


template <typename M>
typename M::mapped_type verboseAt(const M& map, const typename M::key_type& key){
    if(map.count(key) > 0){
	return map.at(key);
    }
    else{
	std::stringstream s;
	s << "The key " << key << " does not exist in the given map!" << std::endl;
	//const auto s = "The key " + std::to_string(key) + " does not exist in the given map!";
	throw std::out_of_range(s.str());
    }
}

template <typename V>
typename V::value_type verboseAt(const V& vec, std::size_t idx){
    if(idx < vec.size()){
	return vec.at(idx);
    }
    else{
	std::stringstream s;
	s << "The index " << idx << " does not exist in the given vector!" << std::endl;
	throw std::out_of_range(s.str());
    }
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//TODO: write documentation!
// Macro to produce a dataframe keeping only the events that pass the cuts.
//
// Functions: defined the at the beginning of each function.
//
//Usage:
//restRootvccv
//.L macroFilename.C
//cut()
///////////////////////////////////////////////////////////////////////////////
//********************************
// Get run duration function
//********************************
int GetRunDuration(ROOT::RDataFrame data_filtered, double binw = 60, int rate_threshold = 10) { //1 minute bins, --> threshold = 10 couns per minute
    int bincount;
    int run_duration = 0;
    // calculate nbins to get desired bin width
    // vector<double> timestamps = data_filtered.Take<double>("timeStamp").GetValue();
    auto result_ts = data_filtered.Take<double>("timeStamp");
    std::vector<double> timestamps = result_ts.GetValue();

    double tmax = data_filtered.Max("timeStamp").GetValue();
    double tmin = data_filtered.Min("timeStamp").GetValue();


    tmax = tmin + (round((tmax - tmin + binw / 2) / binw) * binw);
    // create histogram
    TH1D* ht = new TH1D("timestamps", "timestamps", (tmax - tmin) / binw, tmin, tmax);
    for (auto i : timestamps) {
        ht->Fill(i);
    }
    binw = ht->GetBinWidth(1);
    // TCanvas* c00 = new TCanvas();
    // ht->Draw();
    // check if bins are above threshold
    int nbins = ht->GetNbinsX();
    for (unsigned int i = 0; i < nbins; i++) {
        bincount = ht->GetBinContent(i);
        if (bincount > rate_threshold) {
            run_duration += binw;
        }
    }
    std::cout << "run_duration = " << run_duration << std::endl;
    delete(ht);
    return run_duration;
    // cout << "Bin width: " << binw << endl;
}


enum DataKind { dkBackground = 0, dkXrayTube = 1 };

//****************************************************************
// Function used to calibrate the 6 energy ranges. Calibration factors
//based on X-ray tube data.
//****************************************************************
// int getCalibrationFactor(int energyIdx) {
//     switch(energyIdx) {
//         case 1: return 1356;
//         case 2: return 1472;
//         case 3: return 1512;
//         case 4: return 1478;
//         case 5: return 1458;
//         case 6: return 1446;
//         default:
//             std::cout << "Invalid energy range given: `" << energyIdx
//                       << "`. Only values in [1, 6] allowed" << std::endl;
//             assert(false);
//             return 0;
//     }
// }

//***************************************************************
// Function to get the energy index given an energy. This is used in the
// `applyCut` function to determine which cut parameters to look at.
//***************************************************************

int getEnergyIdx(double energy){
    const std::vector<std::pair<double, double>> energyRanges = {
        { 0.00, 1.75 },
        { 1.75, 2.50 },
        { 2.50, 3.75 },
        { 3.75, 5.25 },
        { 5.25, 7.00 },
        { 7.00, 12.0 } };
    int idx = 1;
    for(auto range : energyRanges){
        if(energy > range.first && energy <= range.second){
            return idx;
        }
        idx++;
    }
    return -1;
}


//***************************************************************
// Function to get the calibrated enery, i.e., to turn ADC into keV.
//***************************************************************
// ROOT::RDF::RNode calcEnergy(ROOT::RDF::RNode df, int energyIdx, DataKind dKind){
//     switch(dKind) {
//         case dkBackground:
//             return df.Define("energy_keV","tckAna_MaxTrackEnergy/1070");
//         case dkXrayTube:
//             auto calcStr = "tckAna_MaxTrackEnergy/" + to_string(getCalibrationFactor(energyIdx));
//             return df.Define("energy_keV", calcStr);
//     }
// }

//***************************************************************
// Function to access the relevant scaling factor or paramter p_i.
//***************************************************************
std::pair<double, double> getScalingFactor(int energyIdx, std::string obs){
    // Given the observable `obs`, return the correct initial scaling factor
    // for the given energy
    std::map<std::string, std::map<int, std::pair<double, double>>> map = {
/*
        { "tckAna_MaxTrackEnergyBalanceXY", {
                { 1,  { 2.27, 1.86 }  },
                { 2,  { 2.45, 2.09 }  },
                { 3,  { 1.48, 1.17 }  },
                { 4,  { 1.30, 1.06 }  },
                { 5,  { 1.0, 1.0 }  },
                { 6,  { 0.94, 0.94 }  },
          }
        },
        { "tckAna_MaxTrack_XZ_GaussSigmaX", {
                { 1,  { 1, 0.83 }  },
                { 2,  { 1, 0.8 }  },
                { 3,  { 0.76, 1.8 }  },
                { 4,  { 2.06, 0.7 }  },
                { 5,  { 1.0, 1.0 }  },
                { 6,  { 2.06, 1.0 }  },
          }
        },
        { "tckAna_MaxTrack_YZ_GaussSigmaY", {
                { 1,  { 1.08, 1.07 }  },
                { 2,  { 0.76, 1.04 }  },
                { 3,  { 0.68, 1.04 }  },
                { 4,  { 1.56, 0.9 }  },
                { 5,  { 1.0, 1.0 }  },
                { 6,  { 1.72, 1.33 }  },
          }
        },
        { "hitsAna_z2Sigma", {
                { 1,  { 0.0, 3.71 }  },//Ar: 3.4 before revising the cuts, 10.5 after --> *factor 3.088 | Xe gauss: 0.23 * 3.088 = 0.74, Xe RMS: 3.71 * 3.088 = 11.5
                { 2,  { 0.0, 5.81 }  }, //Ar: 2.2 before revising the cuts, 7.0 after --> *factor 3.18 | Xe gauss: 0.78 * 3.18 = 2.48, Xe RMS: 5.81 * 3.18 = 18.5
                { 3,  { 0.0, 1.58 }  }, //Ar: 1.8 before revising the cuts, 5.8 after --> *factor 3.22 | Xe gauss: 0.87 * 3.22 = 2.80, Xe RMS: 1.58 * 3.22 = 5.1
                { 4,  { 0.0, 1.06 }  },
                { 5,  { 0.0, 1.0 }  },
                { 6,  { 0.0, 0.74 }  },
          }
        },
        { "tckAna_MaxTrackxySigmaGausBalance", {
                { 1,  { 1.99, 3.84 }  },
                { 2,  { 1.99, 3.84 }  },
                { 3,  { 1.87, 3.84 }  },
                { 4,  { 1.99, 1.70 }  },
                { 5,  { 1.0, 1.0 }  },
                { 6,  { 0.8, 1.05 }  },
          }
        },
        { "hitsAna_xySkew", {
                { 1,  { 0.76, 1.89 }  },
                { 2,  { 1.29, 2.19 }  },
                { 3,  { 1.37, 1.27 }  },
                { 4,  { 1.70, 0.84 }  },
                { 5,  { 1.0, 1.0 }  },
                { 6,  { 0.71, 1.65 }  },
          }
        },
        { "hitsAna_zSkew", {
                { 1,  { 1.0, 1.0 }  },
                { 2,  { 1.0, 1.0  }  },
                { 3,  { 1.0, 1.0 }  },
                { 4,  { 1.0, 1.0  }  },
                { 5,  { 1.0, 1.0 }  },
                { 6,  { 1.0, 1.0 }  },
          }
        },
*/
      { "tckAna_MaxTrackEnergyBalanceXY", {
                { 1,  { 2.27, 2.84 }  },
                { 2,  { 2.45, 3.16 }  },
                { 3,  { 1.3, 1.48 }  },
                { 4,  { 1.12, 1.16 }  },
                { 5,  { 1.0, 1.0 }  },
                { 6,  { 0.82, 0.84 }  },
          }
        },
        { "tckAna_MaxTrack_XZ_GaussSigmaX", {
                { 1,  { 0.49, 1.03 }  },
                { 2,  { 0.33, 0.93 }  },
                { 3,  { 1.0, 0.68 }  },
                { 4,  { 0.44, 0.98 }  },
                { 5,  { 1.0, 1.0 }  },
                { 6,  { 0.59, 1.02 }  },
          }
        },
        { "tckAna_MaxTrack_YZ_GaussSigmaY", {
                { 1,  { 0.46, 1.08 }  },
                { 2,  { 0.32, 1.08 }  },
                { 3,  { 0.95, 0.85 }  },
                { 4,  { 0.41, 0.77 }  },
                { 5,  { 1.0, 1.0 }  },
                { 6,  { 0.56, 0.94 }  },
          }
        },
        { "hitsAna_z2Sigma", {
                { 1,  { 0.0, 11.6 }  }, // Xe RMS: 3.76 * 3.088 = 11.61
                { 2,  { 0.0, 19.7 }  }, // Xe RMS: 6.21 * 3.18 = 19.7
                { 3,  { 0.0, 5 }  }, // Xe RMS: 1.55 * 3.22 = 5
                { 4,  { 0.0, 1.0 }  },
                { 5,  { 0.0, 1.0 }  },
                { 6,  { 0.0, 0.72 }  },
          }
        },
        { "tckAna_MaxTrackxySigmaGausBalance", {
                { 1,  { 2.79, 1.42 }  },
                { 2,  { 3.17, 1.46 }  },
                { 3,  { 0.7, 1.46 }  },
                { 4,  { 2.57, 1.46 }  },
                { 5,  { 1.0, 1.0 }  },
                { 6,  { 1.68, 1.46 }  },
          }
        },
        { "hitsAna_xySkew", {
                { 1,  { 2.05, 1.73 }  },
                { 2,  { 2.14, 2.22 }  },
                { 3,  { 1.19, 1.05 }  },
                { 4,  { 1.1, 1.24 }  },
                { 5,  { 1.0, 1.0 }  },
                { 6,  { 0.9, 2.12 }  },
          }
        },
        { "hitsAna_zSkew", {
                { 1,  { 1.0, 1.0 }  },
                { 2,  { 1.0, 1.0 }  },
                { 3,  { 1.0, 1.0 }  },
                { 4,  { 1.0, 1.0 }  },
                { 5,  { 1.0, 1.0 }  },
                { 6,  { 1.0, 1.0 }  },
          }
        },

    };
    return verboseAt(verboseAt(map, obs), energyIdx);
}


//***************************************************************
// Function to access the inital values of the cuts, both for the 55Fe
// calibrations at CAST and for the 5.9 keV calibration in the X-ray tube.
//***************************************************************
// TODO: think about whether one might want to change individual cut values
// e.g. z2Sigma separately
std::pair<double, double> getInitialCut(DataKind dKind, std::string obs){
    // Given the data kind and observable, return the initial cut.

    std::map<std::string, std::pair<double, double>> castMap = {
        { "tckAna_MaxTrackEnergyBalanceXY",  { -0.16, 0.2 }  },
        { "tckAna_MaxTrack_XZ_GaussSigmaX",  { 0.35, 1.2 }  },
        { "tckAna_MaxTrack_YZ_GaussSigmaY",  { 0.35, 1.2 }  },
        { "hitsAna_z2Sigma",  { 0.0, 0.23 }  },
        { "tckAna_MaxTrackxySigmaGausBalance",  { -0.3, 0.3 }  },
        { "hitsAna_xySkew",  { -0.8, 0.8 }  },
        { "hitsAna_zSkew",  { -10, 10 }  },
    };
    /*
    std::map<std::string, std::pair<double, double>> xrayTubeMap = {
        { "tckAna_MaxTrackEnergyBalanceXY",  { -0.20, 0.21 }  },
        { "tckAna_MaxTrack_XZ_GaussSigmaX",  { 0.17, 1.75 }  },
        { "tckAna_MaxTrack_YZ_GaussSigmaY",  { 0.25, 1.35 }  },
        { "hitsAna_z2Sigma",  { 0.0, 1.55 }  },
        { "tckAna_MaxTrackxySigmaGausBalance",  { -0.5, 0.26 }  },
        { "hitsAna_xySkew",  { -0.98, 0.74 }  },
        { "hitsAna_zSkew",  { -10, 10 }  },
    };
    */
     std::map<std::string, std::pair<double, double>> xrayTubeMap = {
        { "tckAna_MaxTrackEnergyBalanceXY",  { -0.20, 0.15 }  },
        { "tckAna_MaxTrack_XZ_GaussSigmaX",  { 0.39, 1.79 }  },
        { "tckAna_MaxTrack_YZ_GaussSigmaY",  { 0.41, 1.85 }  },
        { "hitsAna_z2Sigma",  { 0.0, 1.45 }  },
        { "tckAna_MaxTrackxySigmaGausBalance",  { -0.28, 0.68 }  },
        { "hitsAna_xySkew",  { -0.42, 0.82 }  },
        { "hitsAna_zSkew",  { -10, 10 }  },
    };

/*
    std::map<std::string, std::pair<double, double>> castMap = {
        { "tckAna_MaxTrackEnergyBalanceXY",  { -100, 100 }  },
        { "tckAna_MaxTrack_XZ_GaussSigmaX",  { -100, 100 }  },
        { "tckAna_MaxTrack_YZ_GaussSigmaY",  { -100, 100 }  },
        { "hitsAna_z2Sigma",  { -100, 100 }  },
        { "tckAna_MaxTrackxySigmaGausBalance",  { -100, 100 }  },
        { "hitsAna_xySkew",  { -100, 100 }  },
        { "hitsAna_zSkew",  { -100, 100 }  },
    };
    std::map<std::string, std::pair<double, double>> xrayTubeMap = {
        { "tckAna_MaxTrackEnergyBalanceXY",  { -100, 100 }  },
        { "tckAna_MaxTrack_XZ_GaussSigmaX",  { -100, 100 }  },
        { "tckAna_MaxTrack_YZ_GaussSigmaY",  { -100, 100 }  },
        { "hitsAna_z2Sigma",  { -100, 100 }  },
        { "tckAna_MaxTrackxySigmaGausBalance",  { -100, 100 }  },
        { "hitsAna_xySkew",  { -100, 100 }  },
        { "hitsAna_zSkew",  { -100, 100 }  },
    };
*/

    switch(dKind) {
        case dkBackground:
            return verboseAt(castMap, obs);
        case dkXrayTube:
            return verboseAt(xrayTubeMap, obs);
        default:
            assert(false);
            return std::make_pair(0.0, 0.0);
    };
}
//***************************************************************
// Function to get the cut scaling applied to each energy range.
//***************************************************************
double getCutScaling(int energyIdx){
    //std::vector<double> vec = { 0.2, 0.18, 0.15, 0.15, 0.0, 0.05 };
    //std::vector<double> vec = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    //std::vector<double> vec = { 0.05, 0.25, 0.2, 0.1, 0.0, 0.05 };
    std::vector<double> vec = { 0.15, 0.2, 0, 0.1, 0.0, 0.0 };
    return vec[energyIdx-1];
}

EnergyCutsVec cutVec; // global variable.
//***************************************************************
// Function to decide if an event pass the cuts or not.
//***************************************************************
bool applyCut(double energyBalance, double gaussSigmaX, double gaussSigmaY, double z2Sigma,
              double sigmaGaussBalance, double xySkew, double zSkew, double energy){
    const int energyIdx = getEnergyIdx(energy);
    if(energyIdx < 0){
        return false;
    }
    const auto vec = cutVec[energyIdx-1];
    // for(auto val : vec){
    //     std::cout << val.first << " : " << val.second << "\t" << std::flush;
    // }
    // std::cout << std::endl;
    if(energyBalance       < vec[0].first || energyBalance         > vec[0].second){
        //std::cout << "energyBalance = " << energyBalance << std::endl;
        return false;
    }
     else if(gaussSigmaX         < vec[1].first || gaussSigmaX           > vec[1].second){
        //std::cout << "gaussSigmaX = " << gaussSigmaX << " vec.first =  " << vec[1].first <<  " vec.second =  " << vec[1].second << std::endl;
        return false;
    }
    else if(gaussSigmaY         < vec[2].first || gaussSigmaY           > vec[2].second){
        //std::cout << "gaussSigmaY = " << gaussSigmaY << std::endl;
        return false;
    }
    else if(z2Sigma             < vec[3].first || z2Sigma               > vec[3].second){
        //std::cout << "z2Sigma = " << z2Sigma << std::endl;
        return false;
    }
    else if(sigmaGaussBalance   < vec[4].first || sigmaGaussBalance     > vec[4].second){
        //std::cout << "sigmaGaussBalance = " << sigmaGaussBalance << std::endl;
        return false;
    }
    else if(xySkew              < vec[5].first || xySkew                > vec[5].second){
        //std::cout << "xySkew = " << xySkew << std::endl;
        return false;
    }
    else if(zSkew               < vec[6].first || zSkew                 > vec[6].second){
        //std::cout << "zSkew = " << zSkew << std::endl;
        return false;
    }
    else {
	//if(z2Sigma > 2.9){
	//    std::cout << "event passsing with " << z2Sigma << " with cut range " << vec[3].first << " to " << vec[3].second << std::endl;
	//}
        return true;
    }
    assert(false); // this can never happen!
    return false;  // this clearly ever happens, because assert will kill program
}

//*****************************************************************
// Function to generate the cut string.
// This cut string will be the input to the main cut() function.
//*****************************************************************
//std::string genCutString(int energyIdx, DataKind dKind, double cutScaling){
//    const auto names = {
//        "tckAna_MaxTrackEnergyBalanceXY",
//        "tckAna_MaxTrack_XZ_GaussSigmaX",
//        "tckAna_MaxTrack_YZ_GaussSigmaY",
//        "hitsAna_z2Sigma",
//        "tckAna_MaxTrackxySigmaGausBalance",
//        "hitsAna_xySkew",
//        "hitsAna_zSkew"};
//
//    std::string result = "";
//    for(auto name : names) {
//        auto initCut = getInitialCut(dKind, name);
//        auto scaling = getScalingFactor(energyIdx, name);
//        // define scaled cuts
//        auto sLow = initCut.first * scaling.first;
//        auto sHigh = initCut.second * scaling.second;
//        auto scaledDifference = abs(sHigh - sLow) * cutScaling / 2.0;
//        cout << "_________________________________________________"  << endl;
//        cout << "low cut " << sLow << "+" << scaledDifference << ", high cut " << sHigh << "-" << scaledDifference <<  " scaling = " << scaledDifference << endl;
//        //cout << "scaled difference: " << scaledDifference << endl;
//        if (sLow != 0) sLow = sLow + scaledDifference;
//        sHigh = sHigh - scaledDifference;
//        cout << "new low cut = " << sLow << ", new high cut = " << sHigh << endl;
//        cout << "_________________________________________________"  << endl;
//        if(result.size() > 0){
//            result = result + "&&";
//        }
//        result = result + name + ">" + to_string(sLow) + "&&" +
//            name + "<" + to_string(sHigh);
//    }
//    return result;
//}

void fillGlobalCutVec(int run, RunCutMap csvData){
    // 1. look up the correct gauss sigma from CSV
    // 2. loop over all energy ranges
    // 3. get all other cut values
    // 4. assign to global vector
    const auto names = {
        "tckAna_MaxTrackEnergyBalanceXY",
        "tckAna_MaxTrack_XZ_GaussSigmaX",
        "tckAna_MaxTrack_YZ_GaussSigmaY",
        "hitsAna_z2Sigma",
        "tckAna_MaxTrackxySigmaGausBalance",
        "hitsAna_xySkew",
        "hitsAna_zSkew"
    };
    cutVec.clear();
    for(size_t energyIdx = 1; energyIdx <= 6; energyIdx++){
        auto cutScaling = getCutScaling(energyIdx);
        std::pair<double, double> cut;
        std::vector<std::pair<double, double>> energyVec;
        int idx = 0;
        for(auto name : names) {
            auto initCut = getInitialCut(dkBackground, name);
            auto scaling = getScalingFactor(energyIdx, name);
            // define scaled cuts
            auto sLow = initCut.first * scaling.first;
            auto sHigh = initCut.second * scaling.second;
            auto scaledDifference = abs(sHigh - sLow) * cutScaling / 2.0;
            if (sLow != 0) sLow = sLow + scaledDifference;
            sHigh = sHigh - scaledDifference;
            cut = std::make_pair( sLow, sHigh );
            if(idx != 1 && idx != 2 && idx != 3){
                energyVec.push_back(cut);
            }
            else if(idx == 1){ //it can be 1, 2 or 3, because we push back the csv data only once.
                //std::cout << "run: " << run << ", energy index: " << energyIdx << std::endl;
                //std::cout << csvData << std::endl;
                auto gaussSigmaCuts = verboseAt(csvData, std::make_pair(run, energyIdx));
                energyVec.push_back( std::get<0>(gaussSigmaCuts));  // corresponds to X cuts
                //std::cout << "gauss sigma cut: " << std::get<0>(gaussSigmaCuts).first << std::endl;
                //std::cout << "gauss sigma cut: " << std::get<0>(gaussSigmaCuts).second << std::endl;
                energyVec.push_back( std::get<1>(gaussSigmaCuts) ); // corresponds to Y cuts
                //std::cout << "gauss sigma cut: " << std::get<1>(gaussSigmaCuts).first << std::endl;
                //std::cout << "gauss sigma cut: " << std::get<1>(gaussSigmaCuts).second << std::endl;
                energyVec.push_back( std::get<2>(gaussSigmaCuts) ); // corresponds to Z cuts
                //std::cout << "gauss sigma cut: " << std::get<2>(gaussSigmaCuts).first << std::endl;
                //std::cout << "gauss sigma cut: " << std::get<2>(gaussSigmaCuts).second << std::endl;
            }
            idx++;
        }
	cutVec.push_back(energyVec);
    }
}

RunCutMap readGaussSigmaCuts(std::string fname, int energyIdx){
    std::string line;
    std::ifstream file(fname);
    RunCutMap result;
    int lineCount = 0;
    if(file.is_open()){
        while(std::getline(file, line, '\n')){ // for each line
            // probly have to replace `line` by a stringstream?
            if(lineCount == 0){
                lineCount++;
                continue; // skip the header line
            }
            std::string cell;
            int col = 0;
            int run = 0;
            std::pair<double, double> x;
            std::pair<double, double> y;
            std::pair<double, double> z;
            std::stringstream sline(line);
                if(line.size() > 0){
            while(std::getline(sline, cell, ',')){
                switch(col){
                    case 0:
                        run = std::stoi(cell);
                        break;
                    case 1:
                        x.first = std::stof(cell);
                        break;
                    case 2:
                        x.second = std::stof(cell);
                        break;
                    case 3:
                        y.first = std::stof(cell);
                        break;
                    case 4:
                        y.second = std::stof(cell);
                        break;
                    case 5:
                        z.first = std::stof(cell);
                        break;
                    case 6:
                        z.second = std::stof(cell);
                        break;
                }
                col++;
            }
            //std::cout <<"energy = " << energyIdx << ", run = " << run << ", x = " << x.first << ", y = " << y.first << ", z = " << z.first << std::endl;
            //result[std::make_pair(run, energyIdx)] = std::make_pair( x, y );
            result[std::make_pair(run, energyIdx)] = std::make_tuple( x, y, z );
            }
        }
    }
    file.close();
    return result;
}

EnergyCalibMap readEnergyCalibMap(std::string fname){
    // TODO: would be nicer to have a single CSV parsing interface, but well
    std::string line;
    std::ifstream file(fname);
    EnergyCalibMap result;
    int lineCount = 0;
    if(file.is_open()){
        while(std::getline(file, line, '\n')){ // for each line
            // probly have to replace `line` by a stringstream?
            if(lineCount == 0){
                lineCount++;
                continue; // skip the header line
            }
            std::string cell;
            int col = 0;
            int run = 0;
            int factor = 0;
                std::stringstream sline(line);
                if(line.size() > 0){
            while(std::getline(sline, cell, ',')){
                switch(col){
                    case 0:
                        run = std::stoi(cell);
                        break;
                    case 1:
                        factor = std::stoi(cell);
                        break;
                }
                col++;
            }
            result[run] = factor;
            }
        }
    }
    file.close();
    return result;
}


// *****************************************************************
// Function to define custom columns and filters for the initial DF.
// *****************************************************************
//ROOT::RDF::RNode calcEnergyForBackground(ROOT::RDataFrame df, int runNumber, EnergyCalibMap energyCalibMap){
ROOT::RDF::RNode calcEnergyForBackground(ROOT::RDF::RNode df, int runNumber, EnergyCalibMap energyCalibMap){
    // 1. look up the correct calibration factor
    const double calibFactor = (double) verboseAt(energyCalibMap, runNumber);
    //auto calibEnergy = [calibFactor](const RVec<double> &vec){
    //    return vec / calibFactor;
    //};
    // auto calibEnergy = [calibFactor](const RVec<double> &vec){
    // 	std::vector<double> result(vec.size());
    //     for(size_t i = 0; i < vec.size(); i++){
    // 	    result[i] = vec[i] / calibFactor;
    // 	}
    // 	return result;
    // };
    auto calibEnergy = [calibFactor](const double &val){
	    return val / calibFactor;
    };
    //std::cout << typeid(calibEnergy).name() << std::endl;
    return df.Define("energy_keV", calibEnergy, {"tckAna_MaxTrackEnergy"});
}

ROOT::RDF::RNode calcEnergyAndFilterTracks(ROOT::RDataFrame df, int runNumber, EnergyCalibMap energyCalibMap){
    ROOT::RDF::RNode result = df; //ROOT::RDataFrame(0);
    // calibrated energy

    result = calcEnergyForBackground(result, runNumber, energyCalibMap);
    // ******************************

    // to initialize filtered dataframe
    //result = result.Filter("eventID > 0");
    result = result.Filter("eventID > 0 && tckAna_nTracks_X == 1 && tckAna_nTracks_Y == 1");

    // tracks cut
    //result = result.Filter("tckAna_nTracks_X == 1");
    //result = result.Filter("tckAna_nTracks_Y == 1");

    // energy cut
    result = result.Filter("energy_keV <= 12");
    //result = result.Filter("hitsAna_xMean*hitsAna_xMean+hitsAna_yMean*hitsAna_yMean < 100 ");

    result.Display()->Print();
    /*
    // Print columns' names
    auto colNames = result.GetColumnNames();
    for (auto &&colName : colNames) std::cout << colName << std::endl;
    */
    auto nEntriesBeforeCuts = result.Count().GetValue();
    std::cout << "Number of 1 track events before X-ray cuts: " << nEntriesBeforeCuts << std::endl;
    return result;
}


std::pair<std::string,int> applyBackgroundCuts(std::string fname, std::string commonPrefix,
                                RunCutMap runCutMap, EnergyCalibMap energyCalibMap){
    // 1. generate DF from the input file
    ROOT::RDataFrame inputDf("AnalysisTree", fname);
    // 1.1 compute the duration of each run
    double binw = 60;
    int rate_threshold = 10;  // counts per minute
    int run_duration = GetRunDuration(inputDf, binw, rate_threshold);
    std::cout << "\n\033[1;35mrun duration: " << run_duration << " s = " << (double)run_duration / 3600
         << " h = " << (double)run_duration / (3600 * 24) << " d" << std::endl;
    // 2. find run number in root file using REST tools
    std::string prefix = commonPrefix + "/R";
    const auto runNumberDigits = 5;
    const int substrStart = fname.find(prefix) + prefix.size();
    const auto sub = fname.substr(substrStart, runNumberDigits);
    std::cout << "The sub is : " << sub << std::endl;
    int runNumber = std::stoi( sub );

    // 3. compute energy of the data & filter to single tracks
    auto bgData = calcEnergyAndFilterTracks(inputDf, runNumber, energyCalibMap);

    // 4. (re-)generate the global cuts vector for this run number
    fillGlobalCutVec(runNumber, runCutMap);
    // 5. apply all cuts required for this run to get a "background" df
    bgData = bgData.Filter( applyCut, { "tckAna_MaxTrackEnergyBalanceXY",
                               "tckAna_MaxTrack_XZ_GaussSigmaX",
                               "tckAna_MaxTrack_YZ_GaussSigmaY",
                               "hitsAna_z2Sigma",
                               "tckAna_MaxTrackxySigmaGausBalance",
                               "hitsAna_xySkew",
                               "hitsAna_zSkew",
                               "energy_keV"} );
    // 6. if events left in background DF:
    //   6a). generate an output filename based on input & run number?
    //   6b). write `Snapshot` of dataframe to file and return filename
    // 7. else: return empty string
    std::string result = "";
    auto cnt = bgData.Count();
    if(bgData.Count().GetValue() > 0){
        result = "background_df_run_" + std::to_string(runNumber) + "_v2.4.0_RMSsigmaZ_symmetricSigmaXY_final_withSC.root";
        bgData.Snapshot("AnalysisTree", result);
    }
    //bgData.Display({"hitsAna_z2Sigma"}, 100)->Print();
    bgData.Display()->Print();
        std::cout << "Number of events left after cuts in run " << runNumber << " = " << *cnt << std::endl;
    return std::make_pair(result,run_duration);
        // Print columns' names
    //auto colNames = bgData.GetColumnNames();
    //for (auto &&colName : colNames) std::cout << colName << std::endl;
}


// Main function
//***************************************************************************************
// Cut function: it is the funciton that is intended to be used by the user.
// The arguments are:
// cutString: it is now deprecated and should be empty, as it will be
// automatically generated by the function genCutString().
// energyIdx: it goes from 1 to 6, and defined the energy range in which we want
// to make the computation.
// cutScaling: in case we want to modify the parameters in getScalingFactor().
// Use example: cut ("",1,0.1). Here, the cutString is empty because it will be auto-generated,
// the energy range is 1, i.e., the one corresponding to the 1.5 keV peak, and 0.1 or 10%
// is how much the cut will be tightened. In this case, it will narrow the original cut window
// by 5% on each side. If the lower cut is 0, the correction is only applied to the
// upper cut.
//***************************************************************************************

void cut(std::string runPrefix) {

    // ******************************
    // Important variables and declarations
    // ******************************

    std::string processName_tckAna = "tckAna_";                       // process prefix
    std::string processName_hitsAna = "hitsAna_";                         // process prefix
    std::string path = "/storage/cast/SR2019/analysis/xenon/official/v2.4.0/";  // data path

    // Introduce the filenames
    // std::vector<std::string> runCutsCsv = {"/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange1.csv", "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange2.csv","/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange3.csv","/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange4.csv","/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange5.csv","/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange6.csv"};
    //std::vector<std::string> runCutsCsv = {"/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange1_v2_scaled.csv", "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange2_v2_scaled.csv","/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange3_v2_scaled.csv","/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange4_v2_scaled.csv","/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange5_v2_scaled.csv","/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange6_v2_scaled.csv"};
    //std::vector<std::string> runCutsCsv = {"/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange1_v2.csv", "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange2_v2.csv","/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange3_v2.csv","/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange4_v2.csv","/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange5_v2.csv","/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange6_v2.csv"};

    //std::vector<std::string> runCutsCsv = { "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange1_Xe_constantEff.csv", "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange2_Xe_constantEff.csv", "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange3_Xe_constantEff.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange4_Xe_constantEff.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange5_Xe.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange6_Xe_constantEff.csv" };
    //std::vector<std::string> runCutsCsv = { "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange1_Xe_maxFOM.csv", "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange2_Xe_maxFOM.csv", "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange3_Xe_maxFOM.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange4_Xe_maxFOM.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange5_Xe.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange6_Xe_maxFOM.csv" };
    //v 2.4.0 with low sigma Z=0 and rescaled sigmas X and Y
    //std::vector<std::string> runCutsCsv = { "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange1_Xe_maxFOM_v2.4.0.csv", "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange2_Xe_maxFOM_v2.4.0.csv", "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange3_Xe_maxFOM_v2.4.0.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange4_Xe_maxFOM_v2.4.0.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange5_Xe_v2.4.0.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange6_Xe_maxFOM_v2.4.0.csv" };

    // v2.4.0 after revising the cuts, which appeared too strict for the 3 keV range.
    std::vector<std::string> runCutsCsv = { "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange1_Xe_maxFOM_v2.4.0_symmetricSigmaXandY.csv", "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange2_Xe_maxFOM_v2.4.0_symmetricSigmaXandY.csv", "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange3_Xe_maxFOM_v2.4.0_symmetricSigmaXandY.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange4_Xe_maxFOM_v2.4.0_symmetricSigmaXandY.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange5_Xe_v2.4.0.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange6_Xe_maxFOM_v2.4.0_symmetricSigmaXandY.csv" };


   //std::vector<std::string> runCutsCsv = { "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange1_Xe_constantEff_lowZscale.csv", "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange2_Xe_constantEff_lowZscale.csv", "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange3_Xe_constantEff_lowZscale.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange4_Xe_constantEff.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange5_Xe.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange6_Xe_constantEff.csv" };

    //std::vector<std::string> runCutsCsv = { "/home/cristina/GitHub/iaxo-quickana/sigmaCuts_Xe_100.csv", "/home/cristina/GitHub/iaxo-quickana/sigmaCuts_Xe_100.csv", "/home/cristina/GitHub/iaxo-quickana/sigmaCuts_Xe_100.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCuts_Xe_100.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCuts_Xe_100.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCuts_Xe_100.csv" };

    //std::vector<std::string> runCutsCsv = { "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange5_Xe.csv", "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange5_Xe.csv", "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange5_Xe.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange5_Xe.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange5_Xe.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange5_Xe.csv" };
    //std::vector<std::string> runCutsCsv = { "sigmaCutsRange1_Xe_constantEff_gaussZ.csv", "sigmaCutsRange2_Xe_constantEff_gaussZ.csv", "sigmaCutsRange3_Xe_constantEff_gaussZ.csv" , "sigmaCutsRange4_Xe_constantEff_gaussZ.csv" , "sigmaCutsRange5_Xe_gaussZ.csv" , "sigmaCutsRange6_Xe_constantEff_gaussZ.csv" };

     //std::vector<std::string> runCutsCsv = { "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange1_Xe_constantEff_Z100.csv", "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange2_Xe_constantEff_Z100.csv", "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange3_Xe_constantEff_Z100.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange4_Xe_constantEff_Z100.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange5_Xe_constantEff_Z100.csv" , "/home/cristina/GitHub/iaxo-quickana/sigmaCutsRange6_Xe_constantEff_Z100.csv" };
    //std::string energyCalibCsv = "/home/cristina/GitHub/iaxo-quickana/calibrationFactorsDataset3_Xe.csv";
    std::string energyCalibCsv = "/home/cristina/GitHub/iaxo-quickana/calibrationFactorsDataset3_Xe_v2.4.0.csv";

    const std::string commonPrefix = "final/SC_new"; //trackAnalysis final v2.3.15/trackAnalysis_2signalThreshold
    RunCutMap runCutMap;
    for(size_t idx = 1; idx <= 6; idx++){
        auto map = readGaussSigmaCuts(verboseAt(runCutsCsv, idx-1), idx);
        runCutMap.insert(map.begin(), map.end());
    }
    const EnergyCalibMap energyCalibMap = readEnergyCalibMap(energyCalibCsv);

    //std::vector<std::string> fileNamesCal;
    //std::vector<std::string> fileNamesBg;
    //std::vector<std::string> fileNamesFe55;


    // ******************************
    // List of Data
    // ******************************
    // Define your data here

    // Background runs
    //fileNamesBg.push_back(path + "trackAnalysis/R10512*Background*root");
    //fileNamesBg.push_back(path + "trackAnalysis/R101*Background*root");
    //fileNamesBg.push_back(path + "trackAnalysis/R10178*Background*root");
    //fileNamesBg.push_back(path + "trackAnalysis/R10182*Background*root"); // source was in
    //fileNamesBg.push_back(path + "trackAnalysis/R10188*Background*root");
    //fileNamesBg.push_back(path + "trackAnalysis/R1019[0246]*Background*root");
    /*
      fileNamesBg.push_back(path + "trackAnalysis/R102*Background*root");
      fileNamesBg.push_back(path + "trackAnalysis/R103*Background*root");
      fileNamesBg.push_back(path + "trackAnalysis/R104*Background*root");
      fileNamesBg.push_back(path + "trackAnalysis/R105*Background*root");
    */

    //Bottle 2, filter 1
    //fileNamesBg.push_back(path + commonPrefix + "/R10999*"); //Be careful!! This is a good calibration run for debugging
    //fileNamesBg.push_back(path + commonPrefix + "/R10982*"); //Be careful!! This is a bad calibration run for debugging
/*
    fileNamesBg.push_back(path + commonPrefix + "/R10850*");
    fileNamesBg.push_back(path + commonPrefix + "/R10853*");
    fileNamesBg.push_back(path + commonPrefix + "/R10857*");
    fileNamesBg.push_back(path + commonPrefix + "/R10861*");
    fileNamesBg.push_back(path + commonPrefix + "/R10863*");
    fileNamesBg.push_back(path + commonPrefix + "/R10865*");
    fileNamesBg.push_back(path + commonPrefix + "/R10867*");
    fileNamesBg.push_back(path + commonPrefix + "/R10869*");
    fileNamesBg.push_back(path + commonPrefix + "/R10871*");
    fileNamesBg.push_back(path + commonPrefix + "/R10873*");
    fileNamesBg.push_back(path + commonPrefix + "/R10875*");
    fileNamesBg.push_back(path + commonPrefix + "/R10877*");
    fileNamesBg.push_back(path + commonPrefix + "/R10879*");
    fileNamesBg.push_back(path + commonPrefix + "/R10895*");
    fileNamesBg.push_back(path + commonPrefix + "/R10897*");
    fileNamesBg.push_back(path + commonPrefix + "/R10900*");
    fileNamesBg.push_back(path + commonPrefix + "/R10902*");
    fileNamesBg.push_back(path + commonPrefix + "/R10904*");
    fileNamesBg.push_back(path + commonPrefix + "/R10906*");
    fileNamesBg.push_back(path + commonPrefix + "/R10908*");
    fileNamesBg.push_back(path + commonPrefix + "/R10910*");
    fileNamesBg.push_back(path + commonPrefix + "/R10912*");
    fileNamesBg.push_back(path + commonPrefix + "/R10914*");
    fileNamesBg.push_back(path + commonPrefix + "/R10918*");
    fileNamesBg.push_back(path + commonPrefix + "/R10920*");
    fileNamesBg.push_back(path + commonPrefix + "/R10922*");
    fileNamesBg.push_back(path + commonPrefix + "/R10924*");
    fileNamesBg.push_back(path + commonPrefix + "/R10927*");
    fileNamesBg.push_back(path + commonPrefix + "/R10932*");
    fileNamesBg.push_back(path + commonPrefix + "/R10934*");
    fileNamesBg.push_back(path + commonPrefix + "/R10936*");
    fileNamesBg.push_back(path + commonPrefix + "/R10938*");
    fileNamesBg.push_back(path + commonPrefix + "/R10941*");
    fileNamesBg.push_back(path + commonPrefix + "/R10943*");
    fileNamesBg.push_back(path + commonPrefix + "/R10945*");
    fileNamesBg.push_back(path + commonPrefix + "/R10947*");
    fileNamesBg.push_back(path + commonPrefix + "/R10949*");
    fileNamesBg.push_back(path + commonPrefix + "/R10951*");
    fileNamesBg.push_back(path + commonPrefix + "/R10953*");
    fileNamesBg.push_back(path + commonPrefix + "/R10955*");
    fileNamesBg.push_back(path + commonPrefix + "/R10957*");
    fileNamesBg.push_back(path + commonPrefix + "/R10959*");
    fileNamesBg.push_back(path + commonPrefix + "/R10961*");
    fileNamesBg.push_back(path + commonPrefix + "/R10963*");
*/
    // ******************************
    // new data
    // ******************************
     // Bottle 1, filter 1
    /*
    fileNamesBg.push_back(path + commonPrefix + "/R10710*");	//Hit channels
    fileNamesBg.push_back(path + commonPrefix + "/R10726*");	//Hit channels
    fileNamesBg.push_back(path + commonPrefix + "/R10730*");	//Hit channels
    fileNamesBg.push_back(path + commonPrefix + "/R10740*");	//Hit channels
    */
    //fileNamesBg.push_back(path + commonPrefix + "/R10746*");
    //fileNamesBg.push_back(path + commonPrefix + "/R10750*");
    /*
    fileNamesBg.push_back(path + commonPrefix + "/R10772*");	//Hit channels
    fileNamesBg.push_back(path + commonPrefix + "/R10777*");	//Hit channels
    */
    // Bottle 2, filter 1
/*
    fileNamesBg.push_back(path + commonPrefix + "/R10850*");
    fileNamesBg.push_back(path + commonPrefix + "/R10853*");
    fileNamesBg.push_back(path + commonPrefix + "/R10857*");
    fileNamesBg.push_back(path + commonPrefix + "/R10861*");
    fileNamesBg.push_back(path + commonPrefix + "/R10863*");
    fileNamesBg.push_back(path + commonPrefix + "/R10865*");
    fileNamesBg.push_back(path + commonPrefix + "/R10867*");
    fileNamesBg.push_back(path + commonPrefix + "/R10869*");
    fileNamesBg.push_back(path + commonPrefix + "/R10871*");
    fileNamesBg.push_back(path + commonPrefix + "/R10873*");
    fileNamesBg.push_back(path + commonPrefix + "/R10875*");
    fileNamesBg.push_back(path + commonPrefix + "/R10877*");
    fileNamesBg.push_back(path + commonPrefix + "/R10879*");
    fileNamesBg.push_back(path + commonPrefix + "/R10895*");
    fileNamesBg.push_back(path + commonPrefix + "/R10897*");
    fileNamesBg.push_back(path + commonPrefix + "/R10900*");
    fileNamesBg.push_back(path + commonPrefix + "/R10902*");
    fileNamesBg.push_back(path + commonPrefix + "/R10904*");
    fileNamesBg.push_back(path + commonPrefix + "/R10906*");
    fileNamesBg.push_back(path + commonPrefix + "/R10908*");
    fileNamesBg.push_back(path + commonPrefix + "/R10910*");
    fileNamesBg.push_back(path + commonPrefix + "/R10912*");
    fileNamesBg.push_back(path + commonPrefix + "/R10914*");
    fileNamesBg.push_back(path + commonPrefix + "/R10918*");
    fileNamesBg.push_back(path + commonPrefix + "/R10920*");
    fileNamesBg.push_back(path + commonPrefix + "/R10922*");
    fileNamesBg.push_back(path + commonPrefix + "/R10924*");
    fileNamesBg.push_back(path + commonPrefix + "/R10927*");
    fileNamesBg.push_back(path + commonPrefix + "/R10932*");
    fileNamesBg.push_back(path + commonPrefix + "/R10934*");
    fileNamesBg.push_back(path + commonPrefix + "/R10936*");
    fileNamesBg.push_back(path + commonPrefix + "/R10938*");
    fileNamesBg.push_back(path + commonPrefix + "/R10941*");
    fileNamesBg.push_back(path + commonPrefix + "/R10943*");
    fileNamesBg.push_back(path + commonPrefix + "/R10945*");
    fileNamesBg.push_back(path + commonPrefix + "/R10947*");
    fileNamesBg.push_back(path + commonPrefix + "/R10949*");
    fileNamesBg.push_back(path + commonPrefix + "/R10951*");
    fileNamesBg.push_back(path + commonPrefix + "/R10953*");
    fileNamesBg.push_back(path + commonPrefix + "/R10955*");
    fileNamesBg.push_back(path + commonPrefix + "/R10957*");
    fileNamesBg.push_back(path + commonPrefix + "/R10959*");
    fileNamesBg.push_back(path + commonPrefix + "/R10961*");
    fileNamesBg.push_back(path + commonPrefix + "/R10963*");

    fileNamesBg.push_back(path + commonPrefix + "/R10965*");
    fileNamesBg.push_back(path + commonPrefix + "/R10967*");
    fileNamesBg.push_back(path + commonPrefix + "/R10969*");
    fileNamesBg.push_back(path + commonPrefix + "/R10971*");
    fileNamesBg.push_back(path + commonPrefix + "/R10973*");
    fileNamesBg.push_back(path + commonPrefix + "/R10975*");
    fileNamesBg.push_back(path + commonPrefix + "/R10977*");
    fileNamesBg.push_back(path + commonPrefix + "/R10979*");
    fileNamesBg.push_back(path + commonPrefix + "/R10981*");
    fileNamesBg.push_back(path + commonPrefix + "/R10984*");
    fileNamesBg.push_back(path + commonPrefix + "/R10987*");
    */

    // Bottle 2, filter 2

/*
    fileNamesBg.push_back(path + commonPrefix + "/R10998*");
    fileNamesBg.push_back(path + commonPrefix + "/R11003*");
    fileNamesBg.push_back(path + commonPrefix + "/R11012*");
    fileNamesBg.push_back(path + commonPrefix + "/R11014*");
    fileNamesBg.push_back(path + commonPrefix + "/R11018*");
    fileNamesBg.push_back(path + commonPrefix + "/R11021*");
    fileNamesBg.push_back(path + commonPrefix + "/R11023*");
    fileNamesBg.push_back(path + commonPrefix + "/R11025*");
    fileNamesBg.push_back(path + commonPrefix + "/R11027*");
    fileNamesBg.push_back(path + commonPrefix + "/R11029*");
    fileNamesBg.push_back(path + commonPrefix + "/R11031*");
    fileNamesBg.push_back(path + commonPrefix + "/R11033*");
    fileNamesBg.push_back(path + commonPrefix + "/R11035*");
    fileNamesBg.push_back(path + commonPrefix + "/R11038*");
    fileNamesBg.push_back(path + commonPrefix + "/R11040*");
    fileNamesBg.push_back(path + commonPrefix + "/R11042*");
    fileNamesBg.push_back(path + commonPrefix + "/R11044*");
    fileNamesBg.push_back(path + commonPrefix + "/R11046*");
    fileNamesBg.push_back(path + commonPrefix + "/R11048*");
    fileNamesBg.push_back(path + commonPrefix + "/R11050*");
    fileNamesBg.push_back(path + commonPrefix + "/R11052*");
    fileNamesBg.push_back(path + commonPrefix + "/R11054*");
*/

    // Bottle 2, filter 3, 2 ln/h
/*
    fileNamesBg.push_back(path + commonPrefix + "/R11073*");
    fileNamesBg.push_back(path + commonPrefix + "/R11078*");
    fileNamesBg.push_back(path + commonPrefix + "/R11081*");
    fileNamesBg.push_back(path + commonPrefix + "/R11083*");
    fileNamesBg.push_back(path + commonPrefix + "/R11086*");
    fileNamesBg.push_back(path + commonPrefix + "/R11089*");
    fileNamesBg.push_back(path + commonPrefix + "/R11091*");
    fileNamesBg.push_back(path + commonPrefix + "/R11093*");
*/

    // Bottle 2, filter 3, 5 ln/h

    /*
    fileNamesBg.push_back(path + commonPrefix + "/R11108*");
    fileNamesBg.push_back(path + commonPrefix + "/R11110*");
    fileNamesBg.push_back(path + commonPrefix + "/R11113*");
    fileNamesBg.push_back(path + commonPrefix + "/R11115*");
    fileNamesBg.push_back(path + commonPrefix + "/R11117*");
    fileNamesBg.push_back(path + commonPrefix + "/R11120*");
    fileNamesBg.push_back(path + commonPrefix + "/R11132*");
    fileNamesBg.push_back(path + commonPrefix + "/R11134*");
    fileNamesBg.push_back(path + commonPrefix + "/R11149*");
    fileNamesBg.push_back(path + commonPrefix + "/R11153*");
    fileNamesBg.push_back(path + commonPrefix + "/R11205*");
    fileNamesBg.push_back(path + commonPrefix + "/R11208*");
    fileNamesBg.push_back(path + commonPrefix + "/R11210*");
    fileNamesBg.push_back(path + commonPrefix + "/R11218*");
    fileNamesBg.push_back(path + commonPrefix + "/R11223*");
    fileNamesBg.push_back(path + commonPrefix + "/R11225*");
    fileNamesBg.push_back(path + commonPrefix + "/R11227*");
    fileNamesBg.push_back(path + commonPrefix + "/R11230*");
    fileNamesBg.push_back(path + commonPrefix + "/R11241*");
    fileNamesBg.push_back(path + commonPrefix + "/R11245*");
    fileNamesBg.push_back(path + commonPrefix + "/R11250*");
    fileNamesBg.push_back(path + commonPrefix + "/R11252*");
    fileNamesBg.push_back(path + commonPrefix + "/R11254*");
    fileNamesBg.push_back(path + commonPrefix + "/R11281*");
    //subrun 4 is not readable! Why?
    fileNamesBg.push_back(path + commonPrefix + "/R11284_0000[0123]*");
    fileNamesBg.push_back(path + commonPrefix + "/R11286*");
    fileNamesBg.push_back(path + commonPrefix + "/R11289*");
    fileNamesBg.push_back(path + commonPrefix + "/R11293*");
    fileNamesBg.push_back(path + commonPrefix + "/R11296*");
    fileNamesBg.push_back(path + commonPrefix + "/R11298*");
    fileNamesBg.push_back(path + commonPrefix + "/R11300*");
    fileNamesBg.push_back(path + commonPrefix + "/R11304*");
    fileNamesBg.push_back(path + commonPrefix + "/R11306*");
    fileNamesBg.push_back(path + commonPrefix + "/R11308*");
    fileNamesBg.push_back(path + commonPrefix + "/R11310*");
    fileNamesBg.push_back(path + commonPrefix + "/R11313*");
    fileNamesBg.push_back(path + commonPrefix + "/R11315*");
    fileNamesBg.push_back(path + commonPrefix + "/R11317*");
    fileNamesBg.push_back(path + commonPrefix + "/R11319*");
    fileNamesBg.push_back(path + commonPrefix + "/R11321*");


    // Bottle 2, filter 3, small volume

    fileNamesBg.push_back(path + commonPrefix + "/R11339*");
    fileNamesBg.push_back(path + commonPrefix + "/R11341*");
    fileNamesBg.push_back(path + commonPrefix + "/R11343*");
    fileNamesBg.push_back(path + commonPrefix + "/R11344*");
    fileNamesBg.push_back(path + commonPrefix + "/R11347*");
    fileNamesBg.push_back(path + commonPrefix + "/R11353*");
    fileNamesBg.push_back(path + commonPrefix + "/R11358*");
    fileNamesBg.push_back(path + commonPrefix + "/R11361*");
    fileNamesBg.push_back(path + commonPrefix + "/R11367*");
    fileNamesBg.push_back(path + commonPrefix + "/R11375*");
    fileNamesBg.push_back(path + commonPrefix + "/R11378*");
    fileNamesBg.push_back(path + commonPrefix + "/R11381*");
    fileNamesBg.push_back(path + commonPrefix + "/R11384*");
    fileNamesBg.push_back(path + commonPrefix + "/R11393*");
    fileNamesBg.push_back(path + commonPrefix + "/R11398*");
    fileNamesBg.push_back(path + commonPrefix + "/R11402*");
    fileNamesBg.push_back(path + commonPrefix + "/R11404*");
    fileNamesBg.push_back(path + commonPrefix + "/R11406*");
    fileNamesBg.push_back(path + commonPrefix + "/R11408*");


    // Bottle 2, filter 3, large volume

    fileNamesBg.push_back(path + commonPrefix + "/R11411*");
    fileNamesBg.push_back(path + commonPrefix + "/R11419*");
    fileNamesBg.push_back(path + commonPrefix + "/R11422*");
    fileNamesBg.push_back(path + commonPrefix + "/R11424*");
    fileNamesBg.push_back(path + commonPrefix + "/R11426*");
    fileNamesBg.push_back(path + commonPrefix + "/R11430*");
    fileNamesBg.push_back(path + commonPrefix + "/R11433*");
    fileNamesBg.push_back(path + commonPrefix + "/R11435*");
    fileNamesBg.push_back(path + commonPrefix + "/R11437*");
    fileNamesBg.push_back(path + commonPrefix + "/R11439*");
    fileNamesBg.push_back(path + commonPrefix + "/R11456*");
    fileNamesBg.push_back(path + commonPrefix + "/R11459*");
    fileNamesBg.push_back(path + commonPrefix + "/R11461*");
    fileNamesBg.push_back(path + commonPrefix + "/R11463*");
    fileNamesBg.push_back(path + commonPrefix + "/R11467*");
    fileNamesBg.push_back(path + commonPrefix + "/R11470*");
    fileNamesBg.push_back(path + commonPrefix + "/R11472*");
    fileNamesBg.push_back(path + commonPrefix + "/R11476*");
    fileNamesBg.push_back(path + commonPrefix + "/R11478*");
    fileNamesBg.push_back(path + commonPrefix + "/R11510*");
    fileNamesBg.push_back(path + commonPrefix + "/R11517*");
    fileNamesBg.push_back(path + commonPrefix + "/R11519*");
    fileNamesBg.push_back(path + commonPrefix + "/R11521*");
    fileNamesBg.push_back(path + commonPrefix + "/R11523*");
    fileNamesBg.push_back(path + commonPrefix + "/R11528*");
    fileNamesBg.push_back(path + commonPrefix + "/R11533*");
    fileNamesBg.push_back(path + commonPrefix + "/R11535*");
    fileNamesBg.push_back(path + commonPrefix + "/R11538*");
    fileNamesBg.push_back(path + commonPrefix + "/R11540*");
    fileNamesBg.push_back(path + commonPrefix + "/R11542*");
    fileNamesBg.push_back(path + commonPrefix + "/R11544*");
    fileNamesBg.push_back(path + commonPrefix + "/R11546*");
    fileNamesBg.push_back(path + commonPrefix + "/R11549*");
    fileNamesBg.push_back(path + commonPrefix + "/R11552*");
    fileNamesBg.push_back(path + commonPrefix + "/R11554*");
    fileNamesBg.push_back(path + commonPrefix + "/R11556*");
    fileNamesBg.push_back(path + commonPrefix + "/R11560*");
    fileNamesBg.push_back(path + commonPrefix + "/R11563*");
    fileNamesBg.push_back(path + commonPrefix + "/R11565*");
    fileNamesBg.push_back(path + commonPrefix + "/R11567*");
*/

    // ******************************
    // Create RDataFrames
    // ******************************

    // ROOT::RDataFrame bgData("AnalysisTree",fileNamesBg);
    // ROOT::RDataFrame calData("AnalysisTree",fileNamesCal);
    // ROOT::RDataFrame Fe55Data("AnalysisTree",fileNamesFe55);

   // this is the file we store our snapshot dataframes in!
    std::ofstream outfile("dataframe_snapshots_55FeCuts_Xe_v2.4.0_withSCdata.txt", std::ofstream::app);
    std::ofstream outfileTime("totalTime_list_55FeCuts_Xe_v2.4.0_withSCdata.txt", std::ofstream::app);
    int totalDuration = 0;
    auto snapshotDurationTuple = applyBackgroundCuts(path+commonPrefix+runPrefix, commonPrefix, runCutMap, energyCalibMap);
    auto snapshotName = snapshotDurationTuple.first;
    totalDuration += snapshotDurationTuple.second;
    if(snapshotName.size() > 0){
        outfile << snapshotName << std::endl;
        outfileTime << totalDuration << std::endl;
    }
    outfile.close();
    outfileTime.close();
    return;
}

std::vector<std::string> parseFile(std::string fname){
    std::vector<std::string> result;
    std::ifstream file(fname.c_str());
    std::string line;
    if(file.is_open()){
        while(std::getline(file, line, '\n')){
            result.push_back(line);
        }
    }
    file.close();
    return result;
}

std::vector<double> asDouble(std::vector<std::string> vec){
    std::vector<double> result;
    for(const auto s : vec){
        result.push_back(std::stod(s));
    }
    return result;
}

double sum(std::vector<double> vec){
    double result = 0;
    for(double s : vec){
        result += s;
    }
    return result;
}

void printBackgroundRate(ROOT::RDataFrame bgData, double duration){
    double area = 36;
    double kev = 12.0;
    auto nEntriesAfterCuts = bgData.Count().GetValue();
    std::cout << "Number of events that pass the cuts: " << nEntriesAfterCuts << std::endl;
    std::cout << "Duration =  " << duration << " s" << " || " << duration / 3600 / 24 << " days " << std::endl;
    auto bgRate = nEntriesAfterCuts / kev / area / duration;
    std::cout << "Background rate = " << bgRate << " ± " << sqrt(nEntriesAfterCuts) /kev / area / duration << " counts/keV/cm²/s between 0 and 12 keV for the full readout area" << std::endl;
}

void analyze(std::string snapshotsFile, std::string durationFile ){
    // parse the input file containing file names that are DF snapshots
    auto snapshots = parseFile(snapshotsFile);
    auto runDuration = asDouble(parseFile(durationFile));
    for( int i = 0; i<snapshots.size(); i++){
        ROOT::RDataFrame bgData("AnalysisTree", snapshots[i]);
        printBackgroundRate(bgData, runDuration[i]);
    }
    ROOT::RDataFrame bgData("AnalysisTree", snapshots);
    printBackgroundRate(bgData, sum(runDuration));
    bgData.Snapshot("AnalysisTree", "background_df_all_runs_Dataset3_Xe_extrapolatedCuts_newSigmas_symmetricSigmaXandY_v2.4.0_SCdata.root");


    //Make some plots and print results

    auto histo1 = bgData.Histo1D({"energy_keV", "Energy background", 48, 0, 12}, "energy_keV");
    TCanvas* c01 = new TCanvas();
    histo1->DrawClone();

    // auto calData2 = applyDefinesAndFilters(calData, energyIdx, dkXrayTube);
    // auto Fe55Data2 = applyDefinesAndFilters(Fe55Data, energyIdx, dkBackground);
    //
    //
    // ************************************************************
    // Estimate run duration
    // ************************************************************

    // double binw = 60;
    // int rate_threshold = 10;  // counts per minute
    // int run_duration = GetRunDuration(bgData, binw, rate_threshold);
    // int run_duration_combined = 0;
    // run_duration_combined += run_duration;

    // cout << "\n\033[1;35mrun duration: " << run_duration << " s = " << (double)run_duration / 3600
    //      << " h = " << (double)run_duration / (3600 * 24) << " d" << endl;

    //
    // // ************************************************************
    // // CUTS: generate and apply the cut strings
    // // ************************************************************
    //
    //
    // // apply cuts
    // auto cutBackground = genCutString(energyIdx, dkBackground, cutScaling);
    // auto cutXrayTube = genCutString(energyIdx, dkXrayTube, cutScaling);
    // auto cutFe55 = genCutString(energyIdx, dkBackground, cutScaling);
    // cout << "\nCuts applied to background data: " << cutBackground << endl;
    // cout << "\nCuts applied to Fe55 data: " << cutFe55 << endl;
    // cout << "\nCuts applied to X-ray tube data: " << cutXrayTube << endl;
    //
    // double n0_bg = bgData2.Count().GetValue();
    // double n0_cal = calData2.Count().GetValue();
    // double n0_fe = Fe55Data2.Count().GetValue();
    //
    // cout << "\n\n\033[1;35m# events in bg DF before: " << n0_bg  << endl;
    // cout << "# events in cal DF before: " << n0_cal  << endl;
    // cout << "# events in Fe55 DF before: " << n0_fe  << endl;
    //
    // bgData2 = bgData2.Filter(cutBackground);
    // Fe55Data2 = Fe55Data2.Filter(cutBackground);
    // calData2 = calData2.Filter(cutXrayTube);
    //
    // double n1_cal;
    // double n1_fe;
    // double n1_bg;
    //
    //
    // n1_bg = bgData2.Count().GetValue();
    // n1_cal = calData2.Count().GetValue();
    // n1_fe = Fe55Data2.Count().GetValue();
    //
    //
    // double FOM0_tube = n0_cal/sqrt(n0_bg);
    // double FOM1_tube = n1_cal/sqrt(n1_bg);
    //
    // double FOM0_Fe55 = n0_fe/sqrt(n0_bg);
    // double FOM1_Fe55 = n1_fe/sqrt(n1_bg);
/*
// Print columns' names
auto colNames = bgData.GetColumnNames();
for (auto &&colName : colNames) std::cout << colName << std::endl;
// Print column type
auto colType = bgData.GetColumnType("veto_MaxPeakAmplitude");
std::cout << "Column veto_MaxPeakAmplitude has type " << colType << std::endl;
*/

    // Draw some histograms of the filtered data

    // auto histo0 = calData2.Histo1D({"tckAna_MaxTrackEnergy", "Energy X-ray tube", 100, 0, 13000}, "tckAna_MaxTrackEnergy");
    // TCanvas* c00 = new TCanvas();
    // histo0->DrawClone();
    //
    // auto histo1 = bgData2.Histo1D({"tckAna_MaxTrackEnergy", "Energy background", 48, 0, 13000}, "tckAna_MaxTrackEnergy");
    // TCanvas* c01 = new TCanvas();
    // //TCanvas c("c","x hist");
    // histo1->DrawClone();
    //
    // auto histo2 = Fe55Data2.Histo1D({"tckAna_MaxTrackEnergy", "Energy Fe55 calibration", 100, 0, 13000}, "tckAna_MaxTrackEnergy");
    // TCanvas* c02 = new TCanvas();
    // histo2->DrawClone();


/*
  auto histo2 = bgData_filtered.Histo1D("eventID");
  TCanvas* c02 = new TCanvas();
  //TCanvas c("c","x hist");
  histo2->DrawClone();

  auto histo3 = bgData_filtered.Histo1D({"timeStamp", "timeStamp", 10000, 0, 17e9}, "timeStamp");
  TCanvas* c03 = new TCanvas();
  //TCanvas c("c","x hist");
  histo3->DrawClone();

  auto histo4 = bgData_filtered.Histo1D("veto_MaxPeakAmplitude.second");
  TCanvas* c04 = new TCanvas();
  //TCanvas c("c","x hist");
  histo4->DrawClone();
*/

    // Print the results in terminal.
    // cout << "\n# events in bg DF after: " <<  n1_bg << endl;
    // cout << "# events in cal DF after: " <<  n1_cal << endl;
    // cout << "# events in Fe55 DF after: " <<  n1_fe << endl;
    // cout << "\nBackground: cut efficiency = " << 100 * n1_bg / n0_bg << " %" << endl;
    // cout << "Calibration: cut efficiency = " << 100 * n1_cal / n0_cal << " %" << endl;
    // cout << "Fe55: cut efficiency = " << 100 * n1_fe / n0_fe << " %" << endl;
    // cout << "\nX-ray tube calibrations" << endl;
    // cout << "\nFOM before cut = " << FOM0_tube << endl;
    // cout << "FOM after cut = " << FOM1_tube << endl;
    // cout << "FOM after cut as eff/eff_bck = " << (100 * n1_cal / n0_cal)/(sqrt(100 * n1_bg / n0_bg )) << endl;
    // cout << "\nFe55 calibrations" << endl;
    // cout << "\nFOM before cut = " << FOM0_Fe55 << endl;
    // cout << "FOM after cut = " << FOM1_Fe55 << endl;
    // cout << "FOM after cut as eff/eff_bck = " << (100 * n1_fe / n0_fe)/(sqrt(100 * n1_bg / n0_bg )) << endl;
    // cout << "\n\n\033[0m";

}

int main(int argc, char **argv) {
    gSystem->Load("/programas/rest/v2.4.0/lib/libRestConnectors.so");
    gSystem->Load("/programas/rest/v2.4.0/lib/libRestDetector.so");
    gSystem->Load("/programas/rest/v2.4.0/lib/libRestFramework.so");
    //gSystem->Load("/programas/rest/v2.4.0/lib/libRestGeant4.so");
    gSystem->Load("/programas/rest/v2.4.0/lib/libRestLegacy.so");
    gSystem->Load("/programas/rest/v2.4.0/lib/libRestRaw.so");
    //gSystem->Load("/programas/rest/v2.4.0/lib/libRestRestG4.so");
    gSystem->Load("/programas/rest/v2.4.0/lib/libRestTrack.so");
    //gSystem->Load("/programas/rest/v2.4.0/lib/libRestWimp.so");

    if(argc != 2){
	std::cout << "Please hand exactly one argument, the glob describing the runs to work on" << std::endl;
	return -1;
    }
    std::string fname = argv[1];
    cut(fname);
    return 0;
}
