// Author: Cristina Margalejo and Sebastian Schmidt, November 2022
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
// #include "ROOT/RVec.hxx"

using namespace ROOT::VecOps;

// custom type alias for a bit easier reading of the map mapping
// each run number to the correct cuts for X / Y gauss sigma.
typedef std::map<std::pair<int, int>, std::pair<std::pair<double, double>, std::pair<double, double>>> RunCutMap;
//               ^--- run number
//                                        ^---  & ^-- X cuts         ^---  ^ ^--- Y cuts
typedef std::map<int, int> EnergyCalibMap;
typedef std::vector<std::vector<std::pair<double, double>>> EnergyCutsVec;
//           ^--- energy index
//                       ^--- observable in specific order
//                                        ^---  & ^--- cut low & high

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//TODO: write documentation!
// Macro to produce a dataframe keeping only the events that pass the cuts.
//
// Functions: defined the at the beginning of each function.
//
//Usage:
//restRoot
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
    vector<double> timestamps = result_ts.GetValue();

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
    cout << "run_duration = " << run_duration << endl;
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
        { "tckAna_MaxTrackEnergyBalanceXY", {
                { 1,  { 2.2, 2.2 }  },
                { 2,  { 1.8, 1.8 }  },
                { 3,  { 1.42, 1.42 }  },
                { 4,  { 1.1, 1.1 }  },
                { 5,  { 1.0, 1.0 }  },
                { 6,  { 0.95, 0.95 }  },
          }
        },
        { "tckAna_MaxTrack_XZ_GaussSigmaX", {
                { 1,  { 0.78, 1.36 }  },
                { 2,  { 0.78, 1.18 }  },
                { 3,  { 0.78, 1.06 }  },
                { 4,  { 0.78, 1.0 }  },
                { 5,  { 1.0, 1.0 }  },
                { 6,  { 1.22, 0.94 }  },
          }
        },
        { "tckAna_MaxTrack_YZ_GaussSigmaY", {
                { 1,  { 0.78, 1.36 }  },
                { 2,  { 0.78, 1.18 }  },
                { 3,  { 0.78, 1.06 }  },
                { 4,  { 0.78, 1.0 }  },
                { 5,  { 1.0, 1.0 }  },
                { 6,  { 1.22, 0.94 }  },
          }
        },
        /*
        { "hitsAna_z2Sigma", {
                { 1,  { 0.0, 10.5 }  }, //7 before revising the cuts
                { 2,  { 0.0, 7.0 }  }, //5
                { 3,  { 0.0, 5.8 }  }, //5
                { 4,  { 0.0, 3.8 }  },
                { 5,  { 0.0, 1.0 }  },
                { 6,  { 0.0, 1.6 }  },
          }
        },
        { "tckAna_MaxTrack_XZ_GaussSigmaZ", {
                { 1,  { 0.0, 10.5 }  }, //7 before revising the cuts
                { 2,  { 0.0, 7.0 }  }, //5
                { 3,  { 0.0, 5.8 }  }, //5
                { 4,  { 0.0, 3.8 }  },
                { 5,  { 0.0, 1.0 }  },
                { 6,  { 0.0, 1.6 }  },
          }
        },
        { "tckAna_MaxTrack_YZ_GaussSigmaZ", {
                { 1,  { 0.0, 10.5 }  }, //7 before revising the cuts
                { 2,  { 0.0, 7.0 }  }, //5
                { 3,  { 0.0, 5.8 }  }, //5
                { 4,  { 0.0, 3.8 }  },
                { 5,  { 0.0, 1.0 }  },
                { 6,  { 0.0, 1.6 }  },
          }
        },
        */
        { "tckAna_MaxTrackxySigmaGausBalance", {
                { 1,  { 1.0, 1.0 }  },
                { 2,  { 1.0, 1.0 }  },
                { 3,  { 1.0, 1.0 }  },
                { 4,  { 1.0, 1.0 }  },
                { 5,  { 1.0, 1.0 }  },
                { 6,  { 1.0, 1.0 }  },
          }
        },
        { "hitsAna_xySkew", {
                { 1,  { 1.92, 1.95 }  },
                { 2,  { 1.62, 1.77 }  },
                { 3,  { 1.62, 1.92 }  },
                { 4,  { 1.15, 1.15 }  },
                { 5,  { 1.0, 1.0 }  },
                { 6,  { 0.85, 0.85 }  },
          }
        },
        { "hitsAna_zSkew", {
                { 1,  { 0.54, 0.56 }  },
                { 2,  { 0.61, 0.72 }  },
                { 3,  { 0.72, 1.02 }  },
                { 4,  { 0.63, 0.86 }  },
                { 5,  { 1.0, 1.0 }  },
                { 6,  { 0.72, 1.18 }  },
          }
        },
    };
    return map[obs][energyIdx];
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
        { "tckAna_MaxTrackEnergyBalanceXY",  { -0.18, 0.18 }  },
        { "tckAna_MaxTrack_XZ_GaussSigmaX",  { 0.35, 1.5 }  },
        { "tckAna_MaxTrack_YZ_GaussSigmaY",  { 0.4, 1.5 }  },
        //{ "hitsAna_z2Sigma",  { 0.0, 1000000 }  },
        //{ "tckAna_MaxTrack_XZ_GaussSigmaZ",  { 0.0, 2 }  },
        //{ "tckAna_MaxTrack_YZ_GaussSigmaZ",  { 0.0, 2 }  },
        { "tckAna_MaxTrackxySigmaGausBalance",  { -0.55, 0.55 }  },
        { "hitsAna_xySkew",  { -0.85, 0.85 }  },
        { "hitsAna_zSkew",  { -2.5, 2.5 }  },
    };
     std::map<std::string, std::pair<double, double>> xrayTubeMap = {
        { "tckAna_MaxTrackEnergyBalanceXY",  { -0.21, 0.21 }  },
        { "tckAna_MaxTrack_XZ_GaussSigmaX",  { 0.225, 0.775 }  },
        { "tckAna_MaxTrack_YZ_GaussSigmaY",  { 0.225, 0.825 }  },
        //{ "hitsAna_z2Sigma",  { 0.0, 1000000 }  },
        //{ "tckAna_MaxTrack_XZ_GaussSigmaZ",  { 0.0, 1.25 }  },
        //{ "tckAna_MaxTrack_YZ_GaussSigmaZ",  { 0.0, 1.25 }  },
        { "tckAna_MaxTrackxySigmaGausBalance",  { -0.5, 0.5 }  },
        { "hitsAna_xySkew",  { -0.52, 0.52 }  },
        { "hitsAna_zSkew",  { -3.5, 3.5 }  },
    };

    switch(dKind) {
        case dkBackground:
            return castMap[obs];
        case dkXrayTube:
            return xrayTubeMap[obs];
        default:
            assert(false);
            return std::make_pair(0.0, 0.0);
    };
}
//***************************************************************
// Function to get the cut scaling applied to each energy range.
//***************************************************************
double getCutScaling(int energyIdx){
    //std::vector<double> vec = { 0.2, 0.18, 0.15, 0.15, 0.0, 0.05 }; //v2.3.X
    // std::vector<double> vec = { 0.2, 0.18, 0.0, 0.0, 0.0, 0.0 }; //V2.4 dataset 1
    std::vector<double> vec = { 0.2, 0.0, 0.0, 0.0, 0.0, 0.0 }; //V2.4 dataset 2
    return vec[energyIdx-1];
}

EnergyCutsVec cutVec; // global variable.
//***************************************************************
// Function to decide if an event pass the cuts or not.
//***************************************************************
bool applyCut(double energyBalance, double gaussSigmaX, double gaussSigmaY,
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
    if(energyBalance            < vec[0].first || energyBalance         > vec[0].second){
        return false;
    }
    else if(gaussSigmaX         < vec[1].first || gaussSigmaX           > vec[1].second){
        return false;
    }
    else if(gaussSigmaY         < vec[2].first || gaussSigmaY           > vec[2].second){
        return false;
    }
    /*
    else if(z2Sigma             < vec[3].first || z2Sigma               > vec[3].second){
        return false;
    }
    else if(gaussSigmaXZ        < vec[4].first || gaussSigmaXZ          > vec[4].second){
        return false;
    }
    else if(gaussSigmaYZ        < vec[5].first || gaussSigmaYZ          > vec[5].second){
        return false;
    }
    */
    else if(sigmaGaussBalance   < vec[3].first || sigmaGaussBalance     > vec[3].second){
        return false;
    }
    else if(xySkew              < vec[4].first || xySkew                > vec[4].second){
        return false;
    }
    else if(zSkew               < vec[5].first || zSkew                 > vec[5].second){
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
        //"hitsAna_z2Sigma",
        //"tckAna_MaxTrack_XZ_GaussSigmaZ",
        //"tckAna_MaxTrack_YZ_GaussSigmaZ",
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
            /*
            if(idx != 1 && idx != 2){
                energyVec.push_back(cut);
            }
            else if(idx == 2){
                auto gaussSigmaCuts = csvData[std::make_pair(run, energyIdx)];
                energyVec.push_back( gaussSigmaCuts.first );  // corresponds to X cuts
                energyVec.push_back( gaussSigmaCuts.second ); // corresponds to Y cuts
            }
            */
            energyVec.push_back(cut);
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
                }
                col++;
            }
            result[std::make_pair(run, energyIdx)] = std::make_pair( x, y );
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
    const double calibFactor = (double)energyCalibMap[runNumber];
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
    result = result.Filter("eventID > 0");

    // tracks cut
    //result = result.Filter("tckAna_nTracks_X == 1");
    //result = result.Filter("tckAna_nTracks_Y == 1");
    //1+1, 1+2, 2+1
    result = result.Filter("tckAna_nTracks_X >= 1 && tckAna_nTracks_Y >= 1 && (tckAna_nTracks_X+tckAna_nTracks_Y<=3)");
    //calData_filtered = calData_filtered.Filter("tckAna_nTracks_X >= 1 && tckAna_nTracks_Y >= 1 && (tckAna_nTracks_X+tckAna_nTracks_Y<=3)");

    // energy cut
    //result = result.Filter("energy_keV <= 12");
    //result = result.Filter("hitsAna_xMean*hitsAna_xMean+hitsAna_yMean*hitsAna_yMean < 100 ");

    result.Display()->Print();

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
    cout << "\n\033[1;35mrun duration: " << run_duration << " s = " << (double)run_duration / 3600
         << " h = " << (double)run_duration / (3600 * 24) << " d" << endl;
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
                               //"hitsAna_z2Sigma",
                               //"tckAna_MaxTrack_XZ_GaussSigmaZ",
                               //"tckAna_MaxTrack_YZ_GaussSigmaZ",
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
        result = "background_df_run_" + std::to_string(runNumber) + ".root";
        bgData.Snapshot("AnalysisTree", result);
    }
    // bgData.Display({"hitsAna_z2Sigma"}, 100)->Print();
    std::cout << "Number of events left after cuts in run " << runNumber << " = " << *cnt << std::endl;
    return std::make_pair(result,run_duration);
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

void cut() {

    // ******************************
    // Important variables and declarations
    // ******************************

    string processName_tckAna = "tckAna_";                       // process prefix
    string processName_hitsAna = "hitsAna_";                         // process prefix
    string path = "/storage/cast/SR2019/analysis/argon/official/v2.4.0/";  // data path

    // Introduce the filenames

    std::string energyCalibCsv = "/home/cristina/GitHub/iaxo-quickana/calibrationFactorsDataset2_v2.4.0.csv";

    const std::string commonPrefix = "final"; //trackAnalysis

    RunCutMap runCutMap;
    /*
    for(size_t idx = 1; idx <= 6; idx++){
        auto map = readGaussSigmaCuts(runCutsCsv[idx-1], idx);
        runCutMap.insert(map.begin(), map.end());
    }
    */
    const EnergyCalibMap energyCalibMap = readEnergyCalibMap(energyCalibCsv);

    vector<string> fileNamesCal;
    vector<string> fileNamesBg;
    vector<string> fileNamesFe55;


    // ******************************
    // List of Data
    // ******************************
    // Define your data here

    // Background runs
/*

    fileNamesBg.push_back(path + "/R10178*");
    fileNamesBg.push_back(path + "/R10182*");
    fileNamesBg.push_back(path + "/R10188*");
    fileNamesBg.push_back(path + "/R10190*");
    fileNamesBg.push_back(path + "/R10192*");
    fileNamesBg.push_back(path + "/R10194*");
    fileNamesBg.push_back(path + "/R10196*");
*/
    // Dataset 2
    /*
    fileNamesBg.push_back(path + commonPrefix + "/R10648*");
    fileNamesBg.push_back(path + commonPrefix + "/R10660*");
    fileNamesBg.push_back(path + commonPrefix + "/R10663*");
    fileNamesBg.push_back(path + commonPrefix + "/R10665*");*/

    fileNamesBg.push_back(path + commonPrefix + "/R10671*");
    fileNamesBg.push_back(path + commonPrefix + "/R10690*");
    fileNamesBg.push_back(path + commonPrefix + "/R10691*");
    fileNamesBg.push_back(path + commonPrefix + "/R10781*");
    fileNamesBg.push_back(path + commonPrefix + "/R10793*");
    fileNamesBg.push_back(path + commonPrefix + "/R10795*");
    fileNamesBg.push_back(path + commonPrefix + "/R10796*");
    fileNamesBg.push_back(path + commonPrefix + "/R10800*");
    fileNamesBg.push_back(path + commonPrefix + "/R10802*");
    fileNamesBg.push_back(path + commonPrefix + "/R10805*");
    fileNamesBg.push_back(path + commonPrefix + "/R10807*");
    fileNamesBg.push_back(path + commonPrefix + "/R10811*");
    fileNamesBg.push_back(path + commonPrefix + "/R10822*");
    fileNamesBg.push_back(path + commonPrefix + "/R10822*");
    fileNamesBg.push_back(path + commonPrefix + "/R10824*");
    fileNamesBg.push_back(path + commonPrefix + "/R10826*");
    fileNamesBg.push_back(path + commonPrefix + "/R10828*");
    fileNamesBg.push_back(path + commonPrefix + "/R10830*");
    fileNamesBg.push_back(path + commonPrefix + "/R10832*");
    fileNamesBg.push_back(path + commonPrefix + "/R10834*");
    fileNamesBg.push_back(path + commonPrefix + "/R10836*");
    fileNamesBg.push_back(path + commonPrefix + "/R10838*");
    fileNamesBg.push_back(path + commonPrefix + "/R10840*");
    fileNamesBg.push_back(path + commonPrefix + "/R10842*");
    fileNamesBg.push_back(path + commonPrefix + "/R10843*");


    // ******************************
    // Create RDataFrames
    // ******************************

    // ROOT::RDataFrame bgData("AnalysisTree",fileNamesBg);
    // ROOT::RDataFrame calData("AnalysisTree",fileNamesCal);
    // ROOT::RDataFrame Fe55Data("AnalysisTree",fileNamesFe55);

    std::vector<std::string> snapshots;
    int total_duration = 0;
    for(auto fname : fileNamesBg){
        //auto snapshotDurationTuple = applyBackgroundCuts(fname, commonPrefix, runCutMap, energyCalibMap);
        auto snapshotDurationTuple = applyBackgroundCuts(fname, commonPrefix, runCutMap, energyCalibMap);
	auto snapshotName = snapshotDurationTuple.first;
	total_duration += snapshotDurationTuple.second;
        if(snapshotName.size() > 0){
            snapshots.push_back(snapshotName);
        }
    }

    ROOT::RDataFrame bgData("AnalysisTree", snapshots);
    bgData.Snapshot("AnalysisTree", "background_df_all_runs_dataset2_subsets2and3_v2.4.0_extrapolatedCuts_withRateAna.root"); //background_df_all_runs_dataset2_subset3_calibrated_3tracks_v2.4.0.root

    //Make some plots and print results
    double area = 36;
    double kev = 12.0;
    auto nEntriesAfterCuts = bgData.Count().GetValue();
    std::cout << "Number of events that pass the cuts: " << nEntriesAfterCuts << std::endl;
    std::cout << "Duration =  " << total_duration << " s" << " || " << total_duration / 3600 / 24 << " days " << std::endl;
    auto bgRate = nEntriesAfterCuts / kev / area / total_duration;
    std::cout << "Background rate = " << bgRate << " ± " << sqrt(nEntriesAfterCuts) /kev / area / total_duration << " counts/keV/cm²/s between 0 and 12 keV for the full readout area" << endl;
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
auto colNames = bgData_filtered.GetColumnNames();
for (auto &&colName : colNames) std::cout << colName << std::endl;
// Print column type
auto colType = bgData_filtered.GetColumnType("veto_MaxPeakAmplitude");
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
