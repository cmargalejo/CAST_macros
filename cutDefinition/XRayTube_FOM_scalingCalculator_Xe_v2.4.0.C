

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

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//TODO: write documentation!
// Macro to define one or several cuts and calculate the resulting FOM.
//
// Purpose: main idea-> start from predefined cuts, apply scaling...
// Function: define the at the beginning of each function.
//
//
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
    TCanvas* c00 = new TCanvas();
    ht->Draw();
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
// Function used to calibrate the 6 energy ranges.
//****************************************************************
int getCalibrationFactor(int energyIdx) {
    switch(energyIdx) {
        case 1: return 438;
        case 2: return 476;
        case 3: return 497;
        case 4: return 504;
        case 5: return 489;
        case 6: return 500;
        default:
            std::cout << "Invalid energy range given: `" << energyIdx
                      << "`. Only values in [1, 6] allowed" << std::endl;
            assert(false);
            return 0;
    }
}

//***************************************************************
// Function to get the energy range. The argument energyIdx is given by the
// user in the main cut() function.
//***************************************************************
std::pair<double, double> getEnergyRange(int energyIdx){
    switch(energyIdx) {
        case 1: return std::make_pair(0.40, 1.75);
        case 2: return std::make_pair(1.75, 2.50);
        case 3: return std::make_pair(2.50, 3.75);
        case 4: return std::make_pair(3.75, 5.25);
        case 5: return std::make_pair(5.25, 7.00);
        case 6: return std::make_pair(7.00, 12.0);
        default:
            std::cout << "Invalid energy range given: `" << energyIdx
                      << "`. Only values in [1, 6] allowed" << std::endl;
            assert(false);
            return std::make_pair(0.0, 0.0);
    }
}

//***************************************************************
// Function to get the calibrated enery, i.e., to turn ADC into keV.
//***************************************************************
ROOT::RDF::RNode calcEnergy(ROOT::RDF::RNode df, int energyIdx, DataKind dKind){
    switch(dKind) {
        case dkBackground:
            return df.Define("energy_keV","tckAna_MaxTrackEnergy/460");
        case dkXrayTube:
            auto calcStr = "tckAna_MaxTrackEnergy/" + to_string(getCalibrationFactor(energyIdx));
            return df.Define("energy_keV", calcStr);
    }
}

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
    return map[obs][energyIdx];
}

//***************************************************************
// Function to access the inital values of the cuts, both for the 55Fe
// calibrations at CAST and for the 5.9 keV calibration in the X-ray tube.
//***************************************************************
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
    std::map<std::string, std::pair<double, double>> xrayTubeMap = {
        { "tckAna_MaxTrackEnergyBalanceXY",  { -0.20, 0.15 }  },
        { "tckAna_MaxTrack_XZ_GaussSigmaX",  { 0.39, 1.79 }  },
        { "tckAna_MaxTrack_YZ_GaussSigmaY",  { 0.41, 1.85 }  },
        { "hitsAna_z2Sigma",  { 0.0, 1.45 }  },
        { "tckAna_MaxTrackxySigmaGausBalance",  { -0.28, 0.68 }  },
        { "hitsAna_xySkew",  { -0.42, 0.82 }  },
        { "hitsAna_zSkew",  { -10, 10 }  },
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

//*****************************************************************
// Function to generate the cut string.
// This cut string will be the input to the main cut() function.
//*****************************************************************
std::string genCutString(int energyIdx, DataKind dKind, double cutScaling){
    const auto names = {
        "tckAna_MaxTrackEnergyBalanceXY",
        "tckAna_MaxTrack_XZ_GaussSigmaX",
        "tckAna_MaxTrack_YZ_GaussSigmaY",
        "hitsAna_z2Sigma",
        "tckAna_MaxTrackxySigmaGausBalance",
        "hitsAna_xySkew",
        "hitsAna_zSkew"};

    std::string result = "";
    for(auto name : names) {
        auto initCut = getInitialCut(dKind, name);
        auto scaling = getScalingFactor(energyIdx, name);
        // define scaled cuts
        auto sLow = initCut.first * scaling.first;
        auto sHigh = initCut.second * scaling.second;
        auto scaledDifference = abs(sHigh - sLow) * cutScaling / 2.0;
        //cout << "_________________________________________________"  << endl;
        //cout << "low cut " << sLow << "+" << scaledDifference << ", high cut " << sHigh << "-" << scaledDifference <<  " scaling = " << scaledDifference << endl;
        //cout << "scaled difference: " << scaledDifference << endl;
        if (sLow != 0) sLow = sLow + scaledDifference;
        sHigh = sHigh - scaledDifference;
        //cout << "new low cut = " << sLow << ", new high cut = " << sHigh << endl;
        //cout << "_________________________________________________"  << endl;
        if(result.size() > 0){
            result = result + "&&";
        }
        result = result + name + ">" + to_string(sLow) + "&&" +
            name + "<" + to_string(sHigh);
    }
    return result;
}

//*****************************************************************
// Function to generate the cutString regarding energy, which will be
// applied with the initial filters in the data frame.
//*****************************************************************
std::string genEnergyCut(int energyIdx){
    std::string result = "";
    auto range = getEnergyRange(energyIdx);
    result = "energy_keV > " + to_string(range.first) + " && " +
	"energy_keV < " + to_string(range.second);
    return result;
}

// *****************************************************************
// Function to define custom columns and filters for the initial DF.
// *****************************************************************
ROOT::RDF::RNode applyDefinesAndFilters(ROOT::RDataFrame df, int energyIdx, DataKind dKind){
    ROOT::RDF::RNode result = ROOT::RDataFrame(0);

    // circle radius
    //result = df.Define("radius",processName_hitsAna + "xMean*" + processName_hitsAna + "xMean+" + processName_hitsAna + "yMean*" + processName_hitsAna + "yMean");
    result = df.Define("radius","hitsAna_xMean*hitsAna_xMean+hitsAna_yMean*hitsAna_yMean");

    // Energy balance
    result = result.Define("energyBalance","(tckAna_MaxTrackEnergy_X - tckAna_MaxTrackEnergy_Y)/(tckAna_MaxTrackEnergy_X + tckAna_MaxTrackEnergy_Y)");

    // MaxTrack_XZ_GaussSigmaX
    result = result.Define("absMaxTrack_XZ_GaussSigmaX","abs(tckAna_MaxTrack_XZ_GaussSigmaX)");

    // MaxTrack_XZ_GaussSigmaY
    result = result.Define("absMaxTrack_YZ_GaussSigmaY","abs(tckAna_MaxTrack_YZ_GaussSigmaY)");

    // one track in X or two in Y and vice verce
    result = result.Define("oneXtwoY","(tckAna_nTracks_X == 1 && tckAna_nTracks_Y == 1) || (tckAna_nTracks_X == 1 && tckAna_nTracks_Y == 2) || (tckAna_nTracks_X == 2 && tckAna_nTracks_Y == 1)");

    // timestamp to date
    //result = bgData.Define("date","timeStamp.GetDateTime()");

    // veto
    //result = bgData.Define("vetoAmplitude","double(veto_MaxPeakAmplitude.second)");

    // calibrated energy
    result = calcEnergy(result, energyIdx, dKind);

    // ******************************

    // to initialize filtered dataframe
    result = result.Filter("eventID > 0");

    // tracks cut
    result = result.Filter("tckAna_nTracks_X == 1");
    result = result.Filter("tckAna_nTracks_Y == 1");

    // radius cut
    result = result.Filter("radius < 100");

    // energy cut
    //result = result.Filter("tckAna_MaxTrackEnergy > 5300 && tckAna_MaxTrackEnergy < 7500"); //5000 to 12000 or 2000 to 4000 or between
    auto cutEnergy= genEnergyCut(energyIdx);
    result = result.Filter(cutEnergy);

    // sigma cut
    result = result.Filter("tckAna_MaxTrack_XZ_GaussSigmaX > -1 && tckAna_MaxTrack_YZ_GaussSigmaY > -1"); //to not take into account the failed fits.

    return result;
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

void cut(std::string cutstring, int energyIdx, double cutScaling) {

    // ******************************
    // Important variables and declarations
    // ******************************

    string processName_tckAna = "tckAna_";                       // process prefix
    string processName_hitsAna = "hitsAna_";                         // process prefix
    string path = "/storage/cast/SR2019/analysis/xenon/official/v2.4.0/";  // data path
    const std::string commonPrefix = "trackAnalysis";
    vector<string> fileNamesCal;
    vector<string> fileNamesBg;
    vector<string> fileNamesFe55;


    // ******************************
    // List of Data
    // ******************************
    // Define your data here

    // Background runs
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
    //fileNamesBg.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.3.15/trackAnalysis_2signalThreshold/R10850*");

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

    // Fe55 runs

    fileNamesFe55.push_back(path + "trackAnalysis/R10852*");

    // X-ray tube runs
    switch(energyIdx){
        case 1:
            fileNamesCal.push_back(path + "XRayTube/R11664*");
            break;
        case 2:
            fileNamesCal.push_back(path + "XRayTube/R11665*");
            break;
        case 3:
            fileNamesCal.push_back(path + "XRayTube/R11666*");
            break;
        case 4:
            fileNamesCal.push_back(path + "XRayTube/R11667*");
            break;
        case 5:
            fileNamesCal.push_back(path + "XRayTube/R11668*");
            break;
        case 6:
            fileNamesCal.push_back(path + "XRayTube/R11659*");
            break;
        default:
            assert(false);
            break;
    }

    // ******************************
    // Create RDataFrames
    // ******************************

    ROOT::RDataFrame bgData("AnalysisTree",fileNamesBg);
    ROOT::RDataFrame calData("AnalysisTree",fileNamesCal);
    ROOT::RDataFrame Fe55Data("AnalysisTree",fileNamesFe55);

    auto bgData2 = applyDefinesAndFilters(bgData, energyIdx, dkBackground);
    auto calData2 = applyDefinesAndFilters(calData, energyIdx, dkXrayTube);
    auto Fe55Data2 = applyDefinesAndFilters(Fe55Data, energyIdx, dkBackground);


    // ************************************************************
    // Estimate run duration
    // ************************************************************

    double binw = 60;
    int rate_threshold = 10;  // counts per minute
    int run_duration = GetRunDuration(bgData, binw, rate_threshold);
    int run_duration_combined = 0;
    run_duration_combined += run_duration;

    cout << "\n\033[1;35mrun duration: " << run_duration << " s = " << (double)run_duration / 3600
         << " h = " << (double)run_duration / (3600 * 24) << " d" << endl;


    // ************************************************************
    // CUTS: generate and apply the cut strings
    // ************************************************************


    // apply cuts
    auto cutBackground = genCutString(energyIdx, dkBackground, cutScaling);
    auto cutXrayTube = genCutString(energyIdx, dkXrayTube, cutScaling);
    auto cutFe55 = genCutString(energyIdx, dkBackground, cutScaling);
    cout << "\nCuts applied to background data: " << cutBackground << endl;
    cout << "\nCuts applied to Fe55 data: " << cutFe55 << endl;
    cout << "\nCuts applied to X-ray tube data: " << cutXrayTube << endl;

    double n0_bg = bgData2.Count().GetValue();
    double n0_cal = calData2.Count().GetValue();
    double n0_fe = Fe55Data2.Count().GetValue();

    cout << "\n\n\033[1;35m# events in bg DF before: " << n0_bg  << endl;
    cout << "# events in cal DF before: " << n0_cal  << endl;
    cout << "# events in Fe55 DF before: " << n0_fe  << endl;

    bgData2 = bgData2.Filter(cutBackground);
    Fe55Data2 = Fe55Data2.Filter(cutBackground);
    calData2 = calData2.Filter(cutXrayTube);

    double n1_cal;
    double n1_fe;
    double n1_bg;


    n1_bg = bgData2.Count().GetValue();
    n1_cal = calData2.Count().GetValue();
    n1_fe = Fe55Data2.Count().GetValue();


    double FOM0_tube = n0_cal/sqrt(n0_bg);
    double FOM1_tube = n1_cal/sqrt(n1_bg);

    double FOM0_Fe55 = n0_fe/sqrt(n0_bg);
    double FOM1_Fe55 = n1_fe/sqrt(n1_bg);
/*
// Print columns' names
auto colNames = bgData_filtered.GetColumnNames();
for (auto &&colName : colNames) std::cout << colName << std::endl;
// Print column type
auto colType = bgData_filtered.GetColumnType("veto_MaxPeakAmplitude");
std::cout << "Column veto_MaxPeakAmplitude has type " << colType << std::endl;
*/

    // Draw some histograms of the filtered data

    auto histo0 = calData2.Histo1D({"tckAna_MaxTrackEnergy", "Energy X-ray tube", 100, 0, 13000}, "tckAna_MaxTrackEnergy");
    TCanvas* c00 = new TCanvas();
    histo0->DrawClone();

    auto histo1 = bgData2.Histo1D({"tckAna_MaxTrackEnergy", "Energy background", 48, 0, 13000}, "tckAna_MaxTrackEnergy");
    TCanvas* c01 = new TCanvas();
    //TCanvas c("c","x hist");
    histo1->DrawClone();

    auto histo2 = Fe55Data2.Histo1D({"tckAna_MaxTrackEnergy", "Energy Fe55 calibration", 100, 0, 13000}, "tckAna_MaxTrackEnergy");
    TCanvas* c02 = new TCanvas();
    histo2->DrawClone();


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
    cout << "\n# events in bg DF after: " <<  n1_bg << endl;
    cout << "# events in cal DF after: " <<  n1_cal << endl;
    cout << "# events in Fe55 DF after: " <<  n1_fe << endl;
    cout << "\nBackground: cut efficiency = " << 100 * n1_bg / n0_bg << " %" << endl;
    cout << "Calibration: cut efficiency = " << 100 * n1_cal / n0_cal << " %" << endl;
    cout << "Fe55: cut efficiency = " << 100 * n1_fe / n0_fe << " %" << endl;
    cout << "\nX-ray tube calibrations" << endl;
    cout << "\nFOM before cut = " << FOM0_tube << endl;
    cout << "FOM after cut = " << FOM1_tube << endl;
    cout << "FOM after cut as eff/eff_bck = " << (100 * n1_cal / n0_cal)/(sqrt(100 * n1_bg / n0_bg )) << endl;
    cout << "\nFe55 calibrations" << endl;
    cout << "\nFOM before cut = " << FOM0_Fe55 << endl;
    cout << "FOM after cut = " << FOM1_Fe55 << endl;
    cout << "FOM after cut as eff/eff_bck = " << (100 * n1_fe / n0_fe)/(sqrt(100 * n1_bg / n0_bg )) << endl;
    cout << "\n\n\033[0m";



}
