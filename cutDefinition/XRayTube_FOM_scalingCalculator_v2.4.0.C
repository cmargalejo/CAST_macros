

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
        case 1: return 1356;
        case 2: return 1472;
        case 3: return 1512;
        case 4: return 1478;
        case 5: return 1458;
        case 6: return 1446;
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
            return df.Define("energy_keV","tckAna_MaxTrackEnergy/1170");
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
        { "hitsAna_balanceXYenergy", {
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
        { "hitsAna_z2Sigma", {
                { 1,  { 0.0, 10.5 }  }, //7 before revising the cuts
                { 2,  { 0.0, 7.0 }  }, //5
                { 3,  { 0.0, 5.8 }  }, //5
                { 4,  { 0.0, 3.8 }  },
                { 5,  { 0.0, 1.0 }  },
                { 6,  { 0.0, 1.6 }  },
          }
        },
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
std::pair<double, double> getInitialCut(DataKind dKind, std::string obs){
    // Given the data kind and observable, return the initial cut.
    std::map<std::string, std::pair<double, double>> castMap = {
        { "hitsAna_balanceXYenergy",  { -0.25, 0.25 }  },
        { "tckAna_MaxTrack_XZ_GaussSigmaX",  { 0.15, 0.75 }  },
        { "tckAna_MaxTrack_YZ_GaussSigmaY",  { 0.15, 0.75 }  },
        { "hitsAna_z2Sigma",  { 0.0, 0.95 }  },
        { "tckAna_MaxTrackxySigmaGausBalance",  { -0.55, 0.55 }  },
        { "hitsAna_xySkew",  { -0.85, 0.85 }  },
        { "hitsAna_zSkew",  { -2.5, 2.5 }  },
    };
    /*
    std::map<std::string, std::pair<double, double>> xrayTubeMap = {
        { "hitsAna_balanceXYenergy",  { -0.165, 0.174 }  },
        { "tckAna_MaxTrack_XZ_GaussSigmaX",  { 0.225, 0.775 }  },
        { "tckAna_MaxTrack_YZ_GaussSigmaY",  { 0.225, 0.875 }  },
        { "hitsAna_z2Sigma",  { 0.0, 1.25 }  },
        { "tckAna_MaxTrackxySigmaGausBalance",  { -0.7, 0.7 }  },
        { "hitsAna_xySkew",  { -0.52, 0.6 }  },
        { "hitsAna_zSkew",  { -2.52, 2.84 }  },
    };
    */

    std::map<std::string, std::pair<double, double>> xrayTubeMap = {
        { "hitsAna_balanceXYenergy",  { -0.21, 0.21 }  },
        { "tckAna_MaxTrack_XZ_GaussSigmaX",  { 0.225, 0.775 }  },
        { "tckAna_MaxTrack_YZ_GaussSigmaY",  { 0.225, 0.825 }  },
        { "hitsAna_z2Sigma",  { 0.0, 1.25 }  },
        { "tckAna_MaxTrackxySigmaGausBalance",  { -0.5, 0.5 }  },
        { "hitsAna_xySkew",  { -0.52, 0.52 }  },
        { "hitsAna_zSkew",  { -3.5, 3.5 }  },
    };

    /*
    std::map<std::string, std::pair<double, double>> xrayTubeMap = {
        { "hitsAna_balanceXYenergy",  { -0.21, 0.21 }  },
        { "tckAna_MaxTrack_XZ_GaussSigmaX",  { 0.225, 0.825 }  },
        { "tckAna_MaxTrack_YZ_GaussSigmaY",  { 0.225, 0.825 }  },
        { "hitsAna_z2Sigma",  { 0.0, 0.95 }  },
        { "tckAna_MaxTrackxySigmaGausBalance",  { -0.55, 0.55 }  },
        { "hitsAna_xySkew",  { -0.52, 0.52 }  },
        { "hitsAna_zSkew",  { -2.5, 2.5 }  },
    };
    */

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
        "hitsAna_balanceXYenergy",
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
    string path = "/storage/cast/SR2019/analysis/argon/official/v2.4.0/";  // data path
    vector<string> fileNamesCal;
    vector<string> fileNamesBg;
    vector<string> fileNamesFe55;


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

    fileNamesBg.push_back(path + "trackAnalysis/R10178*");
    fileNamesBg.push_back(path + "trackAnalysis/R10182*");
    fileNamesBg.push_back(path + "trackAnalysis/R10188*");
    fileNamesBg.push_back(path + "trackAnalysis/R10190*");
    fileNamesBg.push_back(path + "trackAnalysis/R10192*");
    fileNamesBg.push_back(path + "trackAnalysis/R10194*");
    fileNamesBg.push_back(path + "trackAnalysis/R10196*");

    fileNamesBg.push_back(path + "trackAnalysis/R10201*");
    fileNamesBg.push_back(path + "trackAnalysis/R10204*");
    fileNamesBg.push_back(path + "trackAnalysis/R10206*");
    fileNamesBg.push_back(path + "trackAnalysis/R10212*");
    fileNamesBg.push_back(path + "trackAnalysis/R10214*");
    fileNamesBg.push_back(path + "trackAnalysis/R10216*");
    fileNamesBg.push_back(path + "trackAnalysis/R10218*");
    fileNamesBg.push_back(path + "trackAnalysis/R10220*");
    fileNamesBg.push_back(path + "trackAnalysis/R10222*");
    fileNamesBg.push_back(path + "trackAnalysis/R10224*");
    fileNamesBg.push_back(path + "trackAnalysis/R10228*");
    fileNamesBg.push_back(path + "trackAnalysis/R10234*");
    fileNamesBg.push_back(path + "trackAnalysis/R10236*");
    fileNamesBg.push_back(path + "trackAnalysis/R10238*");
    fileNamesBg.push_back(path + "trackAnalysis/R10240*");
    fileNamesBg.push_back(path + "trackAnalysis/R10242*");
    fileNamesBg.push_back(path + "trackAnalysis/R10244*");
    fileNamesBg.push_back(path + "trackAnalysis/R10246*");
    fileNamesBg.push_back(path + "trackAnalysis/R10248*");
    fileNamesBg.push_back(path + "trackAnalysis/R10250*");
    fileNamesBg.push_back(path + "trackAnalysis/R10252*");
    fileNamesBg.push_back(path + "trackAnalysis/R10254*");
    fileNamesBg.push_back(path + "trackAnalysis/R10256*");
    fileNamesBg.push_back(path + "trackAnalysis/R10258*");
    fileNamesBg.push_back(path + "trackAnalysis/R10260*");
    fileNamesBg.push_back(path + "trackAnalysis/R10264*");
    fileNamesBg.push_back(path + "trackAnalysis/R10266*");
    fileNamesBg.push_back(path + "trackAnalysis/R10268*");
    fileNamesBg.push_back(path + "trackAnalysis/R10270*");
    fileNamesBg.push_back(path + "trackAnalysis/R10272*");
    fileNamesBg.push_back(path + "trackAnalysis/R10274*");
    fileNamesBg.push_back(path + "trackAnalysis/R10276*");
    fileNamesBg.push_back(path + "trackAnalysis/R10278*");
    fileNamesBg.push_back(path + "trackAnalysis/R10283*");
    fileNamesBg.push_back(path + "trackAnalysis/R10285*");
    fileNamesBg.push_back(path + "trackAnalysis/R10287*");
    fileNamesBg.push_back(path + "trackAnalysis/R10289*");
    fileNamesBg.push_back(path + "trackAnalysis/R10291*");
    fileNamesBg.push_back(path + "trackAnalysis/R10293*");
    fileNamesBg.push_back(path + "trackAnalysis/R10295*");
    fileNamesBg.push_back(path + "trackAnalysis/R10297*");
    fileNamesBg.push_back(path + "trackAnalysis/R10307*");
    fileNamesBg.push_back(path + "trackAnalysis/R10309*");
    fileNamesBg.push_back(path + "trackAnalysis/R10313*");
    fileNamesBg.push_back(path + "trackAnalysis/R10323*");
    fileNamesBg.push_back(path + "trackAnalysis/R10325*");
    fileNamesBg.push_back(path + "trackAnalysis/R10327*");
    fileNamesBg.push_back(path + "trackAnalysis/R10329*");
    fileNamesBg.push_back(path + "trackAnalysis/R10331*");
    fileNamesBg.push_back(path + "trackAnalysis/R10333*");
    fileNamesBg.push_back(path + "trackAnalysis/R10335*");
    fileNamesBg.push_back(path + "trackAnalysis/R10337*");
    fileNamesBg.push_back(path + "trackAnalysis/R10339*");
    fileNamesBg.push_back(path + "trackAnalysis/R10341*");
    fileNamesBg.push_back(path + "trackAnalysis/R10342*");
    fileNamesBg.push_back(path + "trackAnalysis/R10344*");
    fileNamesBg.push_back(path + "trackAnalysis/R10346*");
    fileNamesBg.push_back(path + "trackAnalysis/R10350*");
    fileNamesBg.push_back(path + "trackAnalysis/R10352*");
    fileNamesBg.push_back(path + "trackAnalysis/R10355*");
    fileNamesBg.push_back(path + "trackAnalysis/R10357*");
    fileNamesBg.push_back(path + "trackAnalysis/R10359*");
    fileNamesBg.push_back(path + "trackAnalysis/R10361*");
    fileNamesBg.push_back(path + "trackAnalysis/R10363*");
    fileNamesBg.push_back(path + "trackAnalysis/R10370*");
    fileNamesBg.push_back(path + "trackAnalysis/R10372*");
    fileNamesBg.push_back(path + "trackAnalysis/R10374*");
    fileNamesBg.push_back(path + "trackAnalysis/R10376*");
    fileNamesBg.push_back(path + "trackAnalysis/R10382*");
    fileNamesBg.push_back(path + "trackAnalysis/R10384*");
    fileNamesBg.push_back(path + "trackAnalysis/R10386*");
    fileNamesBg.push_back(path + "trackAnalysis/R10388*");
    fileNamesBg.push_back(path + "trackAnalysis/R10390*");
    fileNamesBg.push_back(path + "trackAnalysis/R10402*");
    fileNamesBg.push_back(path + "trackAnalysis/R10404*");
    fileNamesBg.push_back(path + "trackAnalysis/R10413*");
    fileNamesBg.push_back(path + "trackAnalysis/R10415*");
    fileNamesBg.push_back(path + "trackAnalysis/R10418*");
    fileNamesBg.push_back(path + "trackAnalysis/R10420*");
    fileNamesBg.push_back(path + "trackAnalysis/R10422*");
    fileNamesBg.push_back(path + "trackAnalysis/R10424*");
    fileNamesBg.push_back(path + "trackAnalysis/R10426*");
    fileNamesBg.push_back(path + "trackAnalysis/R10428*");
    fileNamesBg.push_back(path + "trackAnalysis/R10430*");
    fileNamesBg.push_back(path + "trackAnalysis/R10432*");
    fileNamesBg.push_back(path + "trackAnalysis/R10434*");
    fileNamesBg.push_back(path + "trackAnalysis/R10436*");
    fileNamesBg.push_back(path + "trackAnalysis/R10438*");
    fileNamesBg.push_back(path + "trackAnalysis/R10444*");
    fileNamesBg.push_back(path + "trackAnalysis/R10446*");
    fileNamesBg.push_back(path + "trackAnalysis/R10448*");
    fileNamesBg.push_back(path + "trackAnalysis/R10450*");
    fileNamesBg.push_back(path + "trackAnalysis/R10454*");
    fileNamesBg.push_back(path + "trackAnalysis/R10468*");
    fileNamesBg.push_back(path + "trackAnalysis/R10470*");
    fileNamesBg.push_back(path + "trackAnalysis/R10472*");
    fileNamesBg.push_back(path + "trackAnalysis/R10474*");
    fileNamesBg.push_back(path + "trackAnalysis/R10476*");
    fileNamesBg.push_back(path + "trackAnalysis/R10478*");
    fileNamesBg.push_back(path + "trackAnalysis/R10480*");
    fileNamesBg.push_back(path + "trackAnalysis/R10482*");
    fileNamesBg.push_back(path + "trackAnalysis/R10484*");
    fileNamesBg.push_back(path + "trackAnalysis/R10486*");
    fileNamesBg.push_back(path + "trackAnalysis/R10488*");
    fileNamesBg.push_back(path + "trackAnalysis/R10490*");
    fileNamesBg.push_back(path + "trackAnalysis/R10492*");
    fileNamesBg.push_back(path + "trackAnalysis/R10494*");
    fileNamesBg.push_back(path + "trackAnalysis/R10496*");
    fileNamesBg.push_back(path + "trackAnalysis/R10498*");
    fileNamesBg.push_back(path + "trackAnalysis/R10500*");
    fileNamesBg.push_back(path + "trackAnalysis/R10504*");
    fileNamesBg.push_back(path + "trackAnalysis/R10506*");
    fileNamesBg.push_back(path + "trackAnalysis/R10508*");
    fileNamesBg.push_back(path + "trackAnalysis/R10510*");
    fileNamesBg.push_back(path + "trackAnalysis/R10512*");


    // Fe55 runs

    fileNamesFe55.push_back(path + "trackAnalysis/R10513*");

    // X-ray tube runs
    switch(energyIdx){
        case 1:
            fileNamesCal.push_back(path + "XRayTube/R00444*"); //1700 to 2400 ADC
            break;
        case 2:
            fileNamesCal.push_back(path + "XRayTube/R00446*"); //2600 to 3550
            break;
        case 3:
            fileNamesCal.push_back(path + "XRayTube/R00447*"); //3950 to 5020
            break;
        case 4:
            fileNamesCal.push_back(path + "XRayTube/R00448*"); //5900 to 7200
            break;
        case 5:
            fileNamesCal.push_back(path + "XRayTube/R00449*"); //7800 to 9200
            break;
        case 6:
            fileNamesCal.push_back(path + "XRayTube/R00450*"); //10500 to 12500
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
