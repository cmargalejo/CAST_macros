// Author: Cristina Margalejo, December 2022
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

using namespace ROOT::VecOps;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Macro to define one or several cuts and calculate the resulting FOM
//

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

// Main function
void cut(std::string cutstring) {

    // ******************************
    // Important variables and declarations
    // ******************************

	//string processName_tckAna = "tckAna_";			     // process prefix
    //string processName_hitsAna = "hitsAna_";			     // process prefix
    string path = "/home/cristina/GitHub/iaxo-quickana/";  // data path
	//vector<string> fileNamesCal;
    vector<string> fileNamesBg;


    // ******************************
    // List of Data
    // ******************************
	// Define your data here

	// Background runs

    fileNamesBg.push_back(path + "background_df_all_runs_dataset2_subsets2and3_v2.4.0_extrapolatedCuts.root"); //v2.4.0 final ones: background_df_all_runs_dataset2_subset2_calibrated_1track_v2.4.0.root; background_df_all_runs_dataset2_subset3_calibrated_1track_v2.4.0.root background_df_all_runs_dataset2_subsets2and3_v2.4.0_extrapolatedCuts.root


	// ******************************
	// Create RDataFrames
	// ******************************

	ROOT::RDataFrame  bgData("AnalysisTree",fileNamesBg);
	//ROOT::RDataFrame calData("AnalysisTree",fileNamesCal);

	ROOT::RDF::RNode bgData2 = ROOT::RDataFrame(0);
	//ROOT::RDF::RNode calData2 = ROOT::RDataFrame(0);

	// ******************************

	// ******************************
	// Estimate run duration
	// ******************************

	double binw = 60;
	int rate_threshold = 10;  // counts per minute
	int run_duration = GetRunDuration(bgData, binw, rate_threshold);
    int run_duration_combined = 0;
	run_duration_combined += run_duration;

	//cout << "\n\033[1;31mRun Nr. " << runNumber[f] << ":\033[0m run duration: " << run_duration << " s = " << (double)run_duration / 3600
	//     << " h = " << (double)run_duration / (3600 * 24) << " d" << endl;
	cout << "run duration: " << run_duration << " s = " << (double)run_duration / 3600
	     << " h = " << (double)run_duration / (3600 * 24) << " d" << endl;



	// Define custom columns

	// circle radius
	bgData2 = bgData.Define("radius","hitsAna_xMean*hitsAna_xMean+hitsAna_yMean*hitsAna_yMean");
	//calData2 = calData.Define("radius",hitsAna_xMean*" hitsAna_xMean+" hitsAna_yMean*" hitsAna_yMean");

    //Errors
    //bgData2 = bgData.Define("error","sqrt(energy_keV)");
	// Energy balance
	/*
	bgData2 = bgData2.Define("energyBalance","(tckAna_MaxTrackEnergy_X -tckAna_MaxTrackEnergy_Y)/(tckAna_MaxTrackEnergy_X+tckAna_MaxTrackEnergy_Y)");
	//calData2 = calData2.Define("energyBalance","(" tckAna_MaxTrackEnergy_X -"tckAna_MaxTrackEnergy_Y)/("+processName_tckAna + "MaxTrackEnergy_X+" tckAna_MaxTrackEnergy_Y)");

	// MaxTrack_XZ_GaussSigmaX
	bgData2 = bgData2.Define("absMaxTrack_XZ_GaussSigmaX","abs(tckAna_MaxTrack_XZ_GaussSigmaX)");
	//calData2 = calData2.Define("absMaxTrack_XZ_GaussSigmaX","abs(" tckAna_MaxTrack_XZ_GaussSigmaX)");

	// MaxTrack_XZ_GaussSigmaY
	bgData2 = bgData2.Define("absMaxTrack_YZ_GaussSigmaY","abs(tckAna_MaxTrack_YZ_GaussSigmaY)");
	//calData2 = calData2.Define("absMaxTrack_YZ_GaussSigmaY","abs(" tckAna_MaxTrack_YZ_GaussSigmaY)");

	// one track in X or two in Y and vice verce
	bgData2 = bgData2.Define("oneXtwoY","(tckAna_nTracks_X == 1 && tckAna_nTracks_Y == 1) || (tckAna_nTracks_X == 1 && tckAna_nTracks_Y == 2) || (tckAna_nTracks_X == 2 && tckAna_nTracks_Y == 1)");
	//calData2 = calData2.Define("oneXtwoY","(tckAna_nTracks_X == 1 && tckAna_nTracks_Y == 1) || (tckAna_nTracks_X == 1 && tckAna_nTracks_Y == 2) || (tckAna_nTracks_X == 2 && tckAna_nTracks_Y == 1)");
*/
	// timestamp to date
	//bgData2 = bgData2.Define("date","timeStamp.GetDateTime()");
	//calData2 = calData2.Define("date","timeStamp.GetDate()");

	// veto
	auto extractSecond = [](std::map<int, double> vec){
	    Double_t result = 0;
	    if(vec.size() != 1){
	        std::cout << "not what expected. size is " << vec.size() << std::endl;
	        //abort();
	    }
	    for(auto pair : vec){
	        result = pair.second;
	    }
	    return result;
	};
	bgData2 = bgData2.Define("vetoAmplitude", extractSecond, {"veto_MaxPeakAmplitude"});//"vetoAmplitude","veto_MaxPeakAmplitude.second");
    bgData2 = bgData2.Define("vetoTime", extractSecond, {"veto_PeakTime"});//"vetoAmplitude","veto_MaxPeakAmplitude.second");
	//bgData2 = bgData2.Define("vetoTime","veto_PeakTime.second");
	//calData2 = calData.Define("vetoAmplitude","timeStamp.GetDate()");

	// calibrated energy
	//bgData2 = bgData2.Define("energy_keV","tckAna_MaxTrackEnergy/1070");
	//calData2 = calData2.Define("energy_keV","tckAna_MaxTrackEnergy/1070");

	// ******************************

	// initialize filtered dataframe
	auto bgData_filtered = bgData2.Filter("eventID > 0");
	//auto calData_filtered = calData2.Filter("eventID > 0");

	// tracks cut
	//bgData_filtered = bgData_filtered.Filter("tckAna_nTracks_X == 1 && tckAna_nTracks_Y == 1");
	//calData_filtered = calData_filtered.Filter("tckAna_nTracks_X == 1 && tckAna_nTracks_Y == 1");

	// radius cut
	//bgData_filtered = bgData_filtered.Filter("radius < 100 ");
	//calData_filtered = calData_filtered.Filter("radius < 100");

	//ring cuts
	//bgData_filtered = bgData_filtered.Filter("hitsAna_xMean*hitsAna_xMean+hitsAna_yMean*hitsAna_yMean > 0 && hitsAna_xMean*hitsAna_xMean+hitsAna_yMean*hitsAna_yMean < 25"); //0 to 5 mm
	//bgData_filtered = bgData_filtered.Filter("hitsAna_xMean*hitsAna_xMean+hitsAna_yMean*hitsAna_yMean >= 25 && hitsAna_xMean*hitsAna_xMean+hitsAna_yMean*hitsAna_yMean < 100"); //5 to 10 mm
	//bgData_filtered = bgData_filtered.Filter("hitsAna_xMean*hitsAna_xMean+hitsAna_yMean*hitsAna_yMean >= 100 && hitsAna_xMean*hitsAna_xMean+hitsAna_yMean*hitsAna_yMean < 225"); //10 to 15 mm
	//bgData_filtered = bgData_filtered.Filter("hitsAna_xMean*hitsAna_xMean+hitsAna_yMean*hitsAna_yMean >= 225 && hitsAna_xMean*hitsAna_xMean+hitsAna_yMean*hitsAna_yMean < 400"); //15 to 20 mm
	//bgData_filtered = bgData_filtered.Filter("hitsAna_xMean*hitsAna_xMean+hitsAna_yMean*hitsAna_yMean >= 400 && hitsAna_xMean*hitsAna_xMean+hitsAna_yMean*hitsAna_yMean < 625"); //15 to 25 mm
	//bgData_filtered = bgData_filtered.Filter("hitsAna_xMean*hitsAna_xMean+hitsAna_yMean*hitsAna_yMean >= 625 && hitsAna_xMean*hitsAna_xMean+hitsAna_yMean*hitsAna_yMean < 900"); //25 to 30 mm


	//veto cut
	//bgData_filtered = bgData_filtered.Filter("vetoAmplitude < 200 || (vetoAmplitude >= 200 && (vetoTime >300 || vetoTime < 100))");
	//bgData_filtered = bgData_filtered.Filter("veto_MaxPeakAmplitude.second < 200");

	// energy cut
	//bgData_filtered = bgData_filtered.Filter("tckAna_MaxTrackEnergy > 5617 && tckAna_MaxTrackEnergy < 7490"); //5000 to 12000 or 2000 to 4000 or between, 5000 to 7500
	//calData_filtered = calData_filtered.Filter("tckAna_MaxTrackEnergy > 5617 && tckAna_MaxTrackEnergy < 7490");
	//bgData_filtered = bgData_filtered.Filter("energy_keV > 7.5 && energy_keV < 10");

	// sigma cut
	//bgData_filtered = bgData_filtered.Filter("tckAna_MaxTrack_XZ_GaussSigmaX > -1 && tckAna_MaxTrack_YZ_GaussSigmaY > -1"); //to not take into account the failed fits.
	//calData_filtered = calData_filtered.Filter("tckAna_MaxTrack_XZ_GaussSigmaX > -1 && tckAna_MaxTrack_YZ_GaussSigmaY > -1");

	double n0_bg = bgData_filtered.Count().GetValue();
	//double n0_cal = calData_filtered.Count().GetValue();

	//double n1_cal;
	double n1_bg;

	// ******************************
	// CUTS
	// ******************************

	cout <<"\n\n\033[1;35m***************************************************************************" << endl;
    cout << "           Cut: " << cutstring << endl;
    cout << "***************************************************************************\n\n\033[0m" << endl;

	cout << "\n\n\033[1;35m# events in bg DF before: " << n0_bg  << endl;
	//cout << "\n\n\033[1;35m# events in cal DF before: " << n0_cal  << endl;

	// apply cuts
	bgData_filtered = bgData_filtered.Filter(cutstring);
	//calData_filtered = calData_filtered.Filter(cutstring);

	n1_bg = bgData_filtered.Count().GetValue();
	//n1_cal = calData_filtered.Count().GetValue();

	//double FOM0 = n0_cal/sqrt(n0_bg);
	//double FOM1 = n1_cal/sqrt(n1_bg);

    // Print columns' names
	//auto colNames = bgData_filtered.GetColumnNames();
    //for (auto &&colName : colNames) std::cout << colName << std::endl;

    // Print column type
    /*
    auto colType = bgData_filtered.GetColumnType("veto_MaxPeakAmplitude.second");
    std::cout << "Column veto_MaxPeakAmplitude has type " << colType << std::endl;
    auto colType2 = bgData_filtered.GetColumnType("veto_PeakTime");
    std::cout << "Column veto_PeakTime has type " << colType2 << std::endl;
    auto colType3 = bgData_filtered.GetColumnType("vetoAmplitude");
    std::cout << "Column vetoAmplitude has type " << colType3 << std::endl;
    auto colType4 = bgData_filtered.GetColumnType("energy_keV");
    std::cout << "Column energy_keV has type " << colType4 << std::endl;
    */

    // Draw some histograms of the filtered data

    auto histo7 = bgData_filtered.Histo1D({"energy_keV", "Energy", 48, 0, 12}, "energy_keV");
    TCanvas* c07 = new TCanvas();
    histo7->DrawClone();

    TCanvas c = new TCanvas("c", "canvas", 1200, 400);
    //c.SetCanvasSize(1200,800);
    c.Divide(2,2);
    /*

    auto histo1 = bgData_filtered.Histo1D({"tckAna_MaxTrackEnergy", "Energy", 48, 0, 13000}, "tckAna_MaxTrackEnergy");
    TCanvas* c01 = new TCanvas();
    histo1->DrawClone();
    */
    auto histo2 = bgData_filtered.Histo1D({"energy_keV", "Energy", 48, 0, 12}, "energy_keV");
    //TCanvas* c02 = new TCanvas();
    c.cd(1);
    histo2->DrawClone();
    /*
    auto histo3 = bgData_filtered.Histo1D({"timeStamp", "timeStamp", 10000, 0, 17e9}, "timeStamp");
    TCanvas* c03 = new TCanvas();
    histo3->DrawClone();
    */

    auto histo4 = bgData_filtered.Histo1D("veto_MaxPeakAmplitude.second");
    //TCanvas* c04 = new TCanvas();
    //TCanvas c("c","x hist");
    c.cd(2);
    histo4->DrawClone();

    auto histo5 = bgData_filtered.Histo2D({"hitsMap", "hitsMap", 120,-30,30,120,-30,30},"tckAna_MaxTrack_Xmean_X","tckAna_MaxTrack_Ymean_Y");
    //TCanvas* c05 = new TCanvas();
    //TCanvas c("c","x hist");
    c.cd(3);
    histo5->DrawClone("colz");

    auto histo6 = bgData_filtered.Histo1D({"energyBalance", "energyBalance", 300, -1, 1}, "tckAna_MaxTrackEnergyBalanceXY");
    TCanvas* c06 = new TCanvas();
    histo6->DrawClone();


    c.SaveAs("histograms.pdf");


/*
	// Preparing the RResultPtr<RDisplay> object with all columns and default number of entries
    //auto d1 = bgData_filtered.Display("",128);
    // Preparing the RResultPtr<RDisplay> object with two columns and 128 entries
    auto d2 = bgData_filtered.Display({"timeStamp", "eventID", "energy_keV", "vetoTime", "vetoAmplitude"}, 1280);
    // Printing the short representations, the event loop will run
    //d1->Print();
    d2->Print();
*/
	cout << "# events in bg DF after: " <<  n1_bg << endl;
	//cout << "# events in cal DF after: " <<  n1_cal << endl;
	cout << "\nBackground: cut efficiency = " << 100 * n1_bg / n0_bg << " %" << endl;
	//cout << "Calibration: cut efficiency = " << 100 * n1_cal / n0_cal << " %" << endl;
	//cout << "\n FOM before cut = " << FOM0 << endl;
	//cout << "FOM after cut = " << FOM1 << endl;
	cout << "\n\n\033[0m";



}

