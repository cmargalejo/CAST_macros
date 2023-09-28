#include <TH1D.h>
#include <vector>
#include <string>
#include <iostream>
#include <ROOT/RDataFrame.hxx>

using namespace std;

int GetRunDuration(ROOT::RDF::RNode data_filtered, double binw = 60, int rate_threshold = 10) {
    int bincount;
    int run_duration = 0;

    auto result_ts = data_filtered.Take<double>("timeStamp");
    vector<double> timestamps = result_ts.GetValue();

    double tmax = data_filtered.Max("timeStamp").GetValue();
    double tmin = data_filtered.Min("timeStamp").GetValue();

    tmax = tmin + (round((tmax - tmin + binw / 2) / binw) * binw);

    TH1D* ht = new TH1D("timestamps", "timestamps", (tmax - tmin) / binw, tmin, tmax);
    for (auto i : timestamps) {
        ht->Fill(i);
    }

    int nbins = ht->GetNbinsX();
    for (unsigned int i = 0; i < nbins; i++) {
        bincount = ht->GetBinContent(i);
        if (bincount > rate_threshold) {
            run_duration += binw;
        }
    }

    delete ht;
    return run_duration;
}

void runDuration() {
    vector<string> fileNames;
    // Add your filenames to this vector

    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10178*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10182*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10188*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10190*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10192*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10194*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10196*");

    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10201*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10204*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10206*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10212*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10214*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10216*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10218*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10220*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10222*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10224*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10228*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10234*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10236*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10238*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10240*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10242*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10244*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10246*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10248*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10250*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10252*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10254*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10256*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10258*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10260*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10264*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10266*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10268*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10270*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10272*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10274*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10276*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10278*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10283*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10285*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10287*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10289*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10291*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10293*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10295*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10297*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10307*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10309*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10313*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10323*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10325*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10327*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10329*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10331*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10333*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10335*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10337*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10339*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10341*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10342*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10344*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10346*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10350*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10352*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10355*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10357*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10359*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10361*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10363*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10370*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10372*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10374*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10376*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10382*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10384*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10386*"); //Christmas break, no SC data
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10388*"); //Christmas break, no SC data
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10390*"); //Christmas break, no SC data
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10402*"); //Christmas break, no SC data
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10404*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10413*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10415*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10418*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10420*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10422*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10424*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10426*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10428*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10430*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10432*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10434*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10436*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10438*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10444*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10446*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10448*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10450*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10454*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10468*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10470*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10472*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10474*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10476*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10478*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10480*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10482*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10484*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10486*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10488*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10490*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10492*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10494*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10496*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10498*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10500*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10504*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10506*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10508*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10510*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10512*");
    //from dataset 2
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10648*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10660*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10663*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/SC_new/R10665*");


    double binw = 60;
    int rate_threshold = 10;

    int total_duration = 0;
    int total_duration_0 = 0;
    int total_duration_positive = 0;

    for (const auto& file : fileNames) {
        ROOT::RDataFrame df("AnalysisTree", file);

        // Entire DataFrame
        int duration = GetRunDuration(df, binw, rate_threshold);
        total_duration += duration;

        // scTck_SolarTracking == 0
        auto df_SolarTracking0 = df.Filter("scTck_SolarTracking == 0");
        int duration_0 = GetRunDuration(df_SolarTracking0, binw, rate_threshold);
        total_duration_0 += duration_0;

        // scTck_SolarTracking > 0
        auto df_SolarTrackingPositive = df.Filter("scTck_SolarTracking > 0");
        int duration_positive = GetRunDuration(df_SolarTrackingPositive, binw, rate_threshold);
        total_duration_positive += duration_positive;
        /*
        cout << "\n\033[1;35mFile: " << file << endl;
        cout << "Run duration: " << duration << " s = " << static_cast<double>(duration) / 3600 << " h = " << static_cast<double>(duration) / (3600 * 24) << " d" << endl;
        cout << "Run duration (scTck_SolarTracking == 0): " << duration_0 << " s = " << static_cast<double>(duration_0) / 3600 << " h = " << static_cast<double>(duration_0) / (3600 * 24) << " d" << endl;
        cout << "Run duration (scTck_SolarTracking > 0): " << duration_positive << " s = " << static_cast<double>(duration_positive) / 3600 << " h = " << static_cast<double>(duration_positive) / (3600 * 24) << " d" << endl;
        */
    }

    cout << "==================================================================" << endl;
    cout << "\n\033[1;35mTotal run duration: " << total_duration << " s = " << static_cast<double>(total_duration) / 3600 << " h = " << static_cast<double>(total_duration) / (3600 * 24) << " d" << endl;
    cout << "Total run duration (scTck_SolarTracking == 0): " << total_duration_0 << " s = " << static_cast<double>(total_duration_0) / 3600 << " h = " << static_cast<double>(total_duration_0) / (3600 * 24) << " d" << endl;
    cout << "Total run duration (scTck_SolarTracking > 0): " << total_duration_positive << " s = " << static_cast<double>(total_duration_positive) / 3600 << " h = " << static_cast<double>(total_duration_positive) / (3600 * 24) << " d" << endl;
}

