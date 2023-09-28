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

    //Ar dataset 1
    /*
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
   */
   //Ar dataset 2
   /*
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10671*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10690*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10691*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10781*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10793*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10795*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10796*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10800*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10802*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10805*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10807*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10811*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10822*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10822*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10824*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10826*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10828*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10830*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10832*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10834*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10836*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10838*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10840*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10842*");
    fileNames.push_back("/storage/cast/SR2019/analysis/argon/official/v2.4.0/final/R10843*");
    */
    // Xe dataset
    //Bottle 1, filter 1
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10710*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10726*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10730*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10740*");

    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10746*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10750*");

    //Hit channels
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10772*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10777*");

    //Back to all channels
    // Bottle 2, filter 1
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10850*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10853*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10857*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10861*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10863*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10865*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10867*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10869*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10871*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10873*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10875*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10877*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10879*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10895*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10897*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10900*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10902*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10904*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10906*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10908*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10910*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10912*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10914*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10918*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10920*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10922*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10924*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10927*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10932*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10934*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10936*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10938*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10941*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10943*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10945*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10947*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10949*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10951*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10953*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10955*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10957*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10959*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10961*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10963*");

    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10965*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10967*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10969*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10971*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10973*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10975*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10977*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10979*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10981*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10984*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10987*");
    // Bottle 2, filter 2
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R10998*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11003*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11012*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11014*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11018*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11021*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11023*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11025*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11027*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11029*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11031*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11033*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11035*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11038*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11040*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11042*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11044*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11046*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11048*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11050*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11052*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11054*");
    // Bottle 2, filter 3, 2 ln/h
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11073*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11078*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11081*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11083*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11086*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11089*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11091*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11093*");
    // Bottle 2, filter 3, 5 ln/h
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11108*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11110*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11113*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11115*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11117*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11120*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11132*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11134*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11149*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11153*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11205*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11208*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11210*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11218*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11223*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11225*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11227*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11230*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11241*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11245*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11250*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11252*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11254*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11281*");
    // subrun 4 is not readable! Why?
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11284_0000[0123]*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11286*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11289*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11293*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11296*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11298*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11300*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11304*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11306*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11308*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11310*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11313*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11315*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11317*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11319*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11321*");
    //Bottle 2, filter 3, small volume
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11339*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11341*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11343*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11344*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11347*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11353*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11358*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11361*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11367*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11375*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11378*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11381*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11384*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11393*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11398*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11402*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11404*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11406*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11408*");
    // Bottle 2, filter 3, large volume
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11411*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11419*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11422*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11424*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11426*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11430*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11433*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11435*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11437*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11439*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11456*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11459*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11461*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11463*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11467*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11470*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11472*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11476*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11478*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11510*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11517*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11519*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11521*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11523*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11528*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11533*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11535*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11538*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11540*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11542*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11544*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11546*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11549*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11552*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11554*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11556*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11560*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11563*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11565*");
    fileNames.push_back("/storage/cast/SR2019/analysis/xenon/official/v2.4.0/final/SC_new/R11567*");




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
        /*");
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

