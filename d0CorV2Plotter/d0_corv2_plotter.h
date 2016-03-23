/***********************************************************************************
 *
 * d0_corv2_plotter
 *
 * Author: Leon He
 ***********************************************************************************
 *
 * Description: 
 *
 ***********************************************************************************
 *
 * Log:
 *
 ***********************************************************************************/

class TString;
class ofstream;
class TH1D;
class TCanvas;
class TFile;
class vector;

using namespace std;

class D0CorV2Plotter
{
  public:
   D0CorV2Plotter() {cout<<"D0 Correlation v2 plotter is running"<<endl;}
    ~D0CorV2Plotter() {}

    void Init(TString, bool, bool);
    vector<vector<double> > getD0V2(vector<double> &hadronV2, TH1D *candOverSignal, int indexBkg);
    vector<vector<double> > getHadronV2();
    pair<double,double> fit_hist(TH1D *histo, TCanvas *cfg, int iptbin ,double nSigma,double fitArray[3]);
    vector<vector<double> > fitMass();

  private:
    TString mInputFileName;
    int mIndexEtaGap;
    TString mFlattenNames[5];
    TFile *mOutputFile;
    TFile *mInputFile;
    ofstream mLog;
    bool mIsYieldSystematic;
    bool mIsSystematic;
};
