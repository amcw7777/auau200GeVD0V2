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
    ~D0CorV2Plotter() {mLog.close();}

    void Init(TString, bool, bool);
    vector<vector<double> > getD0V2(vector<double> &hadronV2, TH1D *candOverSignal, int indexBkg);
    vector<vector<double> > getHadronV2();
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
    pair<double,double> fit_hist(TH1D *histo, TCanvas *cfg, int iptbin ,double nSigma,double fitArray[3]);
    TH1D *getSignalV2(TH1D *candOverSgn,TH1D *bkgV2,TH1D *candV2);
    double getV2Error(double matr[]);
};

inline double D0CorV2Plotter::getV2Error(double matr[6])
{
  return sqrt( pow((matr[2]-matr[1])*matr[3],2) + pow(matr[0]*matr[5],2) + pow((1-matr[0])*matr[4],2) );
}
