#include "d0_corv2_plotter.cc"
#include "algorithm"

using namespace std;
// vector<vector<double> > getD0V2(bool isSystematic,bool isYieldSystematic,TString inputFileName);
// vector<double> getYieldSys();
// vector<double> getEfficiencySys();
// vector<double> getBackgroundSys();

void runD0CorV2Plotter(TString tempFileName="test.root")
{
  getD0V2(false,false,tempFileName.Data());//fill in centra value;
 
}
void runD0CorV2Plotter()
{
  // getD0V2(false,false,"auau2014D0CorV2.root");//fill in centra value;

  vector<double> yieldSys = getYieldSys();
  vector<double> effSys= getEfficiencySys();
  vector<double> bkgSys = getBackgroundSys();
  ofstream outputFile("toM.txt",ofstream::app);
  outputFile<<"====================Systematics Uncertainty====================="<<endl;
  outputFile<<"pT:\t\t\t\tYieldSys\t\t\t\tEff. Sys.\t\t\t\tBkg. Sys.\t\t\t\t"<<endl;
  for(int i=0;i<6;i++)
    outputFile<<i<<"\t\t\t\t"<<yieldSys[i]<<"\t\t\t\t"<<effSys[i]<<"\t\t\t\t"<<bkgSys[i]<<endl;
}

vector<vector<double> > getD0V2(bool inputIsSystematic = false,bool inputIsYieldSystematic = false,TString inputFileName = "auau2014D0CorV2.root")
{
  D0CorV2Plotter *d0V2Plotter = new D0CorV2Plotter();
  d0V2Plotter->Init(inputFileName,inputIsSystematic,inputIsYieldSystematic);
  vector<vector<double> > hadronResult = d0V2Plotter->getHadronV2();
  vector<double> hadronV2Values = hadronResult[0];
  vector<vector<double> > fitResult = d0V2Plotter->fitMass();
  float xbin[7] = {0,1,2,3,4,5,10};
  TH1D *inputCandOverSignal = new TH1D("inputCandOverSignal","",6,xbin);
  for(int ipt=0;ipt<6;ipt++)
  {
    inputCandOverSignal->SetBinContent(ipt+1,fitResult[0][ipt]);
    inputCandOverSignal->SetBinError(ipt+1,fitResult[1][ipt]);
  }
  vector<vector<double> >d0V2;
  if(!inputIsSystematic || inputIsYieldSystematic)
  {
    cout<<"central value mode:"<<endl;
    d0V2.push_back((d0V2Plotter->getD0V2(hadronV2Values, inputCandOverSignal, 7))[0]);
  }
  else
  {
    cout<<"systematic mode:"<<endl;
    d0V2.push_back((d0V2Plotter->getD0V2(hadronV2Values, inputCandOverSignal, 7))[0]);
    d0V2.push_back((d0V2Plotter->getD0V2(hadronV2Values, inputCandOverSignal, 1))[0]);
    d0V2.push_back((d0V2Plotter->getD0V2(hadronV2Values, inputCandOverSignal, 2))[0]);
    d0V2.push_back((d0V2Plotter->getD0V2(hadronV2Values, inputCandOverSignal, 6))[0]);
  }
  delete d0V2Plotter;
  return d0V2;
}

vector<double> getYieldSys()
{
  cout<<"=================================Get Yield Systematic==================================="<<endl;
  vector<double> yieldSys;
  vector<vector<double> > yieldSysContainer;
  cout<<"first call"<<endl;
  yieldSysContainer.push_back( (getD0V2(false,false,"auau2014D0CorV2.root"))[0]);
  cout<<"second call"<<endl;
  yieldSysContainer.push_back( (getD0V2(false,true,"auau2014D0CorV2.root"))[0]);
  for(int i=0;i<6;i++)
    yieldSys.push_back(fabs(yieldSysContainer[1][i]-yieldSysContainer[0][i]));
  return yieldSys;
}

vector<double> getEfficiencySys()
{
  cout<<"=================================Get Efficiency Systematic==================================="<<endl;
  vector<double> efficiencySys;
  vector<vector<double> > efficiencySysContainer;
  efficiencySysContainer.push_back( (getD0V2(true,false,"auau2014D0CorV2.root"))[0]);
  efficiencySysContainer.push_back( (getD0V2(true,false,"no-eff.root"))[0]);
  for(int i=0;i<6;i++)
    efficiencySys.push_back(fabs(efficiencySysContainer[1][i]-efficiencySysContainer[0][i]));
  return efficiencySys;
}
vector<double> getBackgroundSys()
{
  cout<<"=================================Get Background Systematic==================================="<<endl;
  vector<double> backgroundSys;
  vector<vector<double> > backgroundSysContainer;
  vector<vector<double> > tempContainer;
  tempContainer.clear();
  tempContainer = getD0V2(true,false,"auau2014D0CorV2.root");
  for(int i=0;i<4;i++)
    backgroundSysContainer.push_back(tempContainer[i]);

  tempContainer.clear();
  tempContainer = getD0V2(true,false,"eff_50.root");
  for(int i=0;i<4;i++)
    backgroundSysContainer.push_back(tempContainer[i]);
  tempContainer.clear();
  tempContainer = getD0V2(true,false,"eff_150.root");
  for(int i=0;i<4;i++)
    backgroundSysContainer.push_back(tempContainer[i]);

  vector<double> bkgSysMatrix[6];

  for(int isample=1;isample<12;isample++)
    for(int ipt=0;ipt<6;ipt++)
      bkgSysMatrix[ipt].push_back(fabs(backgroundSysContainer[isample][ipt]-backgroundSysContainer[0][ipt]));

  for(int ipt=0;ipt<6;ipt++)
    backgroundSys.push_back((*max_element(bkgSysMatrix[ipt].begin(),bkgSysMatrix[ipt].end()))/sqrt(12));
      
  return backgroundSys;
}
