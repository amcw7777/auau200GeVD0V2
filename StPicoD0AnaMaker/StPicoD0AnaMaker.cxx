#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TLorentzVector.h"

#include "phys_constants.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoD0EventMaker/StPicoD0Event.h"
#include "StPicoD0EventMaker/StKaonPion.h"
#include "StPicoD0AnaMaker.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StCuts.h"
#include "../StPicoPrescales/StPicoPrescales.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
//////Refit include lib
#include "PhysicalConstants.h"
#include "StThreeVectorF.hh"
#include "StLorentzVectorD.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TFile.h"
#include "StEvent/StDcaGeometry.h"
//
#include <vector>
#include <stdio.h>
#include <time.h>
#include <algorithm>

ClassImp(StPicoD0AnaMaker)

  StPicoD0AnaMaker::StPicoD0AnaMaker(char const * name,char const * inputFilesList, 
      char const * outName,StPicoDstMaker* picoDstMaker,StRefMultCorr* grefmultCorrUtil): 
    StMaker(name),mPicoDstMaker(picoDstMaker),mPicoD0Event(NULL), mGRefMultCorrUtil(grefmultCorrUtil),
    mOutFileName(outName), mInputFileList(inputFilesList),mOutputFile(NULL), mChain(NULL), mEventCounter(0){}

Int_t StPicoD0AnaMaker::Init()
{
  mPicoD0Event = new StPicoD0Event();

  mChain = new TChain("T");
  std::ifstream listOfFiles(mInputFileList.Data());
  if (listOfFiles.is_open())
  {
    std::string file;
    while (getline(listOfFiles, file))
    {
      LOG_INFO << "StPicoD0AnaMaker - Adding :" << file <<endm;
      mChain->Add(file.c_str());
      LOG_INFO<<" Entries = "<<mChain->GetEntries()<< endm; 
    }
  }
  else
  {
    LOG_ERROR << "StPicoD0AnaMaker - Could not open list of files. ABORT!" << endm;
    return kStErr;
  }

  mPrescales = new StPicoPrescales(mycuts::prescalesFilesDirectoryName);

  mChain->GetBranch("dEvent")->SetAutoDelete(kFALSE);
  mChain->SetBranchAddress("dEvent", &mPicoD0Event);

  mOutputFile = new TFile(mOutFileName.Data(), "RECREATE");
  mOutputFile->cd();
  mHadronTuple = new TNtuple("mHadronTuple","","sum1:sin1:cos1:sum2:sin2:cos2:cent:reweight");
  etaPhi = new TH2F("etaPhi","",200,-6.29,6.29,100,0.5,1.5);
  etaPhi_D = new TH2F("etaPhiD","",100,-3.1416,3.1416,100,-2,2);
  etaPhi_Hadron = new TH2F("etaPhiHadron","",100,-3.1416,3.1416,100,-2,2);
  etaPhi_Hadron_all = new TH2F("etaPhiHadronAll","",100,-3.1416,3.1416,100,-2,2);
  dEtaDHadron = new TH1F("dEtaDHadron","",1000,0,10);
  hEtaD = new TH1F("hEtaD","",1000,0,10);
  hEtaHadron = new TH1F("hEtaHadron","",1000,0,10);
  vtxz = new TH1F("vtxz","",100,-10,10);
  TString flatten[5];
  flatten[0] = "v2";
  flatten[1] = "cosD";
  flatten[2] = "sinD";
  flatten[3] = "cosHadron";
  flatten[4] = "sinHadron";
  TString sb[8] = {"s1like","s3like","hSBlike","lSBlike","s1unlike","s3unlike","hSBunlike","lSBunlike"};
  float xbin[7] = {0,1,2,3,4,5,10};
  float binMass[2001];
  float binPhi[2001];
  for(int i=0;i<2001;i++)
    binPhi[i] = 0.005*i-5;
  for(int i=0;i<2001;i++)
    binMass[i] = 0.01*i;
  massPt = new TH2D("massPt","",2000,binMass,6,xbin);
  massPtLike = new TH2D("massPtLike","",250,0,2.5,100,0,10);
  massLike = new TH2D("massLike","",250,0,2.5,100,0,10);
  massLike->Sumw2();
  massUnlike = new TH2D("massUnlike","",2000,binMass,6,xbin);
  massUnlike->Sumw2();
  for(int i=0;i!=8;i++)
  {
    for(int k=0;k!=3;k++)
    {
      for(int j=0;j!=5;j++)
      {
        TString name = sb[i]+flatten[j]+Form("_%i",k);
        profV2[i][j][k] = new TProfile(name.Data(),"",6,xbin);
        profV2[i][j][k]->Sumw2();
      }
      TString weightName = sb[i]+Form("_%i_weigth",k);
      float xWeight[10];
      for(int ii=0;ii<10;ii++)
        xWeight[ii] = ii;
      v2Weight[i][k] = new TH2D(weightName.Data(),"",9,xWeight,6,xbin);
      v2Weight[i][k]->Sumw2();

      TString namehPhi = "hadronPhi_"+sb[i]+Form("_%i",k);
      TString nameDPhi = "DPhi_"+sb[i]+Form("_%i",k);
      hPhiHadron[i][k] = new TH2F(namehPhi.Data(),"",2000,binPhi,6,xbin);
      hPhiD[i][k]= new TH2F(nameDPhi.Data(),"",2000,binPhi,6,xbin);
      hPhiD[i][k]->Sumw2();
      hPhiHadron[i][k]->Sumw2();
    }
  }

  float ptbin1[12] = {0.225,0.375,0.525,0.675,0.825,0.975,1.12,1.27,1.42,1.58,1.73,1.88};
  float ptbin2[11];
  for(int i=0;i<11;i++)
    ptbin2[i] = 0.5*(ptbin1[i]+ptbin1[i+1]);
  for(int k=0;k<3;k++)
  {
    for(int i=0;i<5;i++)
    {
      hadronV2[i][k] = new TH1D(Form("hadron_%s_%i",flatten[i].Data(),k),"",9,0,9);
      hadronV2[i][k]->Sumw2();
      hadronV2_sum[i][k] = new TH1D(Form("hadronsum_%s_%i",flatten[i].Data(),k),"",9,0,9);
      hadronV2_sum[i][k]->Sumw2();
      for(int j=0;j<9;j++)
      {
        hadronV2_excl[i][j][k] = new TH1D(Form("hadron_%s_cent%i_%i",flatten[i].Data(),j,k),"",10,ptbin2);
        hadronV2_excl[i][j][k]->Sumw2();
        hadronV2_excl_sum[i][j][k] = new TH1D(Form("hadronsum_%s_cent%i_%i",flatten[i].Data(),j,k),"",10,ptbin2);
        hadronV2_excl_sum[i][j][k]->Sumw2();
      }
    }
  }

  double fitmean[6] = {1.85921,1.8633,1.86403,1.86475,1.86252,1.86534};
  double fitsigma[6] = {0.018139,0.0139476,0.0158346,0.0169282,0.0199567,0.0189131};
  ifstream ifs("efficiency.txt");
  for(int i=0; i<6; i++)
    for(int j=0; j<4; j++)
      ifs>>efficiency[j][i];
  for(int i=0;i<6;i++)//pt bin
  {
    for(int j=0;j<5;j++)//flatten
    {
      TString massName[2];
      massName[0] = Form("likeMass%i",i)+flatten[j];
      massName[1] = Form("unlikeMass%i",i)+flatten[j];
      for(int k=0;k<2;k++)
      {
        V2Mass[k][i][j] = new TProfile(massName[k].Data(),"",18,fitmean[i]-9*fitsigma[i],fitmean[i]+9*fitsigma[i]);
        V2Mass[k][i][j]->Sumw2();
      }
    }
  }

  mOutputFile->cd();


  // -------------- USER VARIABLES -------------------------
  mGRefMultCorrUtil = new StRefMultCorr("grefmult");

  return kStOK;
}
//-----------------------------------------------------------------------------
StPicoD0AnaMaker::~StPicoD0AnaMaker()
{
  /*  */
  delete mGRefMultCorrUtil;
}
//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Finish()
{
  LOG_INFO << " StPicoD0AnaMaker - writing data and closing output file " <<endm;
  fout.close();
  fout1.close();
  mOutputFile->cd();
  // save user variables here
  massPt->Write();
  massPtLike->Write();
  vtxz->Write();
  for(int i=0;i!=8;i++)
  {
    for(int k=0;k!=3;k++)
    {
      for(int j=0;j!=5;j++)
      {
        profV2[i][j][k]->Write();
      }
      v2Weight[i][k]->Write();
      hPhiHadron[i][k]->Write();
      hPhiD[i][k]->Write();
    }
  }
  for(int i=0;i<6;i++)
  {
    for(int j=0;j<5;j++)
    {
      for(int k=0;k<2;k++)
        V2Mass[k][i][j]->Write();
    }
  }
  massLike->Write();
  massUnlike->Write();
  for(int k=0;k<3;k++)
  {
    for(int i=0;i<5;i++)
    {
      hadronV2[i][k]->Write();
      hadronV2_sum[i][k]->Write();
    }
  }

  mOutputFile->Close();
  delete mPrescales;

  return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Make()
{
  readNextEvent();
  if (!mPicoDstMaker)
  {
    LOG_WARN << " StPicoD0AnaMaker - No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  picoDst = mPicoDstMaker->picoDst();

  if (!picoDst)
  {
    LOG_WARN << "StPicoD0AnaMaker - No PicoDst! Skip! " << endm;
    return kStWarn;
  }
  if(mPicoD0Event->runId() != picoDst->event()->runId() ||
      mPicoD0Event->eventId() != picoDst->event()->eventId())
  {
    LOG_ERROR <<" StPicoD0AnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!"<<endm;
    LOG_ERROR <<" StPicoD0AnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoD0Event are not in sync."<<endm;
    exit(1);
  }

  // -------------- USER ANALYSIS -------------------------
  TClonesArray const * aKaonPion = mPicoD0Event->kaonPionArray();


  StThreeVectorF pVtx(-999.,-999.,-999.);
  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  if(!(isGoodEvent()) || !event->isMinBias())//minBias trigger requires
  {
    LOG_WARN << " Not Good Event! Skip! " << endm;
    return kStWarn;
  }
  if(event) {
    pVtx = event->primaryVertex();
  }
  vtxz->Fill(pVtx.z());
  if(!mGRefMultCorrUtil) {
    LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
    return kStWarn;
  }
  mGRefMultCorrUtil->init(picoDst->event()->runId());
  mGRefMultCorrUtil->initEvent(picoDst->event()->grefMult(),pVtx.z(),picoDst->event()->ZDCx()) ;

  for(int k=0;k<3;k++)
    getHadronCorV2(k);//Calculate Hadron V2
  int centrality  = mGRefMultCorrUtil->getCentralityBin9();
  int centBin = 0;
  if(centrality>=7) centBin=1;
  else if(centrality>=4)  centBin=2;
  else centBin=3;

  double reweight = mGRefMultCorrUtil->getWeight();
  for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx)
  {
    StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
    StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
    StPicoTrack const* pion = picoDst->track(kp->pionIdx());

    if (!isGoodTrack(kaon) || !isGoodTrack(pion)) continue;
    if (!isTpcPion(pion)) continue;
    bool tpcKaon = isTpcKaon(kaon,&pVtx);
    float kBeta = getTofBeta(kaon,&pVtx);
    bool tofAvailable = kBeta>0;
    bool tofKaon = tofAvailable && isTofKaon(kaon,kBeta);
    bool goodKaon = (tofAvailable && tofKaon) || (!tofAvailable && tpcKaon);
    if(!goodKaon) continue;
    int charge=0;
    float d0Pt = kp->pt();
    double dMass = kp->m();
    if(d0Pt>10) continue;
    int fitindex = 5;
    if(d0Pt<5)
      fitindex = static_cast<int>(d0Pt);
    double reweight_eff = (efficiency[0][fitindex]/efficiency[centBin][fitindex]);

    if((charge=isD0Pair150(kp))!=0 )
    {
      getCorV2(idx,reweight*reweight_eff);//Fill D-hadron 2PC v2 plots
      if(charge==-1)
        massPt->Fill(dMass,d0Pt,reweight*reweight_eff);
      if(charge>0)
        massPtLike->Fill(dMass,d0Pt,reweight*reweight_eff);

    }//D loop

  }

  return kStOK;
}
//-----------------------------------------------------------------------------

bool StPicoD0AnaMaker::isGoodPair(StKaonPion const* const kp) const
{
  if(!kp) return false;

  StPicoTrack const* kaon = mPicoDstMaker->picoDst()->track(kp->kaonIdx());
  StPicoTrack const* pion = mPicoDstMaker->picoDst()->track(kp->pionIdx());

  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  StThreeVectorF pVtx(-999.,-999.,-999.);
  pVtx = event->primaryVertex();
  int charge = kaon->charge() * pion->charge();
  bool pairCuts = kp->m()>1.6 && kp->m()<2.1 &&
    charge==-1;

  return (isTpcKaon(kaon,&pVtx) && isTpcPion(pion) && 
      pairCuts);
}


int StPicoD0AnaMaker::isD0Pair(StKaonPion const* const kp) const
{

  StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
  StPicoTrack const* pion = picoDst->track(kp->pionIdx());
  bool pairCuts = false;
  if(kp->pt()<1)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0061 &&
      kp->pionDca() > 0.0110 && kp->kaonDca() > 0.0103 &&
      kp->dcaDaughters() < 0.0084 && kp->decayLength()>0.0145;  
  }
  else if(kp->pt()<2)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0049 &&
      kp->pionDca() > 0.0111 && kp->kaonDca() > 0.0091 &&
      kp->dcaDaughters() < 0.0066 && kp->decayLength()>0.0181;  
  }
  else if(kp->pt()<3)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0038 &&
      kp->pionDca() > 0.0086 && kp->kaonDca() > 0.0095 &&
      kp->dcaDaughters() < 0.0057 && kp->decayLength()>0.0212;  
  }
  else if(kp->pt()<5)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0038 &&
      kp->pionDca() > 0.0081 && kp->kaonDca() > 0.0079 &&
      kp->dcaDaughters() < 0.0050 && kp->decayLength()>0.0247;  
  }
  else 
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0040 &&
      kp->pionDca() > 0.0062 && kp->kaonDca() > 0.0058 &&
      kp->dcaDaughters() < 0.0060 && kp->decayLength()>0.0259;  
  }

  int charge = kaon->charge() * pion->charge();
  if(charge>0)
    charge = kaon->charge()>0 ? 1:2;


  if(pairCuts)
    return charge;
  else
    return 0;
}
/*
*/

int StPicoD0AnaMaker::isD0Pair50(StKaonPion const* const kp) const
{

  StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
  StPicoTrack const* pion = picoDst->track(kp->pionIdx());
  bool pairCuts = false;
  if(kp->pt()<1)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0044 &&
      kp->pionDca() > 0.0120 && kp->kaonDca() > 0.0119 &&
      kp->dcaDaughters() < 0.0069 && kp->decayLength()>0.0144;  
  }
  else if(kp->pt()<2)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0036 &&
      kp->pionDca() > 0.0102 && kp->kaonDca() > 0.0110 &&
      kp->dcaDaughters() < 0.0048 && kp->decayLength()>0.0204;  
  }
  else if(kp->pt()<3)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0031 &&
      kp->pionDca() > 0.0118 && kp->kaonDca() > 0.0109 &&
      kp->dcaDaughters() < 0.0044 && kp->decayLength()>0.0242;  
  }
  else if(kp->pt()<5)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0026 &&
      kp->pionDca() > 0.0109 && kp->kaonDca() > 0.0106 &&
      kp->dcaDaughters() < 0.0049 && kp->decayLength()>0.0245;  
  }
  else 
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0032 &&
      kp->pionDca() > 0.0096 && kp->kaonDca() > 0.0080 &&
      kp->dcaDaughters() < 0.0047 && kp->decayLength()>0.0300;  
  }

  int charge = kaon->charge() * pion->charge();
  if(charge>0)
    charge = kaon->charge()>0 ? 1:2;


  if(pairCuts)
    return charge;
  else
    return 0;
}

int StPicoD0AnaMaker::isD0Pair150(StKaonPion const* const kp) const
{

  StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
  StPicoTrack const* pion = picoDst->track(kp->pionIdx());
  bool pairCuts = false;
  if(kp->pt()<1)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0072 &&
      kp->pionDca() > 0.0092 && kp->kaonDca() > 0.0105 &&
      kp->dcaDaughters() < 0.0077 && kp->decayLength()>0.0110;  
  }
  else if(kp->pt()<2)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0053 &&
      kp->pionDca() > 0.0078 && kp->kaonDca() > 0.0068 &&
      kp->dcaDaughters() < 0.0078 && kp->decayLength()>0.0168;  
  }
  else if(kp->pt()<3)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0047 &&
      kp->pionDca() > 0.0086 && kp->kaonDca() > 0.0080 &&
      kp->dcaDaughters() < 0.0074 && kp->decayLength()>0.0187;  
  }
  else if(kp->pt()<5)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0042 &&
      kp->pionDca() > 0.0065 && kp->kaonDca() > 0.0066 &&
      kp->dcaDaughters() < 0.0068 && kp->decayLength()>0.0199;  
  }
  else 
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0062 &&
      kp->pionDca() > 0.0047 && kp->kaonDca() > 0.0041 &&
      kp->dcaDaughters() < 0.0066 && kp->decayLength()>0.0180;  
  }

  int charge = kaon->charge() * pion->charge();
  if(charge>0)
    charge = kaon->charge()>0 ? 1:2;


  if(pairCuts)
    return charge;
  else
    return 0;
}
int StPicoD0AnaMaker::isD0PairOld(StKaonPion const* const kp) const
{

  StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
  StPicoTrack const* pion = picoDst->track(kp->pionIdx());
  bool pairCuts = false;
  if(kp->pt()<1)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0062 &&
      kp->pionDca() > 0.0109 && kp->kaonDca() > 0.0123 &&
      kp->dcaDaughters() < 0.0082 && kp->decayLength()>0.0149;  
  }
  else if(kp->pt()<2)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0047 &&
      kp->pionDca() > 0.0108 && kp->kaonDca() > 0.0097 &&
      kp->dcaDaughters() < 0.0070 && kp->decayLength()>0.0205;  
  }
  else if(kp->pt()<3)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0040 &&
      kp->pionDca() > 0.0100 && kp->kaonDca() > 0.0091 &&
      kp->dcaDaughters() < 0.0056 && kp->decayLength()>0.0216;  
  }
  else if(kp->pt()<5)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0041 &&
      kp->pionDca() > 0.0074 && kp->kaonDca() > 0.0075 &&
      kp->dcaDaughters() < 0.0065 && kp->decayLength()>0.0233;  
  }
  else 
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0042 &&
      kp->pionDca() > 0.067 && kp->kaonDca() > 0.0053 &&
      kp->dcaDaughters() < 0.0065 && kp->decayLength()>0.0282;  
  }

  int charge = kaon->charge() * pion->charge();


  if(pairCuts)
    return charge;
  else
    return 0;
}



bool StPicoD0AnaMaker::isGoodEvent()
{
  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  return (event->triggerWord() & mycuts::triggerWord) &&
    fabs(event->primaryVertex().z()) < mycuts::vz &&
    fabs(event->primaryVertex().z() - event->vzVpd()) < mycuts::vzVpdVz;
  //  return event->triggerWord() & mycuts::triggerWord;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodTrack(StPicoTrack const * const trk) const
{
  // Require at least one hit on every layer of PXL and IST.
  // It is done here for tests on the preview II data.
  // The new StPicoTrack which is used in official production has a method to check this
  return trk->gPt() > mycuts::minPt && trk->nHitsFit() >= mycuts::nHitsFit && trk->isHFTTrack();
  //return  trk->nHitsFit() >= mycuts::nHitsFit;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodHadron(StPicoTrack const * const trk) const
{
  //return trk->pMom().perp() > mycuts::hadronPtMin &&trk->pMom().perp() < mycuts::hadronPtMax && trk->nHitsFit() >= mycuts::nHitsFit &&fabs(trk->pMom().pseudoRapidity())<1.&&fabs(trk->nSigmaElectron())>3 && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52;
  return trk->pMom().perp() > mycuts::hadronPtMin &&trk->pMom().perp() < mycuts::hadronPtMax && trk->nHitsFit() >= mycuts::nHitsFit &&fabs(trk->pMom().pseudoRapidity())<1. && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTpcPion(StPicoTrack const * const trk) const
{
  return fabs(trk->nSigmaPion()) < mycuts::nSigmaPion;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTpcKaon(StPicoTrack const * const trk, StThreeVectorF const* const pVtx) const
{
  return fabs(trk->nSigmaKaon()) < mycuts::nSigmaKaon;
  //      || tofKaon;
}
//-----------------------------------------------------------------------------
float StPicoD0AnaMaker::getTofBeta(StPicoTrack const * const trk, StThreeVectorF const* const pVtx) const
{

  int index2tof = trk->bTofPidTraitsIndex();

  float beta = std::numeric_limits<float>::quiet_NaN();

  if(index2tof >= 0)
  {
    StPicoBTofPidTraits *tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);

    if(tofPid)
    {
      beta = tofPid->btofBeta();

      if (beta < 1e-4)
      {
        StThreeVectorF const btofHitPos = tofPid->btofHitPos();
        StPhysicalHelixD helix = trk->helix();

        float L = tofPathLength(pVtx, &btofHitPos, helix.curvature());
        float tof = tofPid->btof();
        if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
        else beta = std::numeric_limits<float>::quiet_NaN();
      }
    }
  }

  return beta;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTofKaon(StPicoTrack const * const trk, float beta) const
{
  bool tofKaon = false;

  if(beta>0)
  {
    double ptot = trk->dcaGeometry().momentum().mag();
    float beta_k = ptot/sqrt(ptot*ptot+M_KAON_PLUS*M_KAON_PLUS);
    tofKaon = fabs(1/beta - 1/beta_k) < mycuts::kTofBetaDiff ? true : false;
  }

  return tofKaon;
}


bool StPicoD0AnaMaker::getCorHadron(float eta,vector<float> &hadronsPhi, vector<unsigned int> index1, float phi, float etaCut) 
{
  for(unsigned int i=0;i<picoDst->numberOfTracks();++i)
  {
    StPicoTrack const* hadron = picoDst->track(i);
    if(!hadron)  continue; 
    if(hadron->pMom().perp()<0.2) continue;
    etaPhi_Hadron_all->Fill(hadron->pMom().phi(),hadron->pMom().pseudoRapidity());
    vector<unsigned int>::iterator it_index;
    it_index = find(index1.begin(),index1.end(),i);
    if(it_index!=index1.end())  continue;
    if(!isGoodHadron(hadron)) continue;
    float dEta = fabs(hadron->pMom().pseudoRapidity() - eta);
    float dPhi = (hadron->pMom().phi() - phi);
    if(etaCut<0.001)
    {
      dEtaDHadron->Fill(dEta);
      hEtaD->Fill(eta);
      hEtaHadron->Fill(hadron->pMom().pseudoRapidity());
    }
    //if(dPhi>3.1416) dPhi = 2*3.1416-dPhi;
    if(dEta< etaCut|| dEta > mycuts::corDetaMax)  continue;
    etaPhi->Fill(dPhi,dEta);
    etaPhi_Hadron->Fill(hadron->pMom().phi(),hadron->pMom().pseudoRapidity());
    hadronsPhi.push_back(hadron->pMom().phi());
  }
  //  fixPhi(hadronsPhi);
  return true;

}

float StPicoD0AnaMaker::sumCos(float phi,vector<float> &hadronsPhi) 
{
  float sumOfCos = 0;
  for(unsigned int i=0;i<hadronsPhi.size();++i)
  {
    sumOfCos += cos(2*(phi-hadronsPhi[i]));
  }
  return sumOfCos;
}

bool StPicoD0AnaMaker::fixPhi(vector<float> &phi) 
{
  if(phi.size() == 0) return false;
  float sumPhi = 0;
  for(unsigned int i=0;i<phi.size();i++)
    sumPhi+=phi[i];
  float meanPhi = sumPhi/phi.size();
  for(unsigned int i=0;i<phi.size();i++)
    phi[i] = phi[i]-meanPhi;  
  return true;
}


bool StPicoD0AnaMaker::getHadronCorV2(int idxGap)
{
  double etaGap[3] = {0,0.15,0.5};
  double mEtaGap = etaGap[idxGap];
  float hadronFill[7] = {0};
  const double reweight = mGRefMultCorrUtil->getWeight();
  int centrality  = mGRefMultCorrUtil->getCentralityBin9();
  for(unsigned int i=0;i<picoDst->numberOfTracks();++i)
  {
    StPicoTrack const* hadron = picoDst->track(i);
    if(hadron->pMom().perp()<0.2) continue;
    if(!isGoodHadron(hadron)) continue;
    float etaHadron = hadron->pMom().pseudoRapidity();
    float phiHadron = hadron->pMom().phi();
    if(etaHadron<-0.5*mEtaGap)//backward sample 
    {
      hadronFill[0]++;
      hadronFill[1] += sin(2 * phiHadron);
      hadronFill[2] += cos(2 * phiHadron);
    }			
    if(etaHadron>0.5*mEtaGap)//forward sample
    {
      hadronFill[3]++;
      hadronFill[4] += sin(2 * phiHadron);
      hadronFill[5] += cos(2 * phiHadron);
    }			
  }
  hadronFill[6] = centrality;
  hadronFill[7] = reweight;
  //mHadronTuple->Fill(hadronFill);
  if(hadronFill[0]==0 || hadronFill[3]==0)
    return false; 
  double temp = (hadronFill[1]*hadronFill[4]+hadronFill[2]*hadronFill[5]);
  hadronV2[0][idxGap]->Fill(centrality,temp*reweight);
  hadronV2[1][idxGap]->Fill(centrality,hadronFill[2]*reweight);
  hadronV2[2][idxGap]->Fill(centrality,hadronFill[1]*reweight);
  hadronV2[3][idxGap]->Fill(centrality,hadronFill[5]*reweight);
  hadronV2[4][idxGap]->Fill(centrality,hadronFill[4]*reweight);
  hadronV2_sum[0][idxGap]->Fill(centrality,hadronFill[0]*hadronFill[3]*reweight);
  hadronV2_sum[1][idxGap]->Fill(centrality,hadronFill[0]*reweight);
  hadronV2_sum[2][idxGap]->Fill(centrality,hadronFill[0]*reweight);
  hadronV2_sum[3][idxGap]->Fill(centrality,hadronFill[3]*reweight);
  hadronV2_sum[4][idxGap]->Fill(centrality,hadronFill[3]*reweight);
  //    StPicoTrack const* hadron = picoDst->track(i);
  //  hadronV2_excl[0][centrality]->Fill(hadron->pMom().perp(),temp*reweight);
  //  hadronV2_excl[1][centrality]->Fill(hadron->pMom().perp(),hadronFill[2]*reweight);
  //  hadronV2_excl[2][centrality]->Fill(hadron->pMom().perp(),hadronFill[1]*reweight);
  //  hadronV2_excl[3][centrality]->Fill(hadron->pMom().perp(),hadronFill[5]*reweight);
  //  hadronV2_excl[4][centrality]->Fill(hadron->pMom().perp(),hadronFill[4]*reweight);
  //  hadronV2_excl_sum[0][centrality]->Fill(hadron->pMom().perp(),hadronFill[0]*hadronFill[3]*reweight);
  //  hadronV2_excl_sum[1][centrality]->Fill(hadron->pMom().perp(),hadronFill[0]*reweight);
  //  hadronV2_excl_sum[2][centrality]->Fill(hadron->pMom().perp(),hadronFill[0]*reweight);
  //  hadronV2_excl_sum[3][centrality]->Fill(hadron->pMom().perp(),hadronFill[3]*reweight);
  //  hadronV2_excl_sum[4][centrality]->Fill(hadron->pMom().perp(),hadronFill[3]*reweight);
  return true;
}

//bool StPicoD0AnaMaker::getCorV2(int idxCand,double etaGap, int &flagCor, double weight)
bool StPicoD0AnaMaker::getCorV2(int idxCand,double weight)
{
  int centrality  = mGRefMultCorrUtil->getCentralityBin9();
  TClonesArray const * aKaonPion = mPicoD0Event->kaonPionArray();
  StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idxCand);
  StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
  StPicoTrack const* pion = picoDst->track(kp->pionIdx());
  int charge = kaon->charge() * pion->charge();
  double dMass = kp->m();
  int centBin = 0;
  if(centrality>=7) centBin=1;
  else if(centrality>=4)  centBin=2;
  else centBin=3;
  //double hadronv2[9] = {0.0557348,0.0630083,0.0704378,0.0742203,0.0725183,0.0644258,0.0494283,0.0349057,0.024481};
  double hadronv2[9] = {1,1,1,1,1,1,1,1,1};


  if(kp->pt()>10) return false;
  int ptIdx = 5;
  if(kp->pt()<5)
    ptIdx = static_cast<int>(kp->pt());
  double fitmean[6] = {1.85921,1.8633,1.86403,1.86475,1.86252,1.86534};
  double fitsigma[6] = {0.018139,0.0139476,0.0158346,0.0169282,0.0199567,0.0189131};
  double mean = fitmean[ptIdx];
  double sigma = fitsigma[ptIdx];
  bool fillSB[8];
  fillSB[0] =  (charge>0)&& (dMass>(mean-1*sigma)) &&  (dMass<(mean+1*sigma));
  fillSB[1] =  (charge>0)&& (dMass>(mean-3*sigma)) &&  (dMass<(mean+3*sigma));
  fillSB[2] =  (charge>0) && (((dMass>(mean+4*sigma)) &&  (dMass<(mean+9*sigma))) ||((dMass>(mean-9*sigma)) &&  (dMass<(mean-4*sigma))));
  fillSB[4] = (charge==-1)&& (((dMass>(mean+5.5*sigma)) &&  (dMass<(mean+7.5*sigma))) ||((dMass>(mean-7.5*sigma)) &&  (dMass<(mean-5.5*sigma))));
  fillSB[5] = (charge==-1)&& (dMass>(mean-3*sigma)) &&  (dMass<(mean+3*sigma));
  fillSB[6] = (charge==-1)&& (((dMass>(mean+4*sigma)) &&  (dMass<(mean+9*sigma))) ||((dMass>(mean-9*sigma)) &&  (dMass<(mean-4*sigma))));
  fillSB[3] = fillSB[1] || fillSB[2];
  fillSB[7] = fillSB[1] || fillSB[2] || fillSB[6];
  double etaGap[3] = {0,0.15,0.5};

  for(int k=0;k<3;k++)
  {
    double corFill[7] = {0};
    corFill[0] = 1 ;
    corFill[1] = sin(2* kp->phi())/sqrt(hadronv2[centBin]);
    corFill[2] = cos(2* kp->phi())/sqrt(hadronv2[centBin]);
    int chargeIdx = charge>0 ? 0:1;
    if(k==0)
    {
      V2Mass[chargeIdx][ptIdx][1]->Fill(kp->m(),corFill[1],weight);
      V2Mass[chargeIdx][ptIdx][2]->Fill(kp->m(),corFill[2],weight);
    }
    for(int j=0;j<8;j++)
    {
      if(fillSB[j])
      {
        //TH1D *pPhiD = hPhiD[j][k]->ProjectionX("dPhi",ptIdx+1,ptIdx+1);
        //double corDPhi = 1./pPhiD->GetBinContent(pPhiD->FindBin(kp->phi()))*(pPhiD->Integral())/2000;
        profV2[j][1][k]->Fill(kp->pt(),corFill[1]/sqrt(hadronv2[centBin]),weight);
        profV2[j][2][k]->Fill(kp->pt(),corFill[2]/sqrt(hadronv2[centBin]),weight);
        hPhiD[j][k]->Fill(kp->phi(),kp->pt(),weight);
      }
    }
    for(unsigned int i=0; i<picoDst->numberOfTracks();i++)
    {
      StPicoTrack const* hadron = picoDst->track(i);
      if(hadron->pMom().perp()<0.2) continue;
      if(!isGoodHadron(hadron)) continue;
      float etaHadron = hadron->pMom().pseudoRapidity();
      float phiHadron = hadron->pMom().phi();
      if(!isEtaGap(kp->eta(),etaGap[k],etaHadron))  continue;
      //if(etaHadron*kp->eta() > 0) continue;
      //TH1D *pPhiHadron = hPhiHadron[j][k]->ProjectionX("hPhi",ptIdx+1,ptIdx+1);
      //double corhPhi = 1./pPhiHadron->GetBinContent(pPhiHadron->FindBin(phiHadron));
      corFill[3]++;
      //corFill[6]+= corhPhi;
      corFill[4] += sin(2 * phiHadron)/sqrt(hadronv2[centBin]);
      corFill[5] += cos(2 * phiHadron)/sqrt(hadronv2[centBin]);
      if(k==0)
      {
        V2Mass[chargeIdx][ptIdx][3]->Fill(kp->m(),sin(2*phiHadron)/sqrt(hadronv2[centBin]),weight);
        V2Mass[chargeIdx][ptIdx][4]->Fill(kp->m(),cos(2*phiHadron)/sqrt(hadronv2[centBin]),weight);
      }
      for(int j=0;j<8;j++)
      {
        if(fillSB[j])
        {
          //TH1D *pPhiHadron = hPhiHadron[j][k]->ProjectionX("hPhi",ptIdx+1,ptIdx+1);
          //double corhPhi = 1./pPhiHadron->GetBinContent(pPhiHadron->FindBin(phiHadron))*(pPhiHadron->Integral())/2000;
          profV2[j][3][k]->Fill(kp->pt(),sin(2*phiHadron)/sqrt(hadronv2[centBin]),weight);
          profV2[j][4][k]->Fill(kp->pt(),cos(2*phiHadron)/sqrt(hadronv2[centBin]),weight);
          hPhiHadron[j][k]->Fill(phiHadron,kp->pt(),weight);
          //cout<<"wtf = "<<corhPhi<<endl;
        }
      }
    }// Loop over charged tracks
    if(corFill[3]<=0) return false;
    double cumulant = (corFill[1]*corFill[4]+corFill[2]*corFill[5])/(corFill[3]); 
    if(k==0)
    {
      if(charge<0)  massUnlike->Fill(kp->m(),kp->pt(),weight);
      if(charge>0)  massLike->Fill(kp->m(),kp->pt(),weight);
      V2Mass[chargeIdx][ptIdx][0]->Fill(kp->m(),cumulant, corFill[3]*weight);
    }
    for(int j=0;j<8;j++)
    {
      if(fillSB[j])
      {
        //TH1D *pPhiHadron = hPhiHadron[j][k]->ProjectionX("hPhi",ptIdx+1,ptIdx+1);
        //double corhPhi = 1./pPhiHadron->GetBinContent(pPhiHadron->FindBin(phiHadron));
        //TH1D *pPhiD = hPhiD[j][k]->ProjectionX("dPhi",ptIdx+1,ptIdx+1);
        //double corDPhi = 1./pPhiD->GetBinContent(pPhiD->FindBin(kp->phi()))*(pPhiD->Integral())/2000;
        profV2[j][0][k]->Fill(kp->pt(),cumulant, corFill[3]*weight);
        v2Weight[j][k]->Fill(centrality,kp->pt(),corFill[3]*weight);
      }
    }
  }//Loop over different eta gap (k)
  return true;
}

bool StPicoD0AnaMaker::isEtaGap(double dEta,double mGap,double hEta)
{
  if(mGap == 0) return true;
  double range =  2. - mGap*2;
  /*
     if(dEta>mGap/2)
     return hEta<(dEta-mGap) && hEta>(dEta-mGap-range);
     else if(dEta<-1.*mGap/2)
     return hEta>(dEta+mGap) && hEta<(dEta+mGap+range);
     else 
     return ((hEta>(dEta+mGap) && hEta<(dEta+mGap+range/2)) || (hEta<(dEta-mGap) && hEta>(dEta-mGap-range/2)));
     */
  if(dEta> (1.-2*mGap))
    return hEta<(dEta-mGap) && hEta>(dEta-mGap-range);
  else if(dEta<(-1.+2*mGap))
    return hEta>(dEta+mGap) && hEta<(dEta+mGap+range);
  else 
    return (hEta>(dEta+mGap) || hEta<(dEta-mGap));
}







