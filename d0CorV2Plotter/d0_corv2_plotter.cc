#include "fstream"
#include "iostream"
#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "Math/MinimizerOptions.h"
#include <vector>

#include "d0_corv2_plotter.h"
using namespace std;

void D0CorV2Plotter::Init(TString inputFileName,bool isSystematic,bool isYieldSystematic)
{
  mInputFileName = inputFileName;
  mIsYieldSystematic = isYieldSystematic;
  mIsSystematic = isSystematic;
  mIndexEtaGap = 1;
  mFlattenNames[0] = "v2";
  mFlattenNames[1] = "cosD";
  mFlattenNames[2] = "sinD";
  mFlattenNames[3] = "cosHadron";
  mFlattenNames[4] = "sinHadron";
	mInputFile= new TFile(mInputFileName.Data());
  if(!mIsSystematic)
    mLog.open("toM.txt");
  else
    mLog.open("systematicLog.txt");
}

vector<vector<double> > D0CorV2Plotter::getD0V2(vector<double> &hadronV2, TH1D *candOverSignal, int indexBkg) 
{
	TProfile *inputProfiles[8][5]; // i is different mass region; j is for flattening
	TH1D *inputHistograms[8][5];
	TH1D *flattenV2[8];
	TString sidebandNames[8] = {"s1like","s3like","hSBlike","lSBlike","s1unlike","s3unlike","hSBunlike","lSBunlike"};
	TH2D *v2Weight[8];  // weight of cumulants
	for(int i=0;i<8;i++)//different mass region
	{
		TString weightName = Form("%s_%i_weigth",sidebandNames[i].Data(),mIndexEtaGap);// 'weigth' is spelled wrong in source code
		v2Weight[i] = (TH2D *)mInputFile->Get(weightName.Data())->Clone(weightName.Data());
		for(int j=0;j<5;j++)
		{
			// TString name = sidebandNames[i]+mFlattenNames[j]+"_1";
			// TString nameV2 = sidebandNames[i]+mFlattenNames[j]+"1";
			TString profileName = Form("%s%s_%i",sidebandNames[i].Data(),mFlattenNames[j].Data(),mIndexEtaGap);
			TString histogramName = Form("%s%s%i",sidebandNames[i].Data(),mFlattenNames[j].Data(),mIndexEtaGap);
			inputProfiles[i][j] = (TProfile *)mInputFile->Get(profileName.Data())->Clone(profileName.Data());
			inputHistograms[i][j] = inputProfiles[i][j]->ProjectionX(histogramName);
		}
		inputHistograms[i][1]->Multiply(inputHistograms[i][3]); // cosD * cos Hadron
		inputHistograms[i][2]->Multiply(inputHistograms[i][4]); // sinD * sin Hadron
		for(int ipt=1;ipt<6;ipt++)
		{
			inputHistograms[i][1]->SetBinError(ipt,0);//The flattening terms have no stat. error
			inputHistograms[i][2]->SetBinError(ipt,0);
		}
		inputHistograms[i][0]->Add(inputHistograms[i][1],-1);
		inputHistograms[i][0]->Add(inputHistograms[i][2],-1);
		flattenV2[i] = (TH1D *)inputHistograms[i][0]->Clone(sidebandNames[i].Data());// V2 after phi flattening
    // Divide hardron v2 with weight of every cumulant
    // pT and centrality exclusively
		for(int ipt=0;ipt<6;ipt++)
		{
			double hadronv2Weight = 0;
			double dv2Weight = 0;
			for(int icent=0;icent<9;icent++)
			{		
				double mWeight = v2Weight[i]->GetBinContent(icent+1,ipt+1);
				hadronv2Weight += mWeight*hadronV2[icent];
				dv2Weight += mWeight;
			}
			double dv2 = inputHistograms[i][0]->GetBinContent(ipt+1)*dv2Weight/hadronv2Weight;
			double dv2_err = inputHistograms[i][0]->GetBinError(ipt+1)*dv2Weight/hadronv2Weight;
			flattenV2[i]->SetBinContent(ipt+1,dv2);
			flattenV2[i]->SetBinError(ipt+1,dv2_err);
		}
	}//mass region loop done	
/////// Set candidate and background sample
  mLog<<"============Foregroundd V2  ...========================="<<endl;
  mLog<<"pT\t\t\t\tForeground v2\t\t\t\tstat.error"<<endl;
	TH1D *candV2 = (TH1D *)flattenV2[5]->Clone("candV2");
  for(int ipt=1;ipt<7;ipt++)
  {
    double v2Value = candV2->GetBinContent(ipt);
    double v2Error= candV2->GetBinError(ipt);
    mLog<<ipt<<"\t\t\t\t"<<v2Value<<"\t\t\t\t"<<v2Error<<endl;
  }
  mLog<<"============Background V2  ...========================="<<endl;
  mLog<<"pT\t\t\t\tBackground v2\t\t\t\tstat.error"<<endl;
	TH1D *backgroundV2= (TH1D *)flattenV2[indexBkg]->Clone("backgroundV2");
  for(int ipt=1;ipt<7;ipt++)
  {
    double v2Value = backgroundV2->GetBinContent(ipt);
    double v2Error= backgroundV2->GetBinError(ipt);
    mLog<<ipt<<"\t\t\t\t"<<v2Value<<"\t\t\t\t"<<v2Error<<endl;
  }
// Background extraction
	// TH1D *signalV2 = (TH1D *)flattenV2[5]->Clone("signalV2");
	TH1D *signalV2 = getSignalV2(candOverSignal,backgroundV2,candV2);
  vector<double> d0V2Values;
  vector<double> d0V2Errors;
  mLog<<"============D0 v2  ...========================="<<endl;
  mLog<<"pT\t\t\t\tD0 v2\t\t\t\tD0 v2 stat.error"<<endl;
  for(int ipt=1;ipt<7;ipt++)
  {
    double v2Value = signalV2->GetBinContent(ipt);
    double v2Error= signalV2->GetBinError(ipt);
    d0V2Values.push_back(v2Value);
    d0V2Errors.push_back(v2Error);
    mLog<<ipt<<"\t\t\t\t"<<v2Value<<"\t\t\t\t"<<v2Error<<endl;
  }
  vector<vector<double> > d0V2Result;
  d0V2Result.push_back(d0V2Values);
  d0V2Result.push_back(d0V2Errors);
  for(int i=0;i<5;i++)
    for(int j=0;j<8;j++)
      delete inputHistograms[j][i];
  return d0V2Result;
}

vector<vector<double> > D0CorV2Plotter::fitMass()
{
  mLog<<"============Start to fit signal....========================="<<endl;
  TH2D *massPtUnlike = (TH2D *)mInputFile->Get("massUnlike")->Clone("massPtUnlike");
	TH2D *massPtLike = (TH2D *)mInputFile->Get("massLike")->Clone("massPtLike");
	massPtLike->RebinY(10);
	TH1D *massUnlike[6];
	TH1D *massLike[6];
	float xbin[7] = {0,1,2,3,4,5,10};
	TH1D *signalHistogram = new TH1D("signal","",6,xbin);
	TH1D *candidateHistogram = new TH1D("cand","",6,xbin);
	TCanvas *massCheck = new TCanvas();
	massCheck->Divide(3,2);
  mLog<<"pT bin\t\t\t\t#candidate\t\t\t\t\t\t#signal"<<endl;
	for(int ipt=0;ipt<6;ipt++)
	{
		massUnlike[ipt] = massPtUnlike->ProjectionX(Form("massUnlike_%i",ipt),ipt+1,ipt+1);	//unlike sign binning is 0,1,2..,5,10
		if(ipt!=5)
			massLike[ipt] = massPtLike->ProjectionX(Form("massLike_%i",ipt),ipt+1,ipt+1);	// like sign binnning is 0,1,2...,8,9,10
		else
			massLike[ipt] = massPtLike->ProjectionX(Form("massLike_%i",ipt),ipt+1,ipt+5);	
		double fitResult[3];
		massCheck->cd(ipt+1);
		massUnlike[ipt]->Draw();
		pair<double,double> sig_yield = fit_hist(massUnlike[ipt],massCheck,ipt+1,3,fitResult);
		massCheck->cd(ipt+1);
		massLike[ipt]->Draw("pe,same");
		massLike[ipt]->SetLineColor(4);
		//fout<<"i = "<<i<<endl;
		
		signalHistogram->SetBinContent(ipt+1,fitResult[0]);
		signalHistogram->SetBinError(ipt+1,sig_yield.second);
		candidateHistogram->SetBinContent(ipt+1,fitResult[1]);
		candidateHistogram->SetBinError(ipt+1,sqrt(fitResult[1]));
		mLog<<ipt<<"\t\t\t\t"<<fitResult[0]<<"+/-"<<sig_yield.second<<"\t\t\t\t\t\t"<<fitResult[1]<<"+/-"<<fitResult[2]<<endl;
	}
/////calculate gamma = #candidates / #signal
	TH1D *bkgHistogram = (TH1D *)candidateHistogram->Clone("bkgHistogram");
	bkgHistogram->Add(signalHistogram,-1);
	TH1D *candOverSignal = (TH1D *)candidateHistogram->Clone("candOverSignal");
	candOverSignal->Divide(candidateHistogram,signalHistogram,1,1);
  vector<double> candOverSignalValues;
  vector<double> candOverSignalErrors;
  for(int ipt=1;ipt<7;ipt++)
  {
    double cOSValue = candOverSignal->GetBinContent(ipt);
    double cOSError = candOverSignal->GetBinError(ipt);
    candOverSignalValues.push_back(cOSValue);
    candOverSignalErrors.push_back(cOSError);
  }
  vector<vector<double> > candOverSignalResult;
  candOverSignalResult.push_back(candOverSignalValues);
  candOverSignalResult.push_back(candOverSignalErrors);
  return candOverSignalResult;
}


vector<vector<double> > D0CorV2Plotter::getHadronV2()
{
  mLog<<"============Hadron v2 ...========================="<<endl;
	TFile *v2File = new TFile(mInputFileName.Data());
	TH1D *hadronV2Histogram[5][2];
	for(int iflatten=0;iflatten<5;iflatten++)
	{
		TString hadronV2HistogramName[2];
		// hadronV2HistogramName[0] = "hadron_"+mFlattenNames[i]+"_1";
		// hadronV2HistogramName[1] = "hadronsum_"+mFlattenNames[i]+"_1";
		hadronV2HistogramName[0] = Form("hadron_%s_%i",mFlattenNames[iflatten].Data(),mIndexEtaGap);
		hadronV2HistogramName[1] = Form("hadronsum_%s_%i",mFlattenNames[iflatten].Data(),mIndexEtaGap);
		for(int iname=0;iname<2;iname++)
			hadronV2Histogram[iflatten][iname] = (TH1D *)v2File->Get(hadronV2HistogramName[iname].Data())->Clone(hadronV2HistogramName[iname].Data());
		hadronV2Histogram[iflatten][0]->Divide(hadronV2Histogram[iflatten][1]);
	}
	hadronV2Histogram[1][0]->Multiply(hadronV2Histogram[3][0]);
	hadronV2Histogram[2][0]->Multiply(hadronV2Histogram[4][0]);
	hadronV2Histogram[0][0]->Add(hadronV2Histogram[1][0],-1);
	hadronV2Histogram[0][0]->Add(hadronV2Histogram[2][0],-1);
	float xbin[10] = {0,5,10,20,30,40,50,60,70,80};
	TH1D *hadronV2 = new TH1D("hadronV2","hadronV2;cent %;haron v_{2}",9,xbin);
  vector<double> hadronV2Values;
  vector<double> hadronV2Errors;
  mLog<<"centrality\t\t\t\thadron v2\t\t\t\thadron v2 stat.error"<<endl;
	for(int k=0;k<9;k++)
	{
		double mV2 = hadronV2Histogram[0][0]->GetBinContent(k+1);
		double mV2_err = hadronV2Histogram[0][0]->GetBinError(k+1);
		if(mV2>0)
			mV2 = sqrt(mV2);
		else 
			mV2 = 1;
		mV2_err = mV2_err/2/mV2;
		hadronV2->SetBinContent(9-k,mV2);
		hadronV2->SetBinError(9-k,mV2_err);
    hadronV2Values.push_back(mV2);
    hadronV2Errors.push_back(mV2_err);
    mLog<<"["<<xbin[8-k]<<","<<xbin[9-k]<<"]\t\t\t\t";
    mLog<<mV2<<"\t\t\t\t"<<mV2_err<<endl;
	}
  std::vector<std::vector<double> >  hadronV2Result;
  hadronV2Result.push_back(hadronV2Values);
  hadronV2Result.push_back(hadronV2Errors);
  return hadronV2Result;
}

pair<double,double> D0CorV2Plotter::fit_hist(TH1D *histo, TCanvas *cfg, int iptbin ,double nSigma,double fitArray[3])
{

	ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000); 
	cfg->cd(iptbin);
  pair<double,double> fitResult;

	cout << "/////////////////////////////////////********        i             **************         " << iptbin << endl;
	histo->SetMarkerSize(0.8);
	histo->SetLineColor(2);
	histo->SetMarkerColor(2);
	histo->SetMarkerStyle(20);
	histo->GetXaxis()->SetNdivisions(505);
	histo->GetXaxis()->SetTitle("m_{#piK} (GeV)");
	histo->GetYaxis()->SetTitle("Counts");

	double fit_range_low = 1.7;//eff_fit_range_low[iptbin];
	double fit_range_high = 2.0;//eff_fit_range_high[iptbin];

	histo->GetXaxis()->SetRangeUser(fit_range_low, fit_range_high);

	//.. fit with a Gaussian and pol
  TF1 *fit_fun;
  if(!mIsYieldSystematic)
    fit_fun = new TF1("fit_fun", "gausn(0) + pol2(3)", fit_range_low, fit_range_high);
  else 
    fit_fun = new TF1("fit_fun", "gausn(0) + expo(3)", fit_range_low, fit_range_high);
	float max = histo->GetMaximum();
	histo->SetMaximum(1.1 * max);

	float p0 = 1000, p1 = 1.87, p2 = 0.02;
	float p0_L = 0, p1_L = 1.84, p2_L = 0;
	float p0_H = 2*max, p1_H = 1.9, p2_H = 0.05;

	float p3 = -1. * max, p4 = max, p5 = -1. * max;

	int pass = 0;
	int fittingtry = 0;

	char sig_print[100], chi2_print[100], mean_print[100], sigma_print[100],sb_ratio[100];

	while (!pass) {

		fit_fun->SetParameter(0, p0);
		fit_fun->SetParameter(1, p1);
		fit_fun->SetParameter(2, p2);

		//.. fit constraint ..
		fit_fun->SetParLimits(0, p0_L, p0_H);
		fit_fun->SetParLimits(1, p1_L, p1_H);
		fit_fun->SetParLimits(2, p2_L, p2_H);

		//        fit_fun->SetParameter(3, p3);
		//		fit_fun->SetParameter(4, p4);
		//		fit_fun->SetParameter(5, p5);

		if( fittingtry == 0 )
			histo->Fit(fit_fun,"L","", fit_range_low, fit_range_high);
		else 
			histo->Fit(fit_fun,"L","", fit_range_low, fit_range_high);

		//.. draw foreground and background ..
		histo->Draw();

		TF1* fit_fun_1st = (TF1*)fit_fun->Clone("fit_fun_1st");
		fit_fun_1st->SetParameter(3, 0);
		fit_fun_1st->SetParameter(4, 0);
		fit_fun_1st->SetParameter(5, 0);
		//        fit_fun_1st->Draw("same");


		TF1* fit_fun_bg = (TF1*)fit_fun->Clone("fit_fun_bg");
		//        TF1* fit_fun_bg = new TF1("fit_fun_bg", fitfunction, cut_m_low, cut_m_high, 6);
		fit_fun_bg->SetParameter(0, 0);
		fit_fun_bg->SetParameter(1, 0);
		fit_fun_bg->SetParameter(2, 0);
		//		fit_fun_bg->SetParameter(3, fit_fun->GetParameter(3));
		//		fit_fun_bg->SetParameter(4, fit_fun->GetParameter(4));
		//		fit_fun_bg->SetParameter(5, fit_fun->GetParameter(5));


		fit_fun_bg->SetLineColor(8);
		fit_fun_bg->SetLineStyle(2);
		fit_fun_bg->Draw("same");


		fittingtry++;

		//    if( ptbins[iptbin] > lowrange && ptbins[iptbin+1] < highrange )
		{
			float binwidth = 0.01;//histo->GetBinWidth(10);
			//float ptbinwidth = ptbins[iptbin+1] - ptbins[iptbin];
			//counts->SetBinContent( iptbin+1, fit_fun->GetParameter(0)/( binwidth * ptbinwidth ));
			//counts->SetBinError( iptbin+1, fit_fun->GetParError(0)/( binwidth * ptbinwidth ));

			float Nsig = fit_fun->GetParameter(0)/( binwidth );
			float err_Nsig = fit_fun->GetParError(0)/( binwidth );
			float fitchi2 = fit_fun->GetChisquare();
			float fitmeanerror = fit_fun->GetParError(1);
			float fitsigmaerror = fit_fun->GetParError(2);
			int noffreepara = fit_fun->GetNumberFreeParameters();
			int noffitpoints = fit_fun->GetNumberFitPoints();

			float fitmean = fit_fun->GetParameter(1);
			float fitsigma = fit_fun->GetParameter(2);

			//hfg_masssigma->SetBinContent(iptbin+1, fitsigma);
			//hfg_masssigma->SetBinError(iptbin+1, fitsigmaerror);

			cout << " fitchi2: " << fitchi2 << "   noffreepara: " << noffreepara << "  noffitpoints: " << noffitpoints << endl;

			//      if( !isMC )
			//	sprintf( sig_print,"N_{sig}: %7.1f#pm%7.1f", Nsig, err_Nsig);
			//      else
			sprintf( sig_print,"N_{sig}: %7.2f#pm%7.2f", Nsig, err_Nsig);
			sprintf( chi2_print, "#chi^{2}#/d.o.f: %3.2f", fitchi2/( noffitpoints - noffreepara));
			sprintf( mean_print, "mean: %6.4f#pm%6.4f", fitmean,fitmeanerror);
			sprintf( sigma_print, "#sigma: %6.4f#pm%6.4f", fitsigma,fitsigmaerror);
			cout<<fitmean<<"\t"<<fitsigma<<endl;

			TF1 *g = new TF1("g","gausn(0)",1.6,2.1);
			g->SetParameters(Nsig*binwidth,fitmean,fitsigma);
			float inteSig = g->Integral(fitmean-nSigma*fitsigma,fitmean+nSigma*fitsigma) / (binwidth);
			float inteSig_err = err_Nsig * inteSig/Nsig;
      fitResult.first = inteSig;
      fitResult.second = inteSig_err;
			fitArray[0] = inteSig;
			fitArray[1] = fit_fun->Integral(fitmean-nSigma*fitsigma,fitmean+nSigma*fitsigma)/binwidth;
			fitArray[2] = fit_fun->IntegralError(fitmean-nSigma*fitsigma,fitmean+nSigma*fitsigma)/binwidth;
			//fout<<iptbin<<"fitmean = "<<fitmean<<endl;
			//fout<<iptbin<<"fitsigma = "<<fitsigma<<endl;
			double sbratio = fitArray[0]/fitArray[1];
			//sprintf( sb_ratio, "#sb raito: %6.4f", sbratio);


			if (fittingtry == 2)
			{
				TLatex Tl;
				Tl.SetNDC();
				Tl.SetTextAlign(12);
				Tl.SetTextSize(0.06);
				Tl.DrawLatex(0.15,0.8, sig_print);
				Tl.DrawLatex(0.15,0.7, chi2_print);
				Tl.DrawLatex(0.15,0.6, sb_ratio);
				Tl.DrawLatex(0.55,0.8, mean_print);
				Tl.DrawLatex(0.55,0.7, sigma_print);
			}

		}

		if (fittingtry == 2)  
		{
			pass = 1;

		}
		if(!pass) {
			p0 = fit_fun->GetParameter(0);
			p1 = fit_fun->GetParameter(1);
			p2 = fit_fun->GetParameter(2);
			//            p1_L = 1.84, p2_L = 0;
			//            p1_H = 1.9, p2_H = 0.05;
		}
	}
	double fitpar[6];
	for(int i=0;i<6;i++)
	{
		fitpar[i] = fit_fun->GetParameter(i);
	}

	return fitResult;
}
TH1D *D0CorV2Plotter::getSignalV2(TH1D *candOverSgn,TH1D *bkgV2,TH1D *candV2)
{
	TH1D *signalV2 = (TH1D *)candV2->Clone("signalV2");
	signalV2->Add(bkgV2,-1);
	signalV2->Multiply(candOverSgn);
  signalV2->Add(bkgV2,1);
  for(unsigned int ipt=1;ipt<7;ipt++)
  {
    double v_candOverSgn = candOverSgn->GetBinContent(ipt);
    double v_bkgV2 = bkgV2->GetBinContent(ipt);
    double v_candV2 = candV2->GetBinContent(ipt);
    double err_candOverSgn = candOverSgn->GetBinError(ipt);
    double err_bkgV2 = bkgV2->GetBinError(ipt);
    double err_candV2 = candV2->GetBinError(ipt);
    double v2CalMatrix[6] = {v_candOverSgn,v_bkgV2,v_candV2,err_candOverSgn,err_bkgV2,err_candV2};
    signalV2->SetBinError(ipt,getV2Error(v2CalMatrix));
  }
  return signalV2;
}

