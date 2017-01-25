#include <algorithm>
#include <iostream>
#include "TCanvas.h"
#include "TError.h"
#include "TF2.h"
#include "TFile.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

#define NFILE 100000
#define NTRIG 3
#define NPT 6
#define NPHASE 2 
#define NFRAME 2
#define NREBIN 4
#define XNBIN 40
#define YNBIN 40
#define CHIXBIN 100
#define CHIYBIN 100
//#define rebin 0 //1

Int_t PtEdge[NPT+1] = {0,2,3,4,6,8,14};

TString trigName[NTRIG] = {"HT0","HT1","HT2"};
TString trigSet[NTRIG] = {"ht0","ht1","ht2"};

TH2F *rawdata2D[NFILE][NTRIG][NPT][NFRAME];
TH3F *rawdata3D[NFILE][NTRIG][NPT][NFRAME][3]; // unlike, like , unlike-like
TH2F *eff2D[NFILE][NTRIG][NPT][NFRAME][3]; // pass , total , ratio;
TH3F *eff3D[NFILE][NTRIG][NPT][NFRAME][2]; // pass , total ;
TH2F *chi2result[NFILE][NTRIG][NPT][NFRAME];

TFile *efficiencyfile;
TFile *rawdatafile[NFILE][NTRIG];
TFile *rawdatafile1[NFILE][NTRIG];

TCanvas *canvas[NTRIG][NPT][NFRAME];//marker

TH2F* mle[NTRIG][NPT][NFRAME];
TFile* histograms;
TFile* histograms_copy;
TFile* likelihoodhistograms[NTRIG][NPT][NFRAME];

float lambda[2][2]; // theta, phi; initial, terminal;

TF2* sigmafcn;

void rawMCdata(int file,int trig,int pt,int frame){

	gRandom = new TRandom3();
	gRandom->SetSeed(0);

	int rebin=1;
	rawdatafile[file][trig] = new TFile(Form("../rootfile/%s_sys%d.root",trigSet[trig].Data(),0),"read");
	if(frame==0){
		rawdata3D[file][trig][pt][frame][0] = (TH3F*)(rawdatafile[file][trig]->Get("hJpsiCosThetaPhiPt"));
		rawdata3D[file][trig][pt][frame][1] = (TH3F*)(rawdatafile[file][trig]->Get("hJpsiCosThetaPhiPtBG"));
		rawdata3D[file][trig][pt][frame][0]->SetName("rawdata3D_unlike_hx");
		rawdata3D[file][trig][pt][frame][1]->SetName("rawdata3D_like_hx");
	}
	else{
		rawdata3D[file][trig][pt][frame][0] = (TH3F*)(rawdatafile[file][trig]->Get("hJpsiCosThetaPhiPtCS"));
		rawdata3D[file][trig][pt][frame][1] = (TH3F*)(rawdatafile[file][trig]->Get("hJpsiCosThetaPhiPtCSBG"));
		rawdata3D[file][trig][pt][frame][0]->SetName("rawdata3D_unlike_cs");
		rawdata3D[file][trig][pt][frame][1]->SetName("rawdata3D_like_cs");
	}
	rawdata3D[file][trig][pt][frame][0]->Sumw2();
	rawdata3D[file][trig][pt][frame][1]->Sumw2();

	rawdata3D[file][trig][pt][frame][2] = (TH3F*)rawdata3D[file][trig][pt][frame][0]->Clone("rawdata3D_unlike_like");
	if(frame==0) rawdata3D[file][trig][pt][frame][2]->SetName("rawdata3D_unlike_like_hx");
	else rawdata3D[file][trig][pt][frame][2]->SetName("rawdata3D_unlike_like_cs");
	rawdata3D[file][trig][pt][frame][2]->Sumw2();
	rawdata3D[file][trig][pt][frame][2]->Add(rawdata3D[file][trig][pt][frame][1],-1);


	int max = (rawdata3D[file][trig][pt][frame][2]->GetZaxis())->FindBin(PtEdge[pt+1]-0.001);	
	int min = (rawdata3D[file][trig][pt][frame][2]->GetZaxis())->FindBin(PtEdge[pt]+0.001);	
		


	rawdata3D[file][trig][pt][frame][2]->GetZaxis()->SetRange(min,max);


	rawdata2D[file][trig][pt][frame] = (TH2F*)rawdata3D[file][trig][pt][frame][2]->Project3D("xy");
	//					rawdata2D[file][trig][pt][frame]->SetName(Form("rawdata_%d_%s_pT%d_frame%d",file,trigName[trig].Data(),pt,frame));
	rawdata2D[file][trig][pt][frame]->SetName(Form("raw_%d_%d_%d_%d",file,trig,pt,frame));
	rawdata2D[file][trig][pt][frame]->Sumw2();
	if(rebin==1){
		rawdata2D[file][trig][pt][frame]->RebinX(NREBIN);//10*10
		rawdata2D[file][trig][pt][frame]->RebinY(NREBIN);//10*10
	}

	double sample = rawdata2D[file][trig][pt][frame]->Integral();
	cout<<"sample = "<<sample<<endl;

	efficiencyfile = new TFile("../rootfile/OutFile_sys0.root","read"); 
		if(frame==0){
					eff3D[file][trig][pt][frame][0] = (TH3F*)efficiencyfile->Get(Form("h%sJpsiCosThetaPhiPt1",trigName[trig].Data()));
					eff3D[file][trig][pt][frame][0]->SetName(Form("eff3D_pass_%d_%d_%d_%d",file,trig,pt,frame));
					eff3D[file][trig][pt][frame][1] = (TH3F*)efficiencyfile->Get("hJpsiCosThetaPhiPt1");
					eff3D[file][trig][pt][frame][1]->SetName(Form("hJpsiCosThetaPhiPt1_%s",trigName[trig].Data()));
				}
				else{	
					eff3D[file][trig][pt][frame][0] = (TH3F*)efficiencyfile->Get(Form("h%sJpsiCosThetaPhiPtCS1",trigName[trig].Data()));
					eff3D[file][trig][pt][frame][0]->SetName(Form("eff3D_pass_%d_%d_%d_%d",file,trig,pt,frame));
					eff3D[file][trig][pt][frame][1] = (TH3F*)efficiencyfile->Get("hJpsiCosThetaPhiPtCS1");
					eff3D[file][trig][pt][frame][1]->SetName(Form("hJpsiCosThetaPhiPtCS1_%s",trigName[trig].Data()));
				}
				eff3D[file][trig][pt][frame][0]->Sumw2();
				eff3D[file][trig][pt][frame][1]->Sumw2();

				eff3D[file][trig][pt][frame][0]->GetZaxis()->SetRange(min,max);
				eff3D[file][trig][pt][frame][1]->GetZaxis()->SetRange(min,max);
				eff3D[file][trig][pt][frame][0]->SetName(Form("%d_%d_%d_%d_0",file,trig,pt,frame));
				eff3D[file][trig][pt][frame][0]->Sumw2();
				eff3D[file][trig][pt][frame][1]->SetName(Form("%d_%d_%d_%d_1",file,trig,pt,frame));
				eff3D[file][trig][pt][frame][1]->Sumw2();

				eff2D[file][trig][pt][frame][0] = (TH2F*)eff3D[file][trig][pt][frame][0]->Project3D("xy");
				eff2D[file][trig][pt][frame][0]->SetName(Form("eff_pass_%d_%d_%d_%d",file,trig,pt,frame));
				eff2D[file][trig][pt][frame][0]->Sumw2();
				if(rebin==1){
					eff2D[file][trig][pt][frame][0]->RebinX(NREBIN);//10*10
					eff2D[file][trig][pt][frame][0]->RebinY(NREBIN);//10*10
				}
				eff2D[file][trig][pt][frame][1] = (TH2F*)eff3D[file][trig][pt][frame][1]->Project3D("xy");
				eff2D[file][trig][pt][frame][1]->SetName(Form("eff_total_%d_%d_%d_%d",file,trig,pt,frame));
				eff2D[file][trig][pt][frame][1]->Sumw2();
				if(rebin==1){					
					eff2D[file][trig][pt][frame][1]->RebinX(NREBIN);//10*10
					eff2D[file][trig][pt][frame][1]->RebinY(NREBIN);//10*10
				}
				eff2D[file][trig][pt][frame][2] = (TH2F*)eff2D[file][trig][pt][frame][0]->Clone();
				eff2D[file][trig][pt][frame][2]->SetName(Form("eff_ratio_%d_%d_%d_%d",file,trig,pt,frame));
				eff2D[file][trig][pt][frame][2]->Sumw2();

				eff2D[file][trig][pt][frame][2]->Divide(eff2D[file][trig][pt][frame][1]);







	rawdatafile1[file][trig] = new TFile("../zbackupfiles/histograms20161025.root","read");//marker
	sigmafcn = (TF2*)rawdatafile1[file][trig]->Get("sigma_50_50");
	sigmafcn->SetName(Form("sigma_%d_%d_%d_%d",file,trig,pt,frame));

	double phi,costheta;			

	TFile* rawmcfile = new TFile(Form("~/polresults/20160707/MLEmethod/rootfiles/sys_%d/raw_%d_%d_%d_%d.root",file,file,trig,pt,frame),"recreate");
	rawmcfile->cd();

	rawdata2D[file][trig][pt][frame] = new TH2F(Form("raw_%d_%d_%d_%d",file,trig,pt,frame),Form("raw_%d_%d_%d_%d",file,trig,pt,frame),40,-TMath::Pi(),TMath::Pi(),40,-1,1);
	for(int i=0;i<sample;){	
		sigmafcn->GetRandom2(phi,costheta);	
		int xbin = eff2D[file][trig][pt][frame][2]->GetXaxis()->FindBin(phi);
		int ybin = eff2D[file][trig][pt][frame][2]->GetYaxis()->FindBin(costheta);
		if(gRandom->Uniform()<eff2D[file][trig][pt][frame][2]->GetBinContent(xbin,ybin)){
			rawdata2D[file][trig][pt][frame]->Fill(phi,costheta);		
			i++;
		}
	}
	
//	rawdata2D[file][trig][pt][frame]->Scale(sample/10000);
	rawdata2D[file][trig][pt][frame]->SetName(Form("raw_%d_%d_%d_%d",file,trig,pt,frame));
	rawdata2D[file][trig][pt][frame]->Sumw2();
	//					rawdata2D[file][trig][pt][frame]->Scale(sample);

	if(rebin==1){
		rawdata2D[file][trig][pt][frame]->RebinX(NREBIN);//10*10
		rawdata2D[file][trig][pt][frame]->RebinY(NREBIN);//10*10
	}
	rawdata2D[file][trig][pt][frame]->Write();	
}
