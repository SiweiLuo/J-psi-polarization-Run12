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

#define NFILE 1000
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

TFile *efficiencyfile[NFILE];
TFile *rawdatafile[NFILE][NTRIG];
TFile *rawdatafile1[NFILE][NTRIG];

TCanvas *canvas[NTRIG][NPT][NFRAME];//marker

TH2F* mle[NTRIG][NPT][NFRAME];
TFile* histograms;
TFile* histograms_copy;
TFile* likelihoodhistograms[NTRIG][NPT][NFRAME];

float lambda[2][2]; // theta, phi; initial, terminal;

TF2* sigmafcn;

double sample[NTRIG][NPT][NFRAME];

void mlemethod(int ifile =0,int itrig=0, int xsect =11,int ysect=11,int rebin=1,int witheff=1){

	gRandom = new TRandom3();
	gRandom->SetSeed(0);

	histograms = new TFile("../histograms20161206.root","read"); //40*40 bining
	histograms_copy = new TFile("../templates20161206.root","read"); //40*40 bining

	correcteddata(ifile,rebin);

//	for(int ipt=0;ipt<6;ipt++){
	for(int ipt=3;ipt<4;ipt++){
//		for(int iframe=0;iframe<2;iframe++){
		for(int iframe=0;iframe<1;iframe++){
			maxlikelihood(ifile,itrig,ipt,iframe,xsect,ysect,rebin,witheff);	
		}
	}
}

void maxlikelihood(int ifile =0,int itrig =0, int ipt =0, int iframe=0,int xsect =11,int ysect =11,int rebin=0,int witheff =1){
	TH2F* templates[100][100];
	likelihoodhistograms[itrig][ipt][iframe] = new TFile(Form("~/polresults/20160707/MLEmethod/rootfiles/sys_%d/likelihood_%d_%d_%d_%d_%d_%d.root",ifile,ifile,itrig,ipt,iframe,xsect,ysect),"recreate");
	mle[itrig][ipt][iframe] = new TH2F(Form("mle_%d_%d_%d",itrig,ipt,iframe),"MLE;#lambda_{#theta};#lambda_{#phi}",CHIXBIN,-1,1,CHIYBIN,-1,1);
	mle[itrig][ipt][iframe]->Sumw2();

	TH2F* rawdata = rawdata2D[ifile][itrig][ipt][iframe]->Clone();
	rawdata->Sumw2();

	TH2F* eff = (TH2F*)eff2D[ifile][itrig][ipt][iframe][2]->Clone();
	eff->SetName(Form("eff_%d_%d_%d_%d",ifile,itrig,ipt,iframe));
	eff->Sumw2();
//	if(witheff==1)rawdata->Multiply(eff);
//	double rawdatacounts = rawdata->Integral();
//	rawdata->Scale(sample[itrig][ipt][iframe]/rawdatacounts);
//	rawdata->Scale(20*sample[itrig][ipt][iframe]/rawdatacounts);
	likelihoodhistograms[itrig][ipt][iframe]->cd();

	for(int xnbin=xsect;xnbin<xsect+10;xnbin++){
		for(int ynbin=ysect;ynbin<ysect+10;ynbin++){
			double lambda1,lambda2,likelihood=0,sum=0,likelihood2=0;
			lambda1 = mle[itrig][ipt][iframe]->GetXaxis()->GetBinCenter(xnbin);
			lambda2 = mle[itrig][ipt][iframe]->GetYaxis()->GetBinCenter(ynbin);

			templates[xnbin-1][ynbin-1] = (TH2F*)histograms->Get(Form("theta_%d_phi_%d",xnbin-1,ynbin-1));
			templates[xnbin-1][ynbin-1]->SetName(Form("template_%d_%d",xnbin-1,ynbin-1));
			templates[xnbin-1][ynbin-1]->Sumw2();
			if(rebin==1){
				templates[xnbin-1][ynbin-1]->RebinX(4);
				templates[xnbin-1][ynbin-1]->RebinY(4);
			}
			if(witheff==1)templates[xnbin-1][ynbin-1]->Multiply(eff);
			templates[xnbin-1][ynbin-1]->Scale(1./templates[xnbin-1][ynbin-1]->Integral());
			double scale = 0.001/templates[xnbin-1][ynbin-1]->Integral();
			cout<<"scale = "<<scale<<endl;
			//"   "<<1./templates[xnbin-1][ynbin-1]->Integral()<<" or "<<1/templates[xnbin-1][ynbin-1]->Integral()<<endl;
//			templates[xnbin-1][ynbin-1]->Scale(rawdata->Integral()/templates[xnbin-1][ynbin-1]->Integral());

			for(int xbin=1;xbin<XNBIN+1;xbin++){
				for(int ybin=1;ybin<YNBIN+1;ybin++){
					if(templates[xnbin-1][ynbin-1]->GetBinContent(xbin,ybin)<1e-6) continue;
					double costheta,phi;
					costheta = rawdata->GetYaxis()->GetBinCenter(ybin);
					phi = rawdata->GetXaxis()->GetBinCenter(xbin);
					likelihood += -1*rawdata->GetBinContent(xbin,ybin)*TMath::Log(templates[xnbin-1][ynbin-1]->GetBinContent(xbin,ybin));
					likelihood2 += -1*rawdata->GetBinContent(xbin,ybin)*TMath::Log(scale*templates[xnbin-1][ynbin-1]->GetBinContent(xbin,ybin));
				}
			}
			cout<<"xnbin = "<<xnbin<<"ynbin"<<ynbin<<" likelihood - likelihood2 ="<<likelihood-likelihood2<<endl;	
			mle[itrig][ipt][iframe]->Fill(lambda1,lambda2,likelihood);
			mle[itrig][ipt][iframe]->SetBinError(xnbin,ynbin,2);
		}
	}
	rawdata->Write("",TObject::kOverwrite);
	//	sigmafcn->Write();
	eff->Write("",TObject::kOverwrite);
	mle[itrig][ipt][iframe]->Write("",TObject::kOverwrite);
}

void correcteddata(int file = 0,int rebin = 0){

	efficiencyfile[file] = new TFile("../rootfile/OutFile_sys0.root","read"); 
	for(int trig=0;trig<NTRIG;trig++){
		//			if((file>=20 && file<=23) || (file>=29 && file<=30) || file==35) rawdatafile[file][trig] = new TFile(Form("rootfile/%s_sys%d.root",trigSet[trig].Data(),file),"read");//marker
		rawdatafile[file][trig] = new TFile(Form("../rootfile/%s_sys%d.root",trigSet[trig].Data(),0),"read");
		//			else rawdatafile[file][trig] = new TFile("../zbackupfiles/histograms20161025.root","read");//marker
		for(int pt=0;pt<NPT;pt++){
			for(int frame=0;frame<NFRAME;frame++){

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

				int max = (((TH3F*)efficiencyfile[file]->Get("hHT0JpsiCosThetaPhiPt1"))->GetZaxis())->FindBin(PtEdge[pt+1]-0.001);	
				int min = (((TH3F*)efficiencyfile[file]->Get("hHT0JpsiCosThetaPhiPt1"))->GetZaxis())->FindBin(PtEdge[pt]+0.001);	
				if(frame==0){
					eff3D[file][trig][pt][frame][0] = (TH3F*)efficiencyfile[file]->Get(Form("h%sJpsiCosThetaPhiPt1",trigName[trig].Data()));
					eff3D[file][trig][pt][frame][0]->SetName(Form("eff3D_pass_%d_%d_%d_%d",file,trig,pt,frame));
					eff3D[file][trig][pt][frame][1] = (TH3F*)efficiencyfile[file]->Get("hJpsiCosThetaPhiPt1");
					eff3D[file][trig][pt][frame][1]->SetName(Form("hJpsiCosThetaPhiPt1_%s",trigName[trig].Data()));
				}
				else{	
					eff3D[file][trig][pt][frame][0] = (TH3F*)efficiencyfile[file]->Get(Form("h%sJpsiCosThetaPhiPtCS1",trigName[trig].Data()));
					eff3D[file][trig][pt][frame][0]->SetName(Form("eff3D_pass_%d_%d_%d_%d",file,trig,pt,frame));
					eff3D[file][trig][pt][frame][1] = (TH3F*)efficiencyfile[file]->Get("hJpsiCosThetaPhiPtCS1");
					eff3D[file][trig][pt][frame][1]->SetName(Form("hJpsiCosThetaPhiPtCS1_%s",trigName[trig].Data()));
				}
				eff3D[file][trig][pt][frame][0]->Sumw2();
				eff3D[file][trig][pt][frame][1]->Sumw2();

				rawdata3D[file][trig][pt][frame][2]->GetZaxis()->SetRange(min,max);


				eff3D[file][trig][pt][frame][0]->GetZaxis()->SetRange(min,max);
				eff3D[file][trig][pt][frame][1]->GetZaxis()->SetRange(min,max);
				eff3D[file][trig][pt][frame][0]->SetName(Form("%d_%d_%d_%d_0",file,trig,pt,frame));
				eff3D[file][trig][pt][frame][0]->Sumw2();
				eff3D[file][trig][pt][frame][1]->SetName(Form("%d_%d_%d_%d_1",file,trig,pt,frame));
				eff3D[file][trig][pt][frame][1]->Sumw2();

				rawdata2D[file][trig][pt][frame] = (TH2F*)rawdata3D[file][trig][pt][frame][2]->Project3D("xy");
				//					rawdata2D[file][trig][pt][frame]->SetName(Form("rawdata_%d_%s_pT%d_frame%d",file,trigName[trig].Data(),pt,frame));
				rawdata2D[file][trig][pt][frame]->SetName(Form("raw_%d_%d_%d_%d",file,trig,pt,frame));
				rawdata2D[file][trig][pt][frame]->Sumw2();
				if(rebin==1){
					rawdata2D[file][trig][pt][frame]->RebinX(NREBIN);//10*10
					rawdata2D[file][trig][pt][frame]->RebinY(NREBIN);//10*10
				}

				//					double sample = rawdata2D[file][trig][pt][frame]->GetEntries();
				sample[trig][pt][frame] = rawdata2D[file][trig][pt][frame]->Integral();
				cout<<"sample = "<<sample[trig][pt][frame]<<endl;
				double phi,costheta;			
				rawdatafile1[file][trig] = new TFile("../zbackupfiles/histograms20161025.root","read");//marker

				//					rawdata2D[file][trig][pt][frame] = (TH2F*)rawdatafile1[file][trig]->Get("theta_50_phi_50");
				sigmafcn = (TF2*)rawdatafile1[file][trig]->Get("sigma_50_50");
				sigmafcn->SetName(Form("sigma_%d_%d_%d_%d",file,trig,pt,frame));
				cout<<"sigmafcn = "<<sigmafcn<<endl;

				TFile* rawmcfile = new TFile(Form("~/polresults/20160707/MLEmethod/rootfiles/sys_%d/raw_%d_%d_%d_%d.root",file,file,trig,pt,frame),"read");
				if(rawmcfile==0x0) rawmcfile = new TFile(Form("~/polresults/20160707/MLEmethod/rootfiles/sys_%d/raw_%d_%d_%d_%d.root",file,file,trig,pt,frame),"create");

				cout<<"rawmcfile="<<rawmcfile<<endl;

				Bool_t rawMCexists;
				rawMCexists=(rawmcfile->GetListOfKeys())->Contains(Form("raw_%d_%d_%d_%d",file,trig,pt,frame));
				cout<<"rawMCexists = "<<rawMCexists<<endl;
				if(rawMCexists!=1){
					rawmcfile->cd();
					rawdata2D[file][trig][pt][frame] = new TH2F(Form("raw_%d_%d_%d_%d",file,trig,pt,frame),Form("raw_%d_%d_%d_%d",file,trig,pt,frame),40,-TMath::Pi(),TMath::Pi(),40,-1,1);
					for(int i=0;i<sample[trig][pt][frame];i++){	
						sigmafcn->GetRandom2(phi,costheta);	
						rawdata2D[file][trig][pt][frame]->Fill(phi,costheta);		
					}

					rawdata2D[file][trig][pt][frame]->SetName(Form("raw_%d_%d_%d_%d",file,trig,pt,frame));
					rawdata2D[file][trig][pt][frame]->Sumw2();
					//					rawdata2D[file][trig][pt][frame]->Scale(sample);

					if(rebin==1){
						rawdata2D[file][trig][pt][frame]->RebinX(NREBIN);//10*10
						rawdata2D[file][trig][pt][frame]->RebinY(NREBIN);//10*10
					}
					cout<<"writing the mcraw data"<<endl;
					rawdata2D[file][trig][pt][frame]->Write();	
					rawmcfile->Close();
				}
				else {
					rawmcfile = new TFile(Form("~/polresults/20160707/MLEmethod/rootfiles/sys_%d/raw_%d_%d_%d_%d.root",file,file,trig,pt,frame),"read");
					cout<<"reading the mcraw data"<<endl;
					rawdata2D[file][trig][pt][frame]=(TH2F*)rawmcfile->Get(Form("raw_%d_%d_%d_%d",file,trig,pt,frame));
				}


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
			}
		}
	}
}

