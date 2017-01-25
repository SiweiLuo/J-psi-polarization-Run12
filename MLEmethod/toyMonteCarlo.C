#include "TH2F.h"
#include "TCanvas.h"
#include "TF2.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TStyle.h"
#include "TRandom2.h"
#include "TError.h"
#include "TMinuit.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <iostream>
#include <algorithm>

#define NFILE 1000
#define NTRIG 3
#define NPT 6
#define NPHASE 2 
#define NFRAME 2
#define NPLOT 6
#define NREBIN 4
#define XNBIN 10
#define YNBIN 10
#define CHIXBIN 100
#define CHIYBIN 100
#define pi 3.1415926

Int_t PtEdge[NPT+1] = {0,2,3,4,6,8,14};

TString trigSet[NTRIG] = {"ht0","ht1","ht2"};
TString trigName[NTRIG] = {"HT0","HT1","HT2"};
TString frameName[NFRAME] = {"HX","CS"};
TString phaseName[NPHASE+1] = {"#lambda_{#theta}","#lambda_{#phi}","#lambda_{inv}"};

TFile* lambdas;

TH2F *rawdata2D[NFILE][NTRIG][NPT][NFRAME][2];// raw data, corrected data  // marker remove the last dimension
TH3F *rawdata3D[NFILE][NTRIG][NPT][NFRAME][3]; // unlike, like , unlike-like
TH2F *eff2D[NFILE][NTRIG][NPT][NFRAME][3]; // pass , total , ratio;
TH1F *costheta_eff[NFILE][NTRIG][NPT][NFRAME][2]; // eff2D project to costheta pass , total 
TH1F *phi_eff[NFILE][NTRIG][NPT][NFRAME][2]; // eff2D project to phi pass , total 
TH3F *eff3D[NFILE][NTRIG][NPT][NFRAME][3]; // pass , total , ratio;
TH2F *chi2result[NFILE][NTRIG][NPT][NFRAME];

TH1F *costheta[NFILE][NTRIG][NPT][NFRAME][2]; // 1D raw data, 1D corrected data
TH1F *phi[NFILE][NTRIG][NPT][NFRAME][2]; // 1D raw data, 1D corrected data
TGraphAsymmErrors* efficiency1D[NFILE][NTRIG][NPT][NFRAME][NPHASE];
TH1F * eff1D[NFILE][NTRIG][NPT][NFRAME][NPHASE];

TGraphAsymmErrors *efficiency;
TFile *efficiencyfile[NFILE];
TFile *rawdatafile[NFILE][NTRIG];

TF1 *fit1[NPT][NFRAME][2]; //npt frame default/revised
TF1 *fit2[NPT][NFRAME][2];

double lambda_theta[NFILE][NTRIG][NPT][NFRAME];//frame
double lambda_theta_err[NFILE][NTRIG][NPT][NFRAME];//frame
double lambda_phi[NFILE][NTRIG][NPT][NFRAME];//frame
double lambda_phi_err[NFILE][NTRIG][NPT][NFRAME];//frame
double lambda_invariant[NFILE][NTRIG][NPT][NFRAME][3];// central value , statistic error , systematic error

double fit2dtheta[NFILE][NTRIG][NPT][NFRAME];
double fit2dtheta_err[NFILE][NTRIG][NPT][NFRAME];
double fit2dphi[NFILE][NTRIG][NPT][NFRAME];
double fit2dphi_err[NFILE][NTRIG][NPT][NFRAME];
double fitparameters[NFILE][NTRIG][NPT][NFRAME][4];//theta , theta_err , phi , phi_err
double ptsmearing[NTRIG][NPT][NFRAME][NPHASE],weight[NTRIG][NPT][NFRAME][NPHASE],nhitsfit[NTRIG][NPT][NFRAME][NPHASE],nsigma[NTRIG][NPT][NFRAME][NPHASE],dsmadc[NTRIG][NPT][NFRAME][NPHASE],beta[NTRIG][NPT][NFRAME][NPHASE],poe[NTRIG][NPT][NFRAME][NPHASE],pol[NTRIG][NPT][NFRAME][NPHASE];
double systematic_error[NTRIG][NPT][NFRAME][NPHASE+1]; // theta, phi, inv.

Double_t matrix[4][4];
Double_t covariant[NTRIG][NPT][NFRAME];

TCanvas *canvas[NFILE][NTRIG][NPT][NFRAME];//marker
TCanvas *systematic[NTRIG];
TCanvas *check;

double fParamVal[NFILE][NTRIG][NPT][NFRAME][4];
double fParamErr[NFILE][NTRIG][NPT][NFRAME][4];
double arglist[10];
Int_t ierflg=0;
Int_t fitflag[4][NPT][2][4];
TH1F* fitdata[NTRIG][NPT][NFRAME][NPHASE];// theta , phi

TGraphErrors* lambda_parameters[NTRIG][NFRAME][NPHASE+1];// theta , phi , invariant
TGraphErrors* lambda_parameters_sys[NTRIG][NFRAME][NPHASE+1];// theta , phi , invariant

TLegend* legend;

TFile *outputfile;
TGraphErrors* comparisonlambda[NTRIG][NPLOT];

TH2F* chi2;
TH2F* likelihood_contour;

double lambda[2][2];
TF2* mlefcn;
TF2* contfcn;

//double contourlevel[3][6]={
//	{10,10,10,10,10,10},{10,10,10,10,10,10},{10,10,10,10,10,5}
//};

double contourlevel[2][3][6]={
	{{10,10,10,10,10,10},{10,10,10,10,10,10},{10,10,10,10,10,5}},
	{{10,10,10,10,10,10},{10,10,10,10,10,10},{10,10,10,10,10,5}}
};

TFile* inputfile[1000][NTRIG][NPT][NFRAME];
TGraphErrors* getlambda[1000][NTRIG][NPT][NFRAME]; 	
TH2F* toyMC[NTRIG][NPT][NFRAME];
TH1F* toyMC1D[NTRIG][NPT][NFRAME];
TH1F* lambda_theta_uncertainty[NTRIG][NPT][NFRAME];
TH1F* lambda_phi_uncertainty[NTRIG][NPT][NFRAME];

void toyMonteCarlo(int filemin,int filemax,int trig,int pt,int frame){
//void toyMonteCarlo(int filemin,int filemax){

	TH1::AddDirectory(kFALSE);
	TFile* toyMCfile = new TFile(Form("~/polresults/20160707/MLEmethod/combinedrootfiles/toyMonteCarlo/toyMCfile_%d_%d.root",filemin,filemax),"update");
	double x,y;
	
//	for(int trig=0;trig<NTRIG;trig++){
//		for(int pt=0;pt<NPT;pt++){
//			for(int frame=0;frame<NFRAME;frame++){
				toyMC[trig][pt][frame] = new TH2F(Form("toyMC_%d_%d_%d",trig,pt,frame),Form("N times pseudo experiments;#lambda_{#theta};#lambda_{#phi}",trig,pt,frame),100,-1,1,100,-1,1);
				toyMC1D[trig][pt][frame] = new TH1F(Form("toyMC1D_%d_%d_%d",trig,pt,frame),Form("N times pseudo experiments;#lambda_{#theta};#lambda_{#phi}",trig,pt,frame),100,-1,1);
				lambda_theta_uncertainty[trig][pt][frame] = new TH1F(Form("lambda_theta_%d_%d_%d",trig,pt,frame),Form("lambda_theta_%d_%d_%d",trig,pt,frame),100,0,1);
				lambda_phi_uncertainty[trig][pt][frame] = new TH1F(Form("lambda_phi_%d_%d_%d",trig,pt,frame),Form("lambda_phi_%d_%d_%d",trig,pt,frame),100,0,1);

				for(int file=filemin;file<=filemax;file++){
//					inputfile[file][trig][pt][frame] = new TFile(Form("~/polresults/20160707/MLEmethod/combinedrootfiles/sys_%d/lambda_trg%d_pt%d_frame%d.root",file,trig,pt,frame),"read");
//					TFile* openfile = new TFile(Form("~/polresults/20160707/MLEmethod/combinedrootfiles/sys_%d/lambda_trg%d_pt%d_frame%d.root",file,trig,pt,frame),"read");
					TString filename = Form("~/polresults/20160707/MLEmethod/combinedrootfiles/sys_%d/lambda_trg%d_pt%d_frame%d.root",file,trig,pt,frame);
					TFile* openfile = TFile::Open(filename,"read");
					if(openfile==0x0) {
						cout<<filename<<" don't exist"<<endl;
						continue;
					}
//					getlambda[file][trig][pt][frame] = 	(TGraphErrors*)inputfile[file][trig][pt][frame]->Get("lambda");
					getlambda[file][trig][pt][frame] = 	(TGraphErrors*)openfile->Get("lambda");
					if(getlambda[file][trig][pt][frame]==0x0) continue;
					getlambda[file][trig][pt][frame]->GetPoint(0,x,y);
					cout<<"filename = "<<filename<<endl;
					cout<<"x = "<<x<<"+/-"<<getlambda[file][trig][pt][frame]->GetErrorX(0)<<"; y = "<<y<<"+/-"<<getlambda[file][trig][pt][frame]->GetErrorY(0)<<endl;
					toyMC[trig][pt][frame]->Fill(x,y);
					toyMC1D[trig][pt][frame]->Fill(x);
					lambda_theta_uncertainty[trig][pt][frame]->Fill(getlambda[file][trig][pt][frame]->GetErrorX(0));
					lambda_phi_uncertainty[trig][pt][frame]->Fill(getlambda[file][trig][pt][frame]->GetErrorY(0));
//					inputfile[file][trig][pt][frame]->Close();
//					inputfile[file][trig][pt][frame]->Clear();
//					delete inputfile[file][trig][pt][frame];
					delete openfile;
				}
				toyMCfile->cd();
				toyMC[trig][pt][frame]->Write("",TObject::kOverwrite);
				toyMC1D[trig][pt][frame]->Write("",TObject::kOverwrite);
				lambda_theta_uncertainty[trig][pt][frame]->Write("",TObject::kOverwrite);
				lambda_phi_uncertainty[trig][pt][frame]->Write("",TObject::kOverwrite);
//			}
//		}
//	}
	toyMCfile->Close();
}
