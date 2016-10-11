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

#define NFILE 40
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

Int_t gfile,gtrig,gpt,gframe;

Float_t z1[40/NREBIN],z2[40/NREBIN],errorz1[40/NREBIN],errorz2[40/NREBIN],x[40/NREBIN],y[40/NREBIN];
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
double systematic_error[NTRIG][NPT][NFRAME][NPHASE];

Double_t matrix[4][4];
Double_t covariant[NTRIG][NPT][NFRAME];

TCanvas *canvas[NTRIG][NPT][NFRAME];//marker
TCanvas *systematic[NTRIG];
TCanvas *check;

double fParamVal[NFILE][NTRIG][NPT][NFRAME][4];
double fParamErr[NFILE][NTRIG][NPT][NFRAME][4];
double arglist[10];
Int_t ierflg=0;
Int_t fitflag[4][NPT][2][4];
TH1F* fitdata[NTRIG][NPT][NFRAME][NPHASE];// theta , phi

TGraphErrors* lambda_theta_hx;// marker write them as an array
TGraphErrors* lambda_phi_hx; 
TGraphErrors* lambda_theta_cs;
TGraphErrors* lambda_phi_cs; 
TGraphErrors* lambda_parameters[NTRIG][NFRAME][NPHASE+1];// theta , phi , invariant
TGraphErrors* lambda_parameters_sys[NTRIG][NFRAME][NPHASE+1];// theta , phi , invariant

TLegend* legend;

TFile *outputfile;
TGraphErrors* comparisonlambda[NTRIG][NPLOT];

TH2F* chi2;
TH2F* result_mean;
TH2F* result_sigma;

void lambdas(){
	gStyle->SetOptStat(0);
	outputfile = new TFile("~/jpsi/test20160210_Barbara/outputfile.root","read");//marker can be removed after crosschecking	
	int file=0; //marker for loop 
	for(int trig=0;trig<NTRIG;trig++){
		for(int pt=0;pt<NPT;pt++){
			for(int frame=0;frame<NFRAME;frame++){
				scan(file,trig,pt,frame);
			}
		}	
		lambdaparameters(file,trig);
	}
}


void scan(int file = 0, int trig = 0, int pt = 0, int frame = 0){
	TFile* inputfile = new TFile(Form("~/polresultspdsf/20160707/rootcombined/lambda_file%d_trg%d_pt%d_frame%d.root",file,trig,pt,frame),"read");

	TH2F* chi2 = (TH2F*)inputfile->Get(Form("chi2_%d_%d_%d_%d",file,trig,pt,frame));

	if(chi2==0x0 || chi2->GetMean()!=chi2->GetMean() || chi2->GetMean()==0 ){
		lambda_theta[file][trig][pt][frame] = -100;
		lambda_theta_err[file][trig][pt][frame] = 0.;
		lambda_phi[file][trig][pt][frame] = -100;
		lambda_phi_err[file][trig][pt][frame] = 0.;	
		return;
	}

	result_sigma = new TH2F("result_sigma","#sigma_{#chi^{2}};#lambda_{#theta};#lambda_{#phi}",CHIXBIN,-1,1,CHIYBIN,-1,1);
	std::vector<double> lambda_theta_error,lambda_phi_error;

	int xx,yy,zz;
	chi2->GetMinimumBin(xx,yy,zz);
	double minimum = chi2->GetBinContent(xx,yy);
	lambda_theta[file][trig][pt][frame] = chi2->GetXaxis()->GetBinCenter(xx);
	lambda_phi[file][trig][pt][frame] = chi2->GetYaxis()->GetBinCenter(yy);

	for(int i=1;i<CHIXBIN+1;i++){
		for(int j=1;j<CHIYBIN+1;j++){
			if(chi2->GetBinContent(i,j)<=minimum+1 && chi2->GetBinContent(i,j)>=minimum){
				result_sigma->Fill(chi2->GetXaxis()->GetBinCenter(i),chi2->GetYaxis()->GetBinCenter(j),chi2->GetBinContent(i,j));
				lambda_theta_error.push_back(chi2->GetXaxis()->GetBinCenter(i));
				lambda_phi_error.push_back(chi2->GetYaxis()->GetBinCenter(j));
			}
		}
	}
	result_sigma->ProjectionX()->Fit("gaus","Q");
	lambda_theta_err[file][trig][pt][frame] = 2*gaus->GetParameter(2);
	result_sigma->ProjectionY()->Fit("gaus","Q");
	lambda_phi_err[file][trig][pt][frame] = 2*gaus->GetParameter(2);
	cout<<"lambda_theta="<<lambda_theta[file][trig][pt][frame]<<"          "<<"lambda_phi="<<lambda_phi[file][trig][pt][frame]<<endl;
}

void lambdaparameters(int sys3,int trig){// plot and write lambda parameters 

	int selectedfile,file;
	if(sys3>=0 && sys3<=NFILE) file=sys3,selectedfile=sys3+1;
	else if(sys3==100) {
		file=0,selectedfile=NFILE;
	}
	else return;

	double x[NPT],y[NPT],x_err[NPT],y_err[NPT],y_sys[NPT],theta[NPT],thetaerr[NPT],phi[NPT],phierr[NPT];

	for(;file<selectedfile;file++){
		for(int frame=0;frame<NFRAME;frame++){ 
			for(int phase=0;phase<NPHASE+1;phase++){
				for(int pt=0;pt<NPT;pt++) {
					x[pt] = (PtEdge[pt+1]+PtEdge[pt])/2.+0.2*trig;
					x_err[pt] = (PtEdge[pt+1]-PtEdge[pt])/2.;
					theta[pt] = lambda_theta[file][trig][pt][frame];
					thetaerr[pt] = lambda_theta_err[file][trig][pt][frame];
					phi[pt] = lambda_phi[file][trig][pt][frame];
					phierr[pt] = lambda_phi_err[file][trig][pt][frame];

					if(phase==0){
						y[pt] = lambda_theta[file][trig][pt][frame];
						y_err[pt] = lambda_theta_err[file][trig][pt][frame];
					}
					if(phase==1){
						y[pt] = lambda_phi[file][trig][pt][frame];
						y_err[pt] = lambda_phi_err[file][trig][pt][frame];
					}
					if(phase==2){//calculate lambda_invariant and its statistic error
						y[pt] = (theta[pt]+3*phi[pt])/(1-phi[pt]);
						y_err[pt] = TMath::Sqrt(theta[pt]*theta[pt]/((1-phi[pt])*(1-phi[pt]))+phi[pt]*phi[pt]*TMath::Power((3+theta[pt]/((1-phi[pt])*(1-phi[pt]))),2)+2*(3+theta[pt])/TMath::Power((1-phi[pt]),3)*covariant[trig][pt][frame]);
						//							y_sys[pt] = y[pt]-(fitparameters[0][trig][pt][frame][0]+3*fitparameters[0][trig][pt][frame][2])/(1-fitparameters[0][trig][pt][frame][2]);
						//cout<<"invariant ==================="<<y[pt]<<"trig"<<trig<<"   "<<"frame"<<frame<<endl;
					}
				}
				lambda_parameters[trig][frame][phase] = new TGraphErrors(NPT,x,y,x_err,y_err);// marker SetName
				lambda_parameters[trig][frame][phase]->SetName(Form("lambdas_%d_%d_%d",trig,frame,phase));
				lambda_parameters_sys[trig][frame][phase] = new TGraphErrors(NPT,x,y,x_err,y_sys);
				lambda_parameters_sys[trig][frame][phase]->SetName(Form("lambdas_sys_%d_%d_%d",trig,frame,phase));
				if(sys3==100) {
					//			lambda_parameters[trig][frame][phase]->Write();		
					//			lambda_parameters_sys[trig][frame][phase]->Write();		
				}
			}
		}
	}
	drawlambdas(0,trig);
}

void drawlambdas(int drawoptions = 0,int trig){
	TCanvas* lambdacanvas; 
	lambdas = new TFile(Form("~/polresultspdsf/20160707/lambdas/lambdas_trig%d.root",trig),"recreate");
	lambdas->cd();
	if(drawoptions==1){
		lambdacanvas = new TCanvas("lambdacanvas","lambdacanvas",1200,800);
		lambdacanvas->Divide(3,2);
		for(int phase=0;phase<NPHASE+1;phase++){
			for(int frame=0;frame<NFRAME;frame++){ 
				lambdacanvas->cd(frame*3+phase+1);
				legend = new TLegend(0.7,0.7,0.89,0.89);
				legend->SetBorderSize(0);
				lambda_parameters[trig][frame][phase]->GetYaxis()->SetRangeUser(-2,2);
				if(trig==0){
					lambda_parameters[trig][frame][phase]->SetMarkerColor(kRed);
					lambda_parameters[trig][frame][phase]->SetLineColor(kRed);
					lambda_parameters_sys[trig][frame][phase]->SetMarkerColor(kRed);
					lambda_parameters_sys[trig][frame][phase]->SetLineColor(kRed);
					lambda_parameters[trig][frame][phase]->SetMarkerStyle(20);
				}
				else if(trig==1){
					lambda_parameters[trig][frame][phase]->SetMarkerColor(kBlue);
					lambda_parameters[trig][frame][phase]->SetLineColor(kBlue);
					lambda_parameters_sys[trig][frame][phase]->SetMarkerColor(kBlue);
					lambda_parameters_sys[trig][frame][phase]->SetLineColor(kBlue);
					lambda_parameters[trig][frame][phase]->SetMarkerStyle(22);
				}
				else lambda_parameters[trig][frame][phase]->SetMarkerStyle(23);
				lambda_parameters[trig][frame][phase]->SetTitle(Form("%s in %s frame; J/#psi p_{T};%s",phaseName[phase].Data(),frameName[frame].Data(),phaseName[phase].Data()));
				lambda_parameters[trig][frame][phase]->Draw("ap");
				lambda_parameters[trig][frame][phase]->Write();
				compare_lambda(trig,frame,phase);
				legend->AddEntry(lambda_parameters[trig][0][0],Form("%s",trigName[trig].Data()),"p");	
				legend->AddEntry(comparisonlambda[trig][phase+3*frame],Form("old %s",trigName[trig].Data()),"p");
				legend->Draw("same");
				lambda_parameters[trig][frame][phase]->Draw("psame");
			}
		}
		lambdacanvas->SaveAs("~/polresultspdsf/20160707/figures/lambdas_comparison_all.pdf");
		return;
	}

	lambdacanvas = new TCanvas("lambdacanvas","lambdacanvas",1200,800);
	lambdacanvas->Divide(3,2);
	for(int phase=0;phase<NPHASE+1;phase++){
		for(int frame=0;frame<NFRAME;frame++){ 
			lambdacanvas->cd(frame*3+phase+1);
			legend = new TLegend(0.7,0.7,0.89,0.89);
			legend->SetBorderSize(0);
			lambda_parameters[trig][frame][phase]->GetYaxis()->SetRangeUser(-2,2);
			if(trig==0){
				lambda_parameters[trig][frame][phase]->SetMarkerColor(kRed);
				lambda_parameters[trig][frame][phase]->SetLineColor(kRed);
				lambda_parameters_sys[trig][frame][phase]->SetMarkerColor(kRed);
				lambda_parameters_sys[trig][frame][phase]->SetLineColor(kRed);
				lambda_parameters[trig][frame][phase]->SetMarkerStyle(20);
			}
			else if(trig==1){
				lambda_parameters[trig][frame][phase]->SetMarkerColor(kBlue);
				lambda_parameters[trig][frame][phase]->SetLineColor(kBlue);
				lambda_parameters_sys[trig][frame][phase]->SetMarkerColor(kBlue);
				lambda_parameters_sys[trig][frame][phase]->SetLineColor(kBlue);
				lambda_parameters[trig][frame][phase]->SetMarkerStyle(22);
			}
			else lambda_parameters[trig][frame][phase]->SetMarkerStyle(23);
			lambda_parameters[trig][frame][phase]->SetTitle(Form("%s in %s frame; J/#psi p_{T};%s",phaseName[phase].Data(),frameName[frame].Data(),phaseName[phase].Data()));
			lambda_parameters[trig][frame][phase]->Draw("ap");
			lambda_parameters[trig][frame][phase]->Write();
			//				if(phase!=2)
			compare_lambda(trig,frame,phase);
			legend->AddEntry(lambda_parameters[trig][0][0],Form("%s",trigName[trig].Data()),"p");	
			legend->AddEntry(comparisonlambda[trig][phase+3*frame],Form("old %s",trigName[trig].Data()),"p");
			legend->Draw("same");
			lambda_parameters[trig][frame][phase]->Draw("psame");
			//				lambda_parameters_sys[trig][frame][phase]->Draw("samep[]");
		}
	}
	lambdacanvas->SaveAs(Form("~/polresultspdsf/20160707/figures/lambdas_comparison_%d.pdf",trig));
}

void plotaxissetting(int trig,int frame,int phase){
	lambda_parameters[trig][frame][phase]->SetTitle(Form("%s in %s; J/#psi p_{T};#lambda_{#theta}",trigName[trig].Data(),frameName[frame].Data()));
}

double getmaximum(std::vector<double> inputvector){
	double max;
	max = *std::max_element(&inputvector[0],&inputvector[0]+(inputvector.size()));// marker
	return max;
}

double getminimum(std::vector<double> inputvector){
	double min;
	min = *std::min_element(&inputvector[0],&inputvector[0]+(inputvector.size()));
	return min;
}

void compare_lambda(int trig,int frame, int phase){
	if(phase==0) comparisonlambda[trig][phase+3*frame] = (TGraphErrors*)outputfile->Get(Form("plot_%s_theta_%d_1eID",trigName[trig].Data(),frame));
	if(phase==1) comparisonlambda[trig][phase+3*frame] = (TGraphErrors*)outputfile->Get(Form("plot_%s_phi_%d_1eID",trigName[trig].Data(),frame));
	if(phase==2) comparisonlambda[trig][phase+3*frame] = (TGraphErrors*)outputfile->Get(Form("plot_%s_inv_%d_1eID",trigName[trig].Data(),frame));

	if(trig==0){
		comparisonlambda[trig][phase+3*frame]->SetMarkerColor(kRed);
		comparisonlambda[trig][phase+3*frame]->SetLineColor(kRed);
		comparisonlambda[trig][phase+3*frame]->SetMarkerStyle(24);
		comparisonlambda[trig][phase+3*frame]->SetMarkerSize(1.2);
	}
	if(trig==1){
		comparisonlambda[trig][phase+3*frame]->SetMarkerColor(kBlue);
		comparisonlambda[trig][phase+3*frame]->SetLineColor(kBlue);
		comparisonlambda[trig][phase+3*frame]->SetMarkerStyle(26);
		comparisonlambda[trig][phase+3*frame]->SetMarkerSize(1.2);
	}
	if(trig==2){
		comparisonlambda[trig][phase+3*frame]->SetMarkerColor(kBlack);
		comparisonlambda[trig][phase+3*frame]->SetLineColor(kBlack);
		comparisonlambda[trig][phase+3*frame]->SetMarkerStyle(27);
		comparisonlambda[trig][phase+3*frame]->SetMarkerSize(1.2);
	}
	//	if(phase!=2) 
	comparisonlambda[trig][phase+3*frame]->Draw("samep");
}


