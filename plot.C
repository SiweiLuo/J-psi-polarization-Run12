#include "TFile.h"
#include "TCanvas.h"
#include "TPDF.h"

void plot(){

	TPDF *pdf = new TPDF("/Users/siwei/polresults/20160707/pdf/parameters.pdf",111);
	//	TPDF *pdf = new TPDF("parameters.pdf",111);
	TFile* rootfile;
	TFile* combinedfile;
	TH2F* chi2;
	for(int trig=0;trig<3;trig++){
		for(int frame=0;frame<2;frame++){
			TCanvas* c1 = new TCanvas(Form("%d_%d",trig,frame),Form("%d_%d",trig,frame),1200,1200);
			c1->Divide(3,2);
			for(int pt=0;pt<6;pt++){
				combinedfile = new TFile("~/polresults/20160707/rootcombined/lambda_file0_trg0_pt0_frame0.root","read");		
				chi2 = (TH2F*)combinedfile->Get("chi2_0_0_0_0");
				c1->cd(pt+1);
				chi2->Draw("colz");		
			}
//			c1->Update();
		}
	}	
	pdf->Close();
}
