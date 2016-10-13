#include "TFile.h"
#include "TCanvas.h"
#include "TPDF.h"

void plot(){
	TFile* rootfile;
	TFile* combinedfile;
	TH2F* rawdata;
	TH2F* eff;
	TH2F* chi2;
	TH2F* truth;
	TH2F* bestfit;
	int file=0;
	for(int trig=0;trig<3;trig++){
		for(int frame=0;frame<2;frame++){
			for(int pt=0;pt<6;pt++){
				TCanvas* c1 = new TCanvas(Form("trig%d_frame%d_pt%d",trig,frame,pt),Form("trig%d_frame%d_pt%d",trig,frame,pt),2000,1200);
				c1->Divide(3,2);
				combinedfile = new TFile(Form("~/polresults/20160707/rootcombined/lambda_file0_trg%d_pt%d_frame%d.root",trig,pt,frame),"read");		
				rootfile = new TFile(Form("~/polresults/20160707/rootfiles/chi2histograms_%d/chi2histograms_0_%d_%d_%d_1_1.root",file,trig,pt,frame),"read");
				chi2 = (TH2F*)combinedfile->Get(Form("chi2_%d_%d_%d_%d",file,trig,pt,frame));
				int xx,yy,zz;
				chi2->GetMinimumBin(xx,yy,zz);
				bestfit = (TH2F*)combinedfile->Get(Form("template_file%d_trg%d_pt%d_frame%d_x%d_y%d",file,trig,pt,frame,xx,yy));
				truth = (TH2F*)combinedfile->Get(Form("theta_%d_phi_%d",xx-1,yy-1));
				rawdata = (TH2F*)rootfile->Get(Form("rawdata_%d_HT%d_pT%d_frame%d",file,trig,pt,frame));//marker add file
				eff = (TH2F*)rootfile->Get(Form("eff_2dhist_%d_%d_%d_%d",file,trig,pt,frame));	
				
				if(chi2==0x0 || bestfit==0x0 || truth==0x0 || rawdata==0x0 || eff==0x0) continue;
					
				c1->cd(1);
				truth->Draw("colz");
				c1->cd(2);
				eff->Draw("colz");
				c1->cd(3);
				chi2->Draw("colz");		
				c1->cd(4);
				bestfit->Draw("colz");
				c1->cd(5);
				rawdata->Draw("colz");	

				c1->SaveAs(Form("~/polresults/20160707/pdf/trig%d_frame%d_pt%d.pdf",trig,frame,pt));
			}
		}
	}	
}
