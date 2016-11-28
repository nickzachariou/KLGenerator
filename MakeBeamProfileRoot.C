#include <TGraph.h>
#include <TH1D.h>
#include "TFile.h"
#include "TF1.h"

using namespace std;

void MakeBeamProfileRoot(){
    
    int whichroot=1;// 0 to read in KL beam profile data file and create a hist; 1 to create the photon beam profile from function
    
    string loc_NewFileName="BeamProfile.root";
    TFile *locFile = new TFile(loc_NewFileName.c_str(), "RECREATE");

    TF1 *fox3 = new TF1("fox3","1.0/x*(18.3484*(1+(1-x)*(1-x))-17.6817*(1-x)*(2.0/3.0))", 0, 12);
    
    int numbbins=1;
    if (whichroot==0){
        numbbins=6000;
    TH1D *BeamProfile=new TH1D("BeamProfile","K_{L} beam profile; E_{K_{L}} (GeV); counts", numbbins,0, 6);
    }
    else if (whichroot==1) {
        numbbins=12000;
    TH1D *BeamProfile=new TH1D("BeamProfile","#gamma beam profile; E_{#gamma} (GeV); counts", numbbins,0, 12);
    }
    Double_t xval;
    
    if (whichroot==0){
        
        const char nick[]="./kl_mom_prof.dat";
        TGraph *KLbeam=new TGraph(nick,"%lg %lg"," \t");
        
        for (int i=0; i<numbbins; i++){
            xval=6.0/numbbins*i+6.0/(2*numbbins);
            BeamProfile->Fill(xval, KLbeam->Eval(xval));
            cout<<xval<<" "<<KLbeam->Eval(xval)<<endl;
        }
    }
    else if (whichroot==1){
        
        for (int i=0; i<numbbins; i++){
            xval=BeamProfile->GetBinCenter(i);
            BeamProfile->SetBinContent(i, fox3->Eval(xval/12.6));
        }
    }
    BeamProfile->Draw();
        //fox3->Draw("same");
    locFile->Write();

}
