#include <unistd.h>
#include <stdlib.h>
#include <fstream>
#include <stddef.h>
#include <iostream>
#include <signal.h>
#include <string.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TMath.h>

using namespace std;


enum EdistribBeamType_t {mono, plain, histo, kaon, neutron, photon};
enum Jexcept { XnumberArgument, XmaxLowerMin, XfileNotOpen, XhistNotFound };
extern const double kaonmass, neutronmass, photonmass;
extern const string KaonFile, NeutronFile, PhotonFile;
extern const string histname;


class JGenBeamEnergy {
    EdistribBeamType_t typeE;
    EdistribBeamType_t typeB;
    
    
    Double_t Evalue;
    Double_t Emin;
    Double_t Emax;
    TLorentzVector P4beam;
    TRandom3 *r;
    TDirectory* genDir;
    TH1F* hgen;

    
    
    void constructor (EdistribBeamType_t typebeam_, EdistribBeamType_t typedis_, Double_t Emin_, Double_t Emax_);
    void constructor (EdistribBeamType_t typebeam_, string filename, string hstname, Double_t Emin_, Double_t Emax_);
    
    public:
    JGenBeamEnergy (char * copt);
    JGenBeamEnergy (Double_t Efixed) {
        constructor (kaon, mono, Efixed, Efixed);}
    JGenBeamEnergy (Double_t Emin_, Double_t Emax_) {
        constructor (kaon, plain, Emin_, Emax_); }
    JGenBeamEnergy (EdistribBeamType_t typeB_, EdistribBeamType_t typeD_, Double_t Emin_, Double_t Emax_);
    
    Double_t Generate ();
    TLorentzVector GetP4 () { return P4beam; }
};
