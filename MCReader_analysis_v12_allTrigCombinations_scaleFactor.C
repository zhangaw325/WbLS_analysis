#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

using namespace std;

#include <RAT/DS/Run.hh>
#include <RAT/DS/PMTInfo.hh>
#include <RAT/DS/Root.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/MCParticle.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TGraph2D.h>
#include <TGraph.h>
#include <TMath.h>
#include <TApplication.h>
#include <TF1.h>
#include <TLegend.h>
#include <TStyle.h>

//this version is used to test my reconstruction idea with single (known) event.

int main(int argc, char **argv) {

    const int numberofinputfiles = atoi(argv[2]);//150; //148 number of input files

    double npescalefactor[8]; // scale the nep, just like there is an overall efficiency on photon detection 
                              // every PMT may have different scale factor
    for(int i=0;i<8;i++){
        npescalefactor[i] = atof(argv[3+i]);
    }
    int nbadruns = 0;
    //this file is without time information
    string headdir = "./";
    string inputrootfile[numberofinputfiles];

    fstream finnames(argv[1],ios::in);
	
    string onename;
    int cntname = 0;
    while(finnames>>onename){
        inputrootfile[cntname] = onename;
        cntname++;
	if(cntname>=numberofinputfiles) break;
    }
    finnames.close();

    for(int i=0;i<numberofinputfiles;i++) inputrootfile[i] = headdir + inputrootfile[i];
	
    std::string sDigiChMap[8]={"Digi1_CH0_S0","Digi1_CH1_S1","Digi1_CH2_S2","Digi1_CH3_S3",
                               "Digi2_CH0_S4","Digi2_CH1_S5","Digi2_CH2_S6","Digi2_CH3_S7"};
    std::string sDigiCh[8] = {"S0","S1","S2","S3","S4","S5","S6","S7"};
    std::string sTDC_ch_Map[16]={"H0","H1","H2","H3","H4","H5",
                                 "HodoTrig","MultTrig","LedTrig","AllTrig",
                                 "S0","S1","S2","S3","S4","S5"}; // here the "AllTrig" means multi_or_hodo_or_led_

    const int nbTrigs = 27;
    // trigger types:
    // I have 27 types of combinations based on the normal hodotrig defined by (H2 or H0) and (H3 or H1)
    // Using the order H2, H0, H3, H1, H4, H5
    // and 0, 1 for hit or no hit, I list the combinations as following:
    string strig[nbTrigs]={"101000", "101010", "101001",
                           "100100", "100110", "100101",
                           "011000", "011010", "011001",
                           "010100", "010110", "010101",
                           "111000", "111010", "111001",
                           "110100", "110110", "110101",
                           "011100", "011110", "011101",
                           "101100", "101110", "101101",
                           "111100", "111110", "111101"
                          };
     int indexMap_forEvtStartTime[nbTrigs]={2,2,2,2,2,2,  0,0,0,0,0,0,   2,2,2,2,2,2, 0,0,0,  2,2,2,2,2,2};
    // the above scheme, will select exclusive triggers. 
    // I ignore the case where both H4 and H5 are hit (because they are mostly caused by one muon and one scattered electron or gamma)
    // For the case where H4/H5 information should not be considered, the union of H4/H5 "01", "10", "00" will make it ("11" is ignored).

    char temprtfilechar[100]; sprintf(temprtfilechar,"Results_v12_allTrigCombinations_%s_scaleFactor_%s.root", argv[1],argv[3]);
    TFile* foutroot = new TFile(temprtfilechar,"recreate");

    // cut on plastic scintillators, 0.5 MeV
    double edepcut = 0.5;
    // PMT spe calibration
    double converterfactor[8]={146.278, 169.495, 180.849, 173.302, 165.905, 165.331, 140.525, 138.591};
    // the mean values of the pulse start time for muons, I will update the values using the arrival time of first photons in PMTs.    
    float meanPulseStartBin[8]={7.198, 7.098, 7.396, 7.904, 5.509, 4.151, 7.179, 7.026};
    double coincidenceTimeMeanMatrix[8][8]={
              {-100,       -100,       -100,      -100,      -100,     -100,      -100,      -100},
              {-0.194713,   -100,       -100,      -100,      -100,     -100,      -100,      -100},
              {-0.758341,  -0.547047,   -100,      -100,      -100,     -100,      -100,      -100},
              {-1.36467,   -1.17166,  -0.647578,  -100,      -100,     -100,      -100,      -100},
              {0.804999,   1.01294,   1.56821 ,  2.18819,     -100,     -100,      -100,      -100},
              {1.71103,   2.02743,    2.57424, 3.23104,   1.01547, -100,      -100,      -100},
              {3.09198,   -0.981527,   -0.581544, 0.166766,   -2.12371,  -3.31489, -100,      -100},
              {1.36372	, 0.18576,   -0.961526,  -0.460571,   -2.55673, -3.96895, -0.205749, -100}    
                                       };
    double coincidenceTimeWidthMatrix[8][8]={
              {-100,      -100,       -100,      -100,      -100,     -100,      -100,      -100},
              {1.14204,   -100,       -100,      -100,      -100,     -100,      -100,      -100},
              {1.14403,  0.918625,     -100,      -100,      -100,     -100,      -100,      -100},
              {1.14654,   0.882617,    0.962183,   -100,      -100,     -100,      -100,      -100},
              {1.13599,  0.906206,     1.10838,   1.14001,   -100,     -100,      -100,      -100},
              {1.37561,  1.36813,     1.49745,    1.6342,  1.78241,  -100,      -100,      -100},
              {1.02948,  0.985831,      1.00617,   0.900906,   0.923113,  0.902838,   -100,      -100},
              {2.94202,  3.40961,      1.38917,   1.47374,   1.51825,   1.60992,   0.843143,   -100}    
                                       };
    // directories - by trigger type
    foutroot->cd();
    TDirectory* dir[nbTrigs];
    char dirnames[100];
    for(int i=0; i<nbTrigs; i++){
        sprintf(dirnames,"dir_%s",strig[i].c_str());
        dir[i] = foutroot->mkdir(dirnames);
        dir[i]->cd();
    }
    TDirectory* dir_mcwaveforms = foutroot->mkdir("MCWaveforms");
    const int pmtnb = 8;
    char hName[200];
    char tempname[150];
    //input muons's momentum and angles
    TH1F* hMuMomentum[nbTrigs];
    TH1F* hMuPAngle[nbTrigs];
    TH1F* hMuAAngle[nbTrigs];
    TH1F* hDecayZ[nbTrigs];
    TH2F* hEndXY[nbTrigs]; // end x and y coordinates at Z=H4/H5 
    TH2F* hEndXY0[nbTrigs];// end x and y coordinates when entering liquid (under the top acrylic lid)
    TGraph* gNumberOfEvents = new TGraph();
    //TH2F* hEndXY1[nbTrigs];// end x and y coordinates when leaving liquid (under the bottom acrylic lid)
    for(int j=0;j<nbTrigs;j++){
        sprintf(hName,"hMuMomentum_%sTrig",strig[j].c_str());
        hMuMomentum[j] = new TH1F(hName,"",5000,0,50000);
        hMuMomentum[j]->SetXTitle("Momentum (MeV/c)");
        hMuMomentum[j]->SetYTitle("Counts");
        sprintf(hName,"hMuPolarAngle_%sTrig",strig[j].c_str());
        hMuPAngle[j] = new TH1F(hName,"",170,0.5*TMath::Pi(),1.2*TMath::Pi());
        hMuPAngle[j]->SetXTitle("Polar angle (rad.)");
        hMuPAngle[j]->SetYTitle("Counts");
        sprintf(hName,"hMuAzimuthalAngle_%sTrig",strig[j].c_str());
        hMuAAngle[j] = new TH1F(hName,"",100, -1.0*TMath::Pi(), 1.0*TMath::Pi());
        hMuAAngle[j]->SetXTitle("Azimuthal angle (rad.)");
        hMuAAngle[j]->SetYTitle("Counts");

        sprintf(hName,"hDecayZ_%s",strig[j].c_str());
        hDecayZ[j] = new TH1F(hName,"",240, -1200, 1200);
        hDecayZ[j]->SetXTitle("Z coordinate (mm)");
        hDecayZ[j]->SetYTitle("Counts");
        sprintf(hName,"hEndXY_%s_zAtH4H5",strig[j].c_str()); // at z corrdinate where H4 and H5 sit
        hEndXY[j] = new TH2F(hName,"",100, -1000, 1000, 100, -1000, 1000);
        hEndXY[j]->SetXTitle("X (mm)");
        hEndXY[j]->SetYTitle("Y (mm)");
        sprintf(hName,"hEndXY_%s_z-500",strig[j].c_str()); // at z corrdinate where muon enters the liquid (z=800.4)
        hEndXY0[j] = new TH2F(hName,"",100, -1000, 1000, 100, -1000, 1000);
        hEndXY0[j]->SetXTitle("X (mm)");
        hEndXY0[j]->SetYTitle("Y (mm)");
        //sprintf(hName,"hEndXY_%s_z-500",strig[j].c_str()); // at z corrdinate where muon enters the liquid (z=800.4)
        //hEndXY1[j] = new TH2F(hName,"",100, -1000, 1000, 100, -1000, 1000);
        //hEndXY1[j]->SetXTitle("X (mm)");
        //hEndXY1[j]->SetYTitle("Y (mm)");
    }

    // edep and hit time in the 6 plastic scintillators
    TH1F* hEdepInPS[6][nbTrigs];
    TH1F* hHitTimeInPS[6][nbTrigs];
    for(int i=0; i<6; i++){
        for(int j=0;j<nbTrigs;j++){
            sprintf(hName,"hEdepInPs_%d_%sTrig",i,strig[j].c_str());
            hEdepInPS[i][j]=new TH1F(hName,"",500,0,10);
            hEdepInPS[i][j]->SetXTitle("Energy (MeV)");
            hEdepInPS[i][j]->SetYTitle("Counts");
            sprintf(hName,"hHitTimeInPs_%d_%sTrig",i,strig[j].c_str());
            hHitTimeInPS[i][j]=new TH1F(hName,"",500,0,50);
            hHitTimeInPS[i][j]->SetXTitle("Time (ns)");
            hHitTimeInPS[i][j]->SetYTitle("Counts");
        }
    }
    // the event time using the hodoscope detector's time.
    TH1F* hEvtTime[nbTrigs]; 
    TH1F* hTrigTime[nbTrigs];
    for(int j=0;j<nbTrigs;j++){
        sprintf(tempname,"hEvtTime_PS_%s",strig[j].c_str());
        hEvtTime[j] = new TH1F(tempname,"",500,0,10);
        hEvtTime[j]->SetXTitle("Event begin time (ns)");
        hEvtTime[j]->SetYTitle("Number of events");
        sprintf(tempname,"hEvtTriggerTime_%s",strig[j].c_str());
        hTrigTime[j] = new TH1F(tempname,"",500,0,10);
        hTrigTime[j]->SetXTitle("Event begin time (ns)");
        hTrigTime[j]->SetYTitle("Number of events");
    }
	
    //hit time differences in hodoscope detectors
    TH1F* hTDC_H0H2 = new TH1F("hTDC_H0_H2","TDC_H0 - TDC_H2", 1000,-10,10);
    TH1F* hTDC_H1H2 = new TH1F("hTDC_H1_H2","TDC_H1 - TDC_H2", 1000,-10,10);
    TH1F* hTDC_H3H2 = new TH1F("hTDC_H3_H2","TDC_H3 - TDC_H2", 1000,-10,10);
    TH1F* hTDC_H1H0 = new TH1F("hTDC_H1_H0","TDC_H1 - TDC_H0", 1000,-10,10);
    TH1F* hTDC_H3H0 = new TH1F("hTDC_H3_H0","TDC_H3 - TDC_H0", 1000,-10,10);
    TH1F* hTDC_H5H2 = new TH1F("hTDC_H5_H2","TDC_H5 - TDC_H2", 1000,00,20);
    TH1F* hTDC_H4H2 = new TH1F("hTDC_H4_H2","TDC_H4 - TDC_H2", 1000,00,20);
    TH1F* hTDC_H5H0 = new TH1F("hTDC_H5_H0","TDC_H5 - TDC_H0", 1000,00,20);
    TH1F* hTDC_H4H0 = new TH1F("hTDC_H4_H0","TDC_H4 - TDC_H0", 1000,00,20);
    
    TH1F* hTDC_H4H5 = new TH1F("hTDC_H4_H5","TDC_H4 - TDC_H5", 400,-20,20);
    TH2F* hEdep_vs_TDC_H1H2 = new TH2F("hEdep_vs_TDC_H1H2","",1000,-10,10,500,0,10);
    TH2F* hEdepH4_vs_EdepH5 = new TH2F("hEdepH4_vs_EdepH5","",50,0,10,50,0,10);
	
    //charge, time in PMTs one by one
    TH1F* hCharge[pmtnb][nbTrigs]; // the total charge in every pmt within a time window of 2.5us (from event start time, 0 ns, in the MC: muon starts at Z=1200 mm).
    TH1F* hCharge1[pmtnb][nbTrigs]; // the total charge in every pmt within a time window of 50ns. This might be useful
    TH1F* hCharge2[pmtnb][nbTrigs]; // the total charge in every pmt within a time window of 40ns. This might be useful
    TH1F* hTime[pmtnb][nbTrigs];  // time of arrival of the first photon at a pmt
    TH2F* hPMTNpeVsTime[pmtnb][nbTrigs]; // 2D plot of time of photon vs. charge
    TH1F* hPMTPhotonTimeSpan[pmtnb][nbTrigs]; //  time of all photons in each pmt
    for(int i=0;i<pmtnb;i++){
        for(int j=0;j<nbTrigs;j++){
            sprintf(hName,"hCharge_2500ns_S%d_%sTrig",i,strig[j].c_str());
            hCharge[i][j] = new TH1F(hName,"",50,0,50);
            hCharge[i][j]->SetXTitle("Number of detected photons");
            hCharge[i][j]->SetYTitle("Counts");
            sprintf(hName,"hCharge_50ns_S%d_%sTrig",i,strig[j].c_str());
            hCharge1[i][j] = new TH1F(hName,"",50,0,50);
            hCharge1[i][j]->SetXTitle("npe in 50 ns window after the first photon");
            hCharge1[i][j]->SetYTitle("Counts");
            sprintf(hName,"hCharge_40ns_S%d_%sTrig",i,strig[j].c_str());
            hCharge2[i][j] = new TH1F(hName,"",50,0,50);
            hCharge2[i][j]->SetXTitle("npe in 40 ns window after the first photon");
            hCharge2[i][j]->SetYTitle("Counts");
            sprintf(hName,"hTime_S%d_%sTrig",i,strig[j].c_str());
            hTime[i][j]=new TH1F(hName,"",500,0,50); // here I only lock myself with a 2560 ns window
            hTime[i][j]->SetXTitle("Time of the first photon arriving at PMT (ns)");
            hTime[i][j]->SetYTitle("Counts");
            sprintf(hName,"hPMTNpeVsFirstPhotonTime_S%d_%sTrig",i,strig[j].c_str());
            hPMTNpeVsTime[i][j]=new TH2F(hName,"",520,-2,50,2560,-0.5,2560-0.5);
            hPMTNpeVsTime[i][j]->SetXTitle("npe");
            hPMTNpeVsTime[i][j]->SetYTitle("Time of photons arriving at PMT (ns)");	
            sprintf(hName,"hPMTAllPhotonTime_S%d_%sTrig",i,strig[j].c_str());
            hPMTPhotonTimeSpan[i][j]=new TH1F(hName,"",500,0,5000);
            hPMTPhotonTimeSpan[i][j]->SetXTitle("Time difference between last and the first photons (ns)");
            hPMTPhotonTimeSpan[i][j]->SetYTitle("Counts");			
	}
    }
    // total charge from all PMTs
    TH1F* hTotalCharge[nbTrigs]; // total charge of all pmts in npe
    TH1F* hTotalChargeTop[nbTrigs]; // total charge of top pmts in npe
    TH1F* hTotalChargeBottom[nbTrigs]; // total charge of bottom pmts in npe
    TH1F* hNpeRatio[nbTrigs]; // charge ratio: total charge in top PMTs over total charge in ALL PMTs, 
	                          // weighted by number of PMTs.
    for(int j=0;j<nbTrigs;j++){
        sprintf(hName,"hTotalCharge_%sTrig",strig[j].c_str());
        hTotalCharge[j] = new TH1F(hName,"",1000,0,1000);
        hTotalCharge[j]->SetXTitle("Number of detected photons");
        hTotalCharge[j]->SetYTitle("Counts");
        sprintf(hName,"hTotalChargeTop_%sTrig",strig[j].c_str());
        hTotalChargeTop[j] = new TH1F(hName,"",1000,0,1000);
        hTotalChargeTop[j]->SetXTitle("Number of detected photons");
        hTotalChargeTop[j]->SetYTitle("Counts");
        sprintf(hName,"hTotalChargeBottom_%sTrig",strig[j].c_str());
        hTotalChargeBottom[j] = new TH1F(hName,"",1000,0,1000);
        hTotalChargeBottom[j]->SetXTitle("Number of detected photons");
        hTotalChargeBottom[j]->SetYTitle("Counts");
	sprintf(hName,"hNpeRatio_TopBotPMTs_%sTrig",strig[j].c_str());
        hNpeRatio[j] = new TH1F(hName,"",100,0,1);
        hNpeRatio[j]->SetXTitle("Charge ratio: charge in top PMTs over charge in all PMTs");
        hNpeRatio[j]->SetTitle("Weighted by number of PMTs");
    }
	
    TH1F* hPulseTimeDiffToS0[7][nbTrigs]; // the difference in pulse times between other PMTs and S0
    TH1F* hPulseTimeDiffToS1[6][nbTrigs]; // the difference in pulse times between other PMTs and S1
    TH1F* hPulseTimeDiffToS2[5][nbTrigs]; // the difference in pulse times between other PMTs and S2
    TH1F* hPulseTimeDiffToS3[4][nbTrigs]; // the difference in pulse times between other PMTs and S3
    TH1F* hPulseTimeDiffToS4[3][nbTrigs]; // the difference in pulse times between other PMTs and S4
    TH1F* hPulseTimeDiffToS5[2][nbTrigs]; // the difference in pulse times between other PMTs and S5
    TH1F* hPulseTimeDiffToS6[1][nbTrigs]; // the difference in pulse times between other PMTs and S6
    for(int i=0; i<pmtnb; i++){
        for(int j=0;j<nbTrigs;j++){
            sprintf(tempname,"hPulseTimeDiffToS0_%s_%s",sDigiChMap[i].c_str(),strig[j].c_str());
            hPulseTimeDiffToS0[i-1][j] = new TH1F(tempname,"",2560*4,-2560,2560);
            sprintf(tempname,"Pulse time differencce between %s and S0", sDigiChMap[i].c_str());
            hPulseTimeDiffToS0[i-1][j]->SetXTitle(tempname);
            hPulseTimeDiffToS0[i-1][j]->SetYTitle("Counts / 1 ns");
            if(i>1){
              sprintf(tempname,"hPulseTimeDiffToS1_%s_%s",sDigiChMap[i].c_str(),strig[j].c_str());
              hPulseTimeDiffToS1[i-2][j] = new TH1F(tempname,"",2560*4,-2560,2560);
              sprintf(tempname,"Pulse time differencce between %s and S1", sDigiChMap[i].c_str());
              hPulseTimeDiffToS1[i-2][j]->SetXTitle(tempname);
              hPulseTimeDiffToS1[i-2][j]->SetYTitle("Counts / 1 ns");
            }
            if(i>2){
              sprintf(tempname,"hPulseTimeDiffToS2_%s_%s",sDigiChMap[i].c_str(),strig[j].c_str());
              hPulseTimeDiffToS2[i-3][j] = new TH1F(tempname,"",2560*4,-2560,2560);
              sprintf(tempname,"Pulse time differencce between %s and S2", sDigiChMap[i].c_str());
              hPulseTimeDiffToS2[i-3][j]->SetXTitle(tempname);
              hPulseTimeDiffToS2[i-3][j]->SetYTitle("Counts / 1 ns");
            }
            if(i>3){
              sprintf(tempname,"hPulseTimeDiffToS3_%s_%s",sDigiChMap[i].c_str(),strig[j].c_str());
              hPulseTimeDiffToS3[i-4][j] = new TH1F(tempname,"",2560*4,-2560,2560);
              sprintf(tempname,"Pulse time differencce between %s and S3", sDigiChMap[i].c_str());
              hPulseTimeDiffToS3[i-4][j]->SetXTitle(tempname);
              hPulseTimeDiffToS3[i-4][j]->SetYTitle("Counts / 1 ns");
            }
            if(i>4){
              sprintf(tempname,"hPulseTimeDiffToS4_%s_%s",sDigiChMap[i].c_str(),strig[j].c_str());
              hPulseTimeDiffToS4[i-5][j] = new TH1F(tempname,"",2560*4,-2560,2560);
              sprintf(tempname,"Pulse time differencce between %s and S4", sDigiChMap[i].c_str());
              hPulseTimeDiffToS4[i-5][j]->SetXTitle(tempname);
              hPulseTimeDiffToS4[i-5][j]->SetYTitle("Counts / 1 ns");
            }
            if(i>5){
              sprintf(tempname,"hPulseTimeDiffToS5_%s_%s",sDigiChMap[i].c_str(),strig[j].c_str());
              hPulseTimeDiffToS5[i-6][j] = new TH1F(tempname,"",2560*4,-2560,2560);
              sprintf(tempname,"Pulse time differencce between %s and S5", sDigiChMap[i].c_str());
              hPulseTimeDiffToS5[i-6][j]->SetXTitle(tempname);
              hPulseTimeDiffToS5[i-6][j]->SetYTitle("Counts / 1 ns");
            }
            if(i>6){
              sprintf(tempname,"hPulseTimeDiffToS6_%s_%s",sDigiChMap[i].c_str(),strig[j].c_str());
              hPulseTimeDiffToS6[i-7][j] = new TH1F(tempname,"",2560*4,-2560,2560);
              sprintf(tempname,"Pulse time differencce between %s and S6", sDigiChMap[i].c_str());
              hPulseTimeDiffToS6[i-7][j]->SetXTitle(tempname);
              hPulseTimeDiffToS6[i-7][j]->SetYTitle("Counts / 1 ns");
            }
        }
    } // end preparing histograms

    //some variables
    int test_cnt1 = 0, test_cnt2 = 0;
    bool istrigger[nbTrigs];
    int  trigcnts[nbTrigs];
    int trigcnts_hodo = 0;
    int trigcnts_hodo_H5_notH4=0;
    int trigcnts_hodo_notH5_notH4=0;
    int trigcnts_hodo_notH5_H4=0;
    int trigcnts_hodo_H5_H4=0;
    int trigcnts_multi=0;
    int trigcnts_decays = 0;
    int total_number_of_events = 0;
    double eventStartTime[nbTrigs];
    for(int index=0;index<nbTrigs;index++){
      istrigger[index]=false;
      eventStartTime[index]=0;
      trigcnts[index]=0;
    }
    int testcntTtrigtype1[8]; for(int t=0;t<8;t++) testcntTtrigtype1[t]=0;

    //Getting ready to look through ROOT Tree files
    TFile* f;
    TTree* T;
    TTree* runT;
    for(int ffff=0; ffff<numberofinputfiles; ffff++){ // start loop: reading the root tree files
        //open the root file
        f = new TFile(inputrootfile[ffff].c_str(),"read");
        std::cout<<"Reading "<<inputrootfile[ffff]<<std::endl;
        T = (TTree*) f->Get("T");
        if(T==NULL){
            cout<<"bad tree "<<inputrootfile[ffff]<<endl;
            continue;
        }
        runT = (TTree*) f->Get("runT");
        RAT::DS::Run *run = new RAT::DS::Run();
        runT->SetBranchAddress("run", &run);
        if (runT->GetEntries() != 1) {
          cout << "Funny run tree, aborting" << endl;
          return -1;
        }
        runT->GetEntry(0);
        RAT::DS::PMTInfo *pmtinfo = run->GetPMTInfo();
        RAT::DS::Root *ds = new RAT::DS::Root();
        T->SetBranchAddress("ds", &ds);
        int nEvents = T->GetEntries();
        cout << "Reading in " << nEvents << " events" << endl;

        //reading in events
        for (int i = 0; i < nEvents; i++) {  
            total_number_of_events++;
            //if(i%100==0) cout<<"Event "<<i<<endl;
            T->GetEntry(i);
            if (ds->GetEVCount() != 1) {
                cout << "EV " << i << " is multi-trigger, ignoring" << endl;
                continue;
            }
        
            RAT::DS::MC *mc = ds->GetMC();
            RAT::DS::EV *ev = ds->GetEV(0);

            std::vector<int> myTriggerFlag = ev->GetTriggerFlag();
            std::vector<double> myTriggerEdep = ev->GetTriggerEdep();
            std::vector<double> myTriggerHitTime = ev->GetTriggerHitTime();

            double t_hodo_1 = 0, t_hodo_2=0;
            // determine trigger time in H0 and H2 
            if(myTriggerHitTime[0]>0 && myTriggerHitTime[2]>0){t_hodo_1 = (myTriggerHitTime[0]+myTriggerHitTime[2])/2.0; }
            else if(myTriggerHitTime[0]>0 && myTriggerHitTime[2]==0){t_hodo_1 = myTriggerHitTime[0]; }
            else if(myTriggerHitTime[0]==0 && myTriggerHitTime[2]>0){t_hodo_1 = myTriggerHitTime[2]; }
            // determine trigger time in H1 and H3
            if(myTriggerHitTime[1]>0 && myTriggerHitTime[3]>0){t_hodo_2 = (myTriggerHitTime[1]+myTriggerHitTime[3])/2.0; }
            else if(myTriggerHitTime[1]>0 && myTriggerHitTime[3]==0){t_hodo_2 = myTriggerHitTime[1]; }
            else if(myTriggerHitTime[1]==0 && myTriggerHitTime[3]>0){t_hodo_2 = myTriggerHitTime[3]; }
            // take the earlier time in t_hodo_1 and t_hodo_2 as the trigger time
            double eventTime = (t_hodo_1<t_hodo_2)?t_hodo_1:t_hodo_2;
            double delta_t_hodoDet = t_hodo_2 - t_hodo_1;
            // time of the first photon from muon, time of the first photon from decays        
            double t_top_first=0, t_bot_first=0;

            int nPrim = mc->GetMCParticleCount();
            int npe = mc->GetNumPE();
            //MCParticle
            TVector3 evpos = mc->GetMCParticle(0)->GetPosition(); // position in millimeter
            int pdgcode = mc->GetMCParticle(0)->GetPDGCode();     // particle pdg code
            string particlename = mc->GetMCParticle(0)->GetParticleName(); // particle name
            float evtime = mc->GetMCParticle(0)->GetTime(); // particle start time, it is zero for every event!!!

            float ke = mc->GetMCParticle(0)->GetKE(); // particle ke in MeV
            TVector3 momentum = mc->GetMCParticle(0)->GetMomentum(); // particle momentum in MeV/c
            TVector3 polarization = mc->GetMCParticle(0)->GetPolarization(); // particle polarization, not used in my study
            float TrackStartX=evpos.x(), TrackStartY=evpos.y(), TrackStartZ=evpos.z();
            float TrackPx = momentum.x(), TrackPy=momentum.y(), TrackPz=momentum.z();
            float MomentumSquare = TrackPx*TrackPx+TrackPy*TrackPy+TrackPz*TrackPz;
            float Mass = (MomentumSquare-ke*ke)/(2.0*ke);

            float TrackSpeed=TMath::Sqrt(1.0/(1.0+Mass*Mass/MomentumSquare)); // speed of the input particle, I know it since I know its momentum
            float X1=evpos.x(), Y1=evpos.y(), Z1=evpos.z();
            float Xprime=TrackPx/TrackPz, Yprime=TrackPy/TrackPz;
            float Xzero=X1-Xprime*Z1, Yzero=Y1-Yprime*Z1;
            float polarAngle = TMath::Pi() - TMath::ATan(TMath::Sqrt(Xprime*Xprime+Yprime*Yprime));
            float azimuthalAngle = TMath::ATan(Yprime/Xprime);
            // for track reconstruction later
            float Za0 = -500.4+25.4; //800.4-25.4; // Point A at the water surface, where track goes through A
            float Xa0 = Xzero + Xprime*Za0, Ya0 = Yzero + Yprime*Za0;            
            float Za = -1118.8075; // Point B at the surface of H4 and/or H5
            float Xa = Xzero + Xprime*Za, Ya = Yzero + Yprime*Za;
            float distance1 = TMath::Sqrt((TrackStartX-Xa)*(TrackStartX-Xa)+(TrackStartY-Ya)*(TrackStartY-Ya)+(TrackStartZ-Za)*(TrackStartZ-Za));
            float time_stage1 = distance1/(TrackSpeed*300.0); // in ns

            /*
            cout<<"Particle\t"<<particlename<<"\tmass\t"<<Mass<<endl;
            cout<<"Position (mm)\t"<<evpos.x()<<"\t"<<evpos.y()<<"\t"<<evpos.z()<<endl;
            cout<<"Ke (MeV)\t"<<ke<<endl;
            cout<<"Momentum (MeV)\t"<<momentum.x()<<"\t"<<momentum.y()<<"\t"<<momentum.z()<<endl;
            cout<<"Speed (c)\t"<<TrackSpeed<<endl;
            */ 

            //TVector3 energyCentroid = mc->GetMCSummary()->GetEnergyCentroid(); 
            /**
            * Centroid of energy loss.
            *
            * This is the average position of all steps in this event, weighted by
            * the energy lost in that step. Optical photons, and the rest mass of
            * particles when a track terminates are not included.
            */
            //TVector3 energyRMS = mc->GetMCSummary()->GetEnergyRMS();
            //float totalScintEdep = mc->GetMCSummary()->GetTotalScintEdep(); // in MeV
			
            //MCTrack
            int mctrackcnt = mc->GetMCTrackCount();
            double decayElectronEnergy=Mass;
            int mydecayflag = 0;
            double stopx=-1000, stopy=-1000, stopz=-10000, gTime=-1;
	    //cout<<"mctrackcnt="<<mctrackcnt<<endl;
            if(mctrackcnt!=0){ // this means a decay happend

                for(int imctrackcnt=0;imctrackcnt<mctrackcnt;imctrackcnt++){
                    int mctrackstepcnt = mc->GetMCTrack(imctrackcnt)->GetMCTrackStepCount();
					//fstream fout("track_output_test.txt",ios::out | ios::app);
                    //fout<<imctrackcnt<<"\t"<<mc->GetMCTrack(imctrackcnt)->GetID()<<"\t"<<mc->GetMCTrack(imctrackcnt)->GetParentID()<<"\t"<<mctrackstepcnt<<endl;
                    //fout<< "---->Track " << imctrackcnt << " id " << mc->GetMCTrack(imctrackcnt)->GetID() 
                    //    << " pName " << mc->GetMCTrack(imctrackcnt)->GetParticleName() << " parent " << mc->GetMCTrack(imctrackcnt)->GetParentID()<< " stepcnt "<< mctrackstepcnt <<endl;
                    //fout<<"------->stepcnt "<<mctrackstepcnt<<endl;
                    decayElectronEnergy -= mc->GetMCTrack(imctrackcnt)->GetMCTrackStep(0)->GetKE();
                    if(imctrackcnt==0 && mydecayflag==0){
                        stopx = mc->GetMCTrack(0)->GetMCTrackStep(0)->GetEndpoint().x();
                        stopy = mc->GetMCTrack(0)->GetMCTrackStep(0)->GetEndpoint().y();
                        stopz = mc->GetMCTrack(0)->GetMCTrackStep(0)->GetEndpoint().z();
                        gTime = mc->GetMCTrack(0)->GetMCTrackStep(0)->GetGlobalTime();
                        mydecayflag=1;
                    }

                    //for(int jmcstepcnt=0; jmcstepcnt<mctrackstepcnt; jmcstepcnt++){
                    //    TVector3 stepEndPoint = mc->GetMCTrack(imctrackcnt)->GetMCTrackStep(jmcstepcnt)->GetEndpoint();
                    //    fout<< "----------> steps " << jmcstepcnt << " ke "<< mc->GetMCTrack(imctrackcnt)->GetMCTrackStep(jmcstepcnt)->GetKE()
                    //        << " end_points (" <<stepEndPoint[0]<<", "<<stepEndPoint[1]<<", "<<stepEndPoint[2]<<")" << " time " << mc->GetMCTrack(imctrackcnt)->GetMCTrackStep(jmcstepcnt)->GetGlobalTime()<<" process "<< mc->GetMCTrack(imctrackcnt)->GetMCTrackStep(jmcstepcnt)->GetProcess() << " volume " << mc->GetMCTrack(imctrackcnt)->GetMCTrackStep(jmcstepcnt)->GetVolume() <<endl;
					//}
                }
            } //

            // define trigger types
            int trigflag = 0;
            int trigflag_sub = 0;
            int nHits = ev->GetPMTCount();
            //if(nHits>6){ // multiplicity triggers
                //  trigflag = 1; istrigger[trigflag]=true; trigcnts[trigflag]++;
                //hEvtTime_multiTrig->Fill(eventTime);
            //}
            // through going triggers

            if( myTriggerEdep[2]==0 && myTriggerEdep[0]>0 && myTriggerEdep[3]>0 && myTriggerEdep[1]>0 ){
                trigcnts_hodo++;
                if(myTriggerEdep[4]==0 && myTriggerEdep[5]>0.0) trigcnts_hodo_H5_notH4++;
                if(myTriggerEdep[4]>0. && myTriggerEdep[5]==0.) trigcnts_hodo_notH5_H4++;
                if(myTriggerEdep[4]==0. && myTriggerEdep[5]==0.) trigcnts_hodo_notH5_notH4++;
                if(myTriggerEdep[4]>0. && myTriggerEdep[5]>0.){
                    trigcnts_hodo_H5_H4++;
                    if(mydecayflag==1) trigcnts_decays++;
                    hEdepH4_vs_EdepH5->Fill(myTriggerEdep[4],myTriggerEdep[5]);
                    //fstream testout("testout.txt",ios::out|ios::app);
                    //testout<<ffff<<"\t"<<i<<"\t"<<pdgcode<<"\t"<<X1<<"\t"<<Y1<<"\t"<<TrackPx<<"\t"<<TrackPy<<"\t"<<TrackPz<<endl;
                    //testout<<"\t\t Edep in H4, H5: "<<myTriggerEdep[4]<<"\t"<<myTriggerEdep[5]<<endl;
                    //testout.close();
                    hTDC_H4H5->Fill(myTriggerHitTime[4]-myTriggerHitTime[5]);
                }
            }

            //now make my trigger conbinations based on the normal hodo trigger
            int PsTrigFlag[6]={0,0,0,0,0,0};
            for(int pscnt=0; pscnt<6; pscnt++){
              if(myTriggerEdep[pscnt]>0){
                PsTrigFlag[pscnt]=1;
              }
              else PsTrigFlag[pscnt]=0;
            }

            int thistriggerflag =  PsTrigFlag[2]*TMath::Power(2,5)
                                  +PsTrigFlag[0]*TMath::Power(2,4)
                                  +PsTrigFlag[3]*TMath::Power(2,3)
                                  +PsTrigFlag[1]*TMath::Power(2,2)
                                  +PsTrigFlag[4]*TMath::Power(2,1)
                                  +PsTrigFlag[5]*TMath::Power(2,0);
            //cout<<ffff<<"\t"<<i<<"\t"<<PsTrigFlag[2]<<PsTrigFlag[0]<<PsTrigFlag[3]<<PsTrigFlag[1]<<PsTrigFlag[4]<<PsTrigFlag[5]<<endl;
            for(int trigindex=0; trigindex<nbTrigs; trigindex++){
                int testtrignum = std::strtol(strig[trigindex].c_str(), 0, 2);//astring.Atoi();
                //cout<<strig[trigindex]<<"\t"<<testtrignum<<"\t"<<thistriggerflag<<endl;
                if(thistriggerflag == testtrignum){
                    hEndXY[trigindex]->Fill(Xa, Ya); hEndXY0[trigindex]->Fill(Xa0, Ya0);
                    istrigger[trigindex]=true;
                    trigcnts[trigindex]++;
                    eventStartTime[trigindex] = myTriggerHitTime[indexMap_forEvtStartTime[trigindex]];
                    hTrigTime[trigindex]->Fill(eventStartTime[trigindex]);
                    break;
                }
            }

	    bool testTrig = false;			
	    //define triggers 
            if( (myTriggerEdep[2]>0 || myTriggerEdep[0]>0) && (myTriggerEdep[3]>0 || myTriggerEdep[1]>0)  // the normal hodoTrig
              ){
            //select in-time events only, time coincidence is made by constraining time difference between hit times in plastic scintillators
                int index1 = -1, index2=-1;
                if(myTriggerEdep[0]>0 && myTriggerEdep[2]<=0) index1 = 0;
                else index1=2;
                if(myTriggerEdep[3]>0 && myTriggerEdep[1]<=0) index2 = 3;
                else index2=1;
                //if(  myTriggerHitTime[index2] - myTriggerHitTime[index1] > 0 
                //  && myTriggerHitTime[index2] - myTriggerHitTime[index1] < 10
                //  //&& myTriggerHitTime[5]-myTriggerHitTime[index1] > 0
                //  )
                testTrig = true;
                				
                trigflag = 0;
                if(testTrig == true){
                  //istrigger[trigflag]=true; trigcnts[trigflag]++;
                  hEvtTime[trigflag]->Fill(eventTime);
                }

		if( 1 ) {
		    // get time difference between plastic scintillators
   			if(myTriggerEdep[2]>0 && myTriggerEdep[0]>0){
	    			if(myTriggerHitTime[2]<=0 || myTriggerHitTime[0]<=0) cout<<"Interesting mismatch edep>0 and hitTime=0"<<endl;
		    		hTDC_H0H2->Fill(myTriggerHitTime[0] - myTriggerHitTime[2]);
    			}
	    		if(myTriggerEdep[2]>0 && myTriggerEdep[3]>0){
		    		if(myTriggerHitTime[2]<=0 || myTriggerHitTime[3]<=0) cout<<"Interesting mismatch edep>0 and hitTime=0"<<endl;
			    	hTDC_H3H2->Fill(myTriggerHitTime[3] - myTriggerHitTime[2]);
       			}
		    	if(myTriggerEdep[2]>0 && myTriggerEdep[1]>0){
			    	if(myTriggerHitTime[2]<=0 || myTriggerHitTime[1]<=0) cout<<"Interesting mismatch edep>0 and hitTime=0"<<endl;
				    hTDC_H1H2->Fill(myTriggerHitTime[1] - myTriggerHitTime[2]);
					hEdep_vs_TDC_H1H2->Fill(myTriggerHitTime[1] - myTriggerHitTime[2], myTriggerEdep[2]);
    			}
	    		if(myTriggerEdep[2]>0 && myTriggerEdep[4]>0){
		    		if(myTriggerHitTime[2]<=0 || myTriggerHitTime[4]<=0) cout<<"Interesting mismatch edep>0 and hitTime=0"<<endl;
			    	hTDC_H4H2->Fill(myTriggerHitTime[4] - myTriggerHitTime[2]);
    			}
	    		if(myTriggerEdep[2]>0 && myTriggerEdep[5]>0){
		    		if(myTriggerHitTime[2]<=0 || myTriggerHitTime[5]<=0) cout<<"Interesting mismatch edep>0 and hitTime=0"<<endl;
			    	hTDC_H5H2->Fill(myTriggerHitTime[5] - myTriggerHitTime[2]);
			    }
    			if(myTriggerEdep[0]>0 && myTriggerEdep[3]>0){
	    			if(myTriggerHitTime[0]<=0 || myTriggerHitTime[3]<=0) cout<<"Interesting mismatch edep>0 and hitTime=0"<<endl;
		    		hTDC_H3H0->Fill(myTriggerHitTime[3] - myTriggerHitTime[0]);
			    }
    			if(myTriggerEdep[0]>0 && myTriggerEdep[1]>0){
	    			if(myTriggerHitTime[0]<=0 || myTriggerHitTime[1]<=0) cout<<"Interesting mismatch edep>0 and hitTime=0"<<endl;
		    		hTDC_H1H0->Fill(myTriggerHitTime[1] - myTriggerHitTime[0]);
			    }
    			if(myTriggerEdep[0]>0 && myTriggerEdep[4]>0){
	    			if(myTriggerHitTime[0]<=0 || myTriggerHitTime[4]<=0) cout<<"Interesting mismatch edep>0 and hitTime=0"<<endl;
		    		hTDC_H4H0->Fill(myTriggerHitTime[5] - myTriggerHitTime[0]);
    			}
	    		if(myTriggerEdep[0]>0 && myTriggerEdep[5]>0){
		    		if(myTriggerHitTime[0]<=0 || myTriggerHitTime[5]<=0) cout<<"Interesting mismatch edep>0 and hitTime=0"<<endl;
			    	hTDC_H5H0->Fill(myTriggerHitTime[5] - myTriggerHitTime[0]);
			    }
    		}				
            }
	    testTrig=false;
			
            // fill histograms of (1) edep and hittime in plastic scintillators for the different trigger types
            //                    (2) number of hitting PMTs, for the different trigger types
            //                    (3) the initial muons momentum and angles.
            for(int index=0; index<nbTrigs; index++){
                if(istrigger[index]==true){
                    hMuMomentum[index]->Fill(TMath::Sqrt(MomentumSquare));
                    hMuPAngle[index]->Fill(polarAngle);
                    hMuAAngle[index]->Fill(azimuthalAngle);
                    for(int psid=0; psid<6; psid++){
                        hEdepInPS[psid][index]->Fill(myTriggerEdep[psid]);
                        hHitTimeInPS[psid][index]->Fill(myTriggerHitTime[psid]);
                    }
                }
            }

            //Dealing with photons in PMTs. Getting event info from MCphoton
            double totQ1 = 0; // total charge of all pmts in npe, cut photon arrival time to <2500 ns
            double totQ_top1 = 0, totQ_bot1 = 0;
            // use vectors to hold pulses (charge and time) in the PMTs:
	    	// pmtnb=8 columns for the 8 PMTs, each contains 0 photons when initialized
            std::vector< std::vector<double> > pulseTimeInPMTs1; // pulse time
            pulseTimeInPMTs1.resize(pmtnb,std::vector<double>(0,0)); 
            std::vector< std::vector<double> > pulseChargeInPMTs1; // pulse charge
            pulseChargeInPMTs1.resize(pmtnb,std::vector<double>(0,0)); 

            double topcharge_mu=0, topcharge_e=0;
            double botcharge_mu=0, botcharge_e=0;
            double totalcharge_mu=0, totalcharge_e=0;
            int nhit_topPMT=0, nhit_botPMT=0;
            double top_bot_ratio_eSig=0, top_bot_ratio_muSig=0;
            bool hitflag_topPMT=false, hitflag_botPMT=false;

            //cout<<ffff<<"\t"<<i<<"\t"<<mc->GetMCPMTCount()<<endl;
            double thispmt_charge[8], thispmt_charge_50ns[8], thispmt_charge_40ns[8];
            for(int t=0; t<8; t++){ thispmt_charge[t]=0; thispmt_charge_50ns[t]=0; thispmt_charge_40ns[t]=0;}
            //in simulation, if a pmt has no hit it will not record, so mc->GetMCPMTCount()<=8.
            for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++) { // for each pmt
                RAT::DS::MCPMT *mcpmt = mc->GetMCPMT(imcpmt);
                int pmtID = mcpmt->GetID(); // pmt id
                int nPhotons = mcpmt->GetMCPhotonCount(); // number of photons in this pmt
                int firstphotonflag = 0;
		double firstphotontime = 0;
		//See if there's any no-photon-hit PMT in this event.
                //if(nPhotons<1){cout<<"Evt "<<i<<", pmt "<<pmtID<<", no photons."<<endl; continue;}
                double thispmt_charge_afterdecay=0;
                double iTime_firstPhoton = 0;
                for (int iphoton=0; iphoton < nPhotons; iphoton++)  { // for all the photons produced in this pmt
                    double iTime = mcpmt->GetMCPhoton(iphoton)->GetHitTime()-evtime; // global (absolute) time of the photon
                    if(iphoton==0) iTime_firstPhoton = iTime;
                    double iCharge = mcpmt->GetMCPhoton(iphoton)->GetCharge()/converterfactor[pmtID]; // charge of this photon (applied my PMT SPE spectrum)
                    //don't consider photons that arrive later than 2500 ns, 
		    //also I set a 0.5 pe threshold, even if it is a photon, but if its charge is <0.5 pe, don't consider
                    //the charge could be <1 pe because the PMT SPE response is used
                    if(iTime>2500 || (iCharge<0.5 && iphoton==0)) {continue;} 
                    pulseTimeInPMTs1[pmtID].push_back(iTime);
                    pulseChargeInPMTs1[pmtID].push_back(iCharge);
		    totQ1 += iCharge;
    	            thispmt_charge[pmtID] += iCharge;
    	            if(iTime<=(50+iTime_firstPhoton)){
    	                thispmt_charge_50ns[pmtID]+=iCharge;
    	            }
    	            if(iTime<=(40+iTime_firstPhoton)){
    	                thispmt_charge_40ns[pmtID]+=iCharge;
    	            }
	            for(int trigindex=0; trigindex<nbTrigs; trigindex++){
			if(istrigger[trigindex]==true){
			    hPMTPhotonTimeSpan[pmtID][trigindex]->Fill(iTime);
			}
		    }
                    if(pmtID==4||pmtID==5) { // in top PMTs
			totQ_top1 += iCharge;
                    }
                    else{ // in bottom PMTs
                        totQ_bot1 += iCharge;
                    }
                    firstphotonflag++;
                    if(firstphotonflag==1) {
			firstphotontime = iTime;
                        for(int trigindex=0; trigindex<nbTrigs; trigindex++){
                            if(istrigger[trigindex]==true ){ // first photon in this PMT in this event 
                                hTime[pmtID][trigindex]->Fill(iTime);
                            }
                        }
                    }
                }//end of getting photons in one pmt
		for(int trigindex=0; trigindex<nbTrigs; trigindex++){
		    if(istrigger[trigindex]==true){
			if(firstphotonflag==1){
			    hPMTNpeVsTime[pmtID][trigindex]->Fill(thispmt_charge[pmtID], firstphotontime);							
			}
		    }
		}
            }//end of getting photons from all PMTs
            for(int trigindex=0; trigindex<nbTrigs; trigindex++){
                for(int pmtID=0; pmtID<8; pmtID++){
                    if(istrigger[trigindex]==true){
                        hCharge[pmtID][trigindex]->Fill(thispmt_charge[pmtID]);
                        hCharge1[pmtID][trigindex]->Fill(thispmt_charge_50ns[pmtID]*npescalefactor[pmtID]);
                        hCharge2[pmtID][trigindex]->Fill(thispmt_charge_40ns[pmtID]);
                    }
                }
	    }
	    for(int trigindex=0; trigindex<nbTrigs; trigindex++){
		if(istrigger[trigindex]==true){
		    hTotalCharge[trigindex]->Fill(totQ1);
		    hTotalChargeTop[trigindex]->Fill(totQ_top1);
	            hTotalChargeBottom[trigindex]->Fill(totQ_bot1);
		    if(totQ1>0){
			hNpeRatio[trigindex]->Fill(totQ_top1/totQ1*(2/8)); // charge ratio, weighted by number of PMTs.
	            }
		}
	    }
            // get the time difference between pmts, for all photons
	    for(int trigindex=0; trigindex<nbTrigs; trigindex++){
                for(int ngroup = 0; ngroup<7; ngroup++){
                    for(int nbpulseS0=0; nbpulseS0<pulseTimeInPMTs1[ngroup].size(); nbpulseS0++){
                        for(int pmtid=ngroup+1;pmtid<8;pmtid++){
                            for(int nbpulseS1=0; nbpulseS1<pulseTimeInPMTs1[pmtid].size(); nbpulseS1++){
                                double dT = (pulseTimeInPMTs1[pmtid][nbpulseS1]-meanPulseStartBin[pmtid])
            								- (pulseTimeInPMTs1[ngroup][nbpulseS0]-meanPulseStartBin[ngroup]);
                                if(ngroup==0) hPulseTimeDiffToS0[pmtid-1-ngroup][trigindex]->Fill( dT );
                                if(ngroup==1) hPulseTimeDiffToS1[pmtid-1-ngroup][trigindex]->Fill( dT );
                                if(ngroup==2) hPulseTimeDiffToS2[pmtid-1-ngroup][trigindex]->Fill( dT );
                                if(ngroup==3) hPulseTimeDiffToS3[pmtid-1-ngroup][trigindex]->Fill( dT );
                                if(ngroup==4) hPulseTimeDiffToS4[pmtid-1-ngroup][trigindex]->Fill( dT );
                                if(ngroup==5) hPulseTimeDiffToS5[pmtid-1-ngroup][trigindex]->Fill( dT );
                                if(ngroup==6) hPulseTimeDiffToS6[pmtid-1-ngroup][trigindex]->Fill( dT );
                            }
                        }
                    }
                }
	    }
            //reset the flags
            for(int trigindex=0; trigindex<nbTrigs; trigindex++){
                istrigger[trigindex]=false;
                eventStartTime[trigindex]=0;
	    }				
        }  // end of for loop reading ONE root file
        f->Close(); // close current root tree file
    }// end of loop for reading ALL root files

	//Now beging writing histograms in the output root file
    foutroot->cd();
    std::cout<<"Done reading files, now writting to ROOT file ..."<<std::endl;
	dir[0]->cd();
	//plastic scintillators (hodo detectors)
	TCanvas* cEdepInPS = new TCanvas("cEdepInPS","",1200,900);
	TLegend* legEdepInPS = new TLegend(0.6,0.3,0.9,0.5);
	TCanvas* cHitTimeInPS = new TCanvas("cHitTimeInPS","",1200,900);
	TLegend* legHitTimeInPS = new TLegend(0.6,0.3,0.9,0.5);
    for(int i=0; i<6; i++){
        for(int j=0; j<nbTrigs; j++){
            dir[j]->cd();
			cEdepInPS->cd(); 
			hEdepInPS[i][j]->Draw("sames");
			hEdepInPS[i][j]->SetLineColor(i+1);
			legEdepInPS->AddEntry(hEdepInPS[i][j],hEdepInPS[i][j]->GetName(),"l");
			
			cHitTimeInPS->cd(); 
			hHitTimeInPS[i][j]->Draw("sames");
			hHitTimeInPS[i][j]->SetLineColor(i+1);
			legHitTimeInPS->AddEntry(hHitTimeInPS[i][j],hHitTimeInPS[i][j]->GetName(),"l");
			
            hEdepInPS[i][j]->Write();
            hHitTimeInPS[i][j]->Write();
        }
    }
	cEdepInPS->cd(); legEdepInPS->Draw();  cEdepInPS->Write(); cEdepInPS->Close();
	cHitTimeInPS->cd(); legHitTimeInPS->Draw(); cHitTimeInPS->Write(); cHitTimeInPS->Close();
	
	hTDC_H0H2->Write();
	hTDC_H3H2->Write();
	hTDC_H1H2->Write();
	hTDC_H4H2->Write();
	hTDC_H5H2->Write();
	hTDC_H3H0->Write();
	hTDC_H1H0->Write();
	hTDC_H4H0->Write();
	hTDC_H5H0->Write();
	hTDC_H4H5->Write();
	hEdep_vs_TDC_H1H2->Write();
	hEdepH4_vs_EdepH5->Write();

	//input muon information
	TCanvas* cMuonInfo = new TCanvas("cMuonInfo","",1200,900);
	cMuonInfo->Divide(2,1);
	cMuonInfo->cd(2);
	TLegend* legMuonAngles = new TLegend(0.6,0.3,0.9,0.5);
	for(int j=0; j<nbTrigs; j++){
            dir[j]->cd();
	    cMuonInfo->cd(1); 
	    hMuMomentum[j]->Draw();
	    hMuMomentum[j]->Write();
		
	    cMuonInfo->cd(2);
	    hMuPAngle[j]->Draw("sames");
	    hMuPAngle[j]->SetLineColor(1);
	    hMuAAngle[j]->Draw("same");
	    hMuAAngle[j]->SetLineColor(2);
	    legMuonAngles->AddEntry(hMuPAngle[j],hMuPAngle[j]->GetName(),"l");
	    legMuonAngles->AddEntry(hMuAAngle[j],hMuAAngle[j]->GetName(),"l");
            hMuPAngle[j]->Write();
            hMuAAngle[j]->Write();
	    hEvtTime[j]->Write(); // the time defined by hodoDetectors
	    hTrigTime[j]->Write();
	}
	cMuonInfo->cd(2); legMuonAngles->Draw(); cMuonInfo->Write();cMuonInfo->Close();
	

    for(int i=0; i<pmtnb; i++){
        for(int j=0; j<nbTrigs; j++){
            dir[j]->cd();
	    if(i>0) hPulseTimeDiffToS0[i-1][j]->Write();
            if(i>1) hPulseTimeDiffToS1[i-2][j]->Write();
            if(i>2) hPulseTimeDiffToS2[i-3][j]->Write();
            if(i>3) hPulseTimeDiffToS3[i-4][j]->Write();
            if(i>4) hPulseTimeDiffToS4[i-5][j]->Write();
            if(i>5) hPulseTimeDiffToS5[i-6][j]->Write();
            if(i>6) hPulseTimeDiffToS6[i-7][j]->Write();
        }
    }

    for(int i=0; i<pmtnb; i++){
        for(int j=0; j<nbTrigs; j++){
            dir[j]->cd();
            hCharge[i][j]->Write();
            hCharge1[i][j]->Write();
            hCharge2[i][j]->Write();
	    hTime[i][j]->Write();
            hPMTNpeVsTime[i][j]->Write();
            hPMTPhotonTimeSpan[i][j]->Write();
        }
    }
	
	for(int j=0; j<nbTrigs; j++){
	    dir[j]->cd();
	    hTotalCharge[j]->Write();
	    hTotalChargeTop[j]->Write();
	    hTotalChargeBottom[j]->Write();
	    hNpeRatio[j]->Write();	
            hDecayZ[j]->Write();
            hEndXY[j]->Write();
            hEndXY0[j]->Write();
	}

    std::cout<<"Done reading files, now writting to ROOT file ..."<<std::endl;

    cout<<"Trigger type: \tnumber of events\n";
    for(int trigindex=0; trigindex<nbTrigs; trigindex++){
      cout<<strig[trigindex]<<"\t"<<trigcnts[trigindex]<<endl;
      gNumberOfEvents->SetPoint(trigindex,std::strtol(strig[trigindex].c_str(), 0, 2),trigcnts[trigindex]);
    }
    gNumberOfEvents->GetXaxis()->SetTitle("TrigType (converted from binary values)");
    gNumberOfEvents->GetYaxis()->SetTitle("Number Of Events");
    gNumberOfEvents->SetName("gNumberOfEvents_allTrigCombinations");
    TCanvas* cNumberOfEvents = new TCanvas();
    cNumberOfEvents->SetName("cNumberOfEvents");
    gNumberOfEvents->Draw("AP"); gNumberOfEvents->SetMarkerStyle(8);
    foutroot->cd();
    cNumberOfEvents->Write();
    cNumberOfEvents->Close();
	
    cout<<"HodoTrig "<<trigcnts_hodo<<endl;
    cout<<"Hodo + H5 + notH4 "<<trigcnts_hodo_H5_notH4<<endl;
    cout<<"Hodo + notH5 + notH4 "<<trigcnts_hodo_notH5_notH4<<endl;
    cout<<"Hodo + H5 + H4 "<<trigcnts_hodo_H5_H4<<endl;
    cout<<"Hodo + notH5 + H4 "<<trigcnts_hodo_notH5_H4<<endl;
    cout<<" - in which number of decays "<<trigcnts_decays<<endl;


    cout<<"total_number_of_events = "<<total_number_of_events<<endl;

    cout<<"testCntType1 ---- "<<endl;
    for(int t=0;t<8;t++) cout<<testcntTtrigtype1[t]<<"\t";
    cout<<endl;
    foutroot->Close();
    
    return 0; 

}
