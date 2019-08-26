#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TTimeStamp.h>
#include <TDatime.h>
#include <TMath.h>

#include <bits/stdc++.h> 
#include <stdlib.h>  

//using namespace std;
using namespace std;
int main(int argc, char** argv){
    // prepare for filenames to be processed.
    // the logic here is to read ROOT tree files and plot histograms
    // more trees can be read at same time, but only one ROOT file will be saved (histograms get higher statistics).
    // usage:
    //  -- compile this script: g++ -I $ROOTSYS/include -L $ROOTSYS/lib `root-config --libs` Read_Tree_makePlots_*.c -o Read_Tree_makePlots_*
    //  -- check that a filename list exists, for example, filenamelist.txt, in the list, ROOT files are listed, eg., run28656_Tree_v7.root
    //  -- run the program: ./Read_Tree_makePlots_V* NbOfFiles filenamelist.txt
    if(argc==1){
      cout<<"############################################################"<<endl;
      cout<<"# Usage: ./Read_Tree_makePlots_V* NumberOfFiles FileNameList"<<endl;
      cout<<"############################################################"<<endl;
      return 0;
    }
    // ----------------------------------------------------
    // read in the filename list, store the names in array.
    // ----------------------------------------------------
    const int NbOfFilesToRead = atoi(argv[1]);
    std::string filepath = "./";
    std::string rootfilename[NbOfFilesToRead];
    std::string temp_filename = "./";
    temp_filename += argv[2];
    fstream fin_filenamelist(temp_filename.c_str(),ios::in);
    std::string aname="";
    int index = 0; 
    while(fin_filenamelist>>aname){
        rootfilename[index] = aname;
        index++;
        if(index>=NbOfFilesToRead) break;
    }
    fin_filenamelist.close();
    cout<<"In "<<temp_filename<<", "<<index<<" runs ..."<<endl;
    if(index!=NbOfFilesToRead){cout<<"There seems to be a mis-match of the number of files"<<endl; return 0;}
    // ------------------------------
    // prepare a root file for output
    // ------------------------------
    std::string strhelp1 = "./ResultsV16_allTrigCombinations_";
	std::string runBeginStr = rootfilename[0].substr(rootfilename[0].length()-21,8);
	std::string runEndStr = rootfilename[NbOfFilesToRead-1].substr(rootfilename[NbOfFilesToRead-1].length()-21,8);
    std::string outputrootfilename = strhelp1 + runBeginStr + "_" + runEndStr + "_20190605.root"; 

    // -----------------------------
	// the input file containing PMT calibrations
	fstream fin_PMTCalibration("Oneton_PMTCalibration_forInput_20190122.txt",ios::in);
	vector<int> input_runNumberBegin;
	vector<int> input_runNumberEnd;
	vector<vector<double> > input_PMTcalibration;
	int cnt=0;
	double tempvalue=0;
	while(fin_PMTCalibration>>aname){
		std::string temp_run_begin_str = aname.substr(aname.length()-28,5);
		std::string temp_run_end_str = aname.substr(aname.length()-19,5);
		input_runNumberBegin.push_back(atoi(temp_run_begin_str.c_str()));
		input_runNumberEnd.push_back(atoi(temp_run_end_str.c_str()));
		cnt++;
		input_PMTcalibration.resize(cnt);
		for(int i=0;i<8;i++){
			fin_PMTCalibration>>tempvalue;
			input_PMTcalibration[cnt-1].push_back(tempvalue);
		}
	}
	fin_PMTCalibration.close();

    // some mapping inforamtion
    // trigger types: all, multiplicity only, cosmic only, led only, stopped muon candidates, and cosmic triggers after event selection
    // For the moment, I use one trigger type at a single run. (There's problem to run multiple triggers because of memory limit on cluster)
    const int pmtnb = 8;
    const int nbTrig = 27;
    // trigger types:
    // I have 27 types of combinations based on the normal hodotrig defined by (H2 or H0) and (H3 or H1)
    // Using the order H2, H0, H3, H1, H4, H5
    // and 0, 1 for hit or no hit, I list the combinations as following:
    string sTriggerType[nbTrig]={"101000", "101010", "101001",
                           "100100", "100110", "100101",
                           "011000", "011010", "011001",
                           "010100", "010110", "010101",
                           "111000", "111010", "111001",
                           "110100", "110110", "110101",
                           "011100", "011110", "011101",
                           "101100", "101110", "101101",
                           "111100", "111110", "111101"
                          };
    int indexMap_forEvtStartTime[nbTrig]={2,2,2,2,2,2,  0,0,0,0,0,0,   2,2,2,2,2,2, 0,0,0,  2,2,2,2,2,2};

    bool istrigger[nbTrig]; 
    int  triggercnt[nbTrig];
    double eventStartTime[nbTrig];
    for(int index=0; index<nbTrig; index++){istrigger[index]=false; eventStartTime[index]=0; triggercnt[index]=0;} 

    int nhit = 0; // the number of hitting pmts in each event

    std::string sTDC_ch_Map[16]={"H0","H1","H2","H3","H4","H5","CosmicTrig","MultTrig","LedTrig","AllTrig","S0","S1","S2","S3","S4","S5"}; // here the "AllTrig" maybe not really all, but multi_or_hodo_or_led_, I will clarify later
    std::string sQDC1ChMap[8]={"H0","H1","H2","H3","H4","H5","",""};
    std::string sQDC2ChMap[8]={"S0","","S2","S3","S4","S5","S1",""};
    std::string sDigiChMap[8]={"Digi1_CH0_S0","Digi1_CH1_S1","Digi1_CH2_S2","Digi1_CH3_S3","Digi2_CH0_S4","Digi2_CH1_S5","Digi2_CH2_S6","Digi2_CH3_S7"};
    std::string sDigiCh[8] = {"S0","S1","S2","S3","S4","S5","S6","S7"};
    int TDCIndexMap[6] = {10,11,12,13,14,15}; // S0 - TDC10; S1 - TDC11; ..., TDC for S6 and S7 are not avaialable

    double converterfactor[8]={153.345, 173.944, 187.946, 182.012, 176.795, 178.761, 165.962, 148.923}; // spe charge measured by the PMTs
    double converterTDCSlope=(1200.0-140.0)/3840.0; double converterTDCInte=140.0;
    float meanPulseStartBin[8]={150.38, 150.05, 168.97, 151.27, 126.63, 125.17, 117.32, 114.80};
    //float meanPulseStartBin[8]={106.386, 106.028, 124.967, 107.329, 82.5967, 81.1806, 73.2965, 70.7959}; // the mean values of the first pulses' start time in every PMT
    float upperBoundFirstPulseStartBin[8]={200,190,215,190,170,170,155,155};// the upper bound of the start time bin of the first pulse in the PMTs, from my new hodoTrigger
   // the mean values of the muon time between PMTs. since I have 28 groups between PMTs, the upper part I use -100 to indicate not used.
    double coincidenceTimeMeanMatrix[8][8]={
              {-100,      -100,       -100,      -100,      -100,     -100,      -100,      -100},
              {0.0739848, -100,       -100,      -100,      -100,     -100,      -100,      -100},
              {0.0527475, 0.530554,   -100,      -100,      -100,     -100,      -100,      -100},
              {0.43977,   -0.0985247, -0.12697,  -100,      -100,     -100,      -100,      -100},
              {-0.159319, 0.318787,   0.458081,  0.0187623, -100,     -100,      -100,      -100},
              {0.447696,  -0.0542085, 0.0158028, 0.58527,   0.229519, -100,      -100,      -100},
              {0.922964,  0.343448,   0.0859549, 0.814555,  0.42671,  -0.118653, -100,      -100},
              {0.296518,  0.771128,   0.494532,  0.26148,   -0.158229, 0.324322, -0.0807268, -100}    
                                       };
    double coincidenceTimeWidthMatrix[8][8]={
              {-100,      -100,       -100,      -100,      -100,     -100,      -100,      -100},
              {1.2153200, -100,       -100,      -100,      -100,     -100,      -100,      -100},
              {1.6912900, 1.408720,   -100,      -100,      -100,     -100,      -100,      -100},
              {1.67040,   1.51714000, 1.360990,  -100,      -100,     -100,      -100,      -100},
              {5.4553000, 5.379390,   5.348310,  5.2695700, -100,     -100,      -100,      -100},
              {5.130830,  5.02824000, 5.0451100, 4.98300,   5.406260, -100,      -100,      -100},
              {2.775520,  2.744940,   3.0918200, 3.097860,  4.52625,  4.1402900, -100,      -100},
              {2.835060,  2.696970,   2.996930,  3.018390,  4.486300, 4.1450400, 0.7568430, -100}    
                                       };

    // root file to save histograms
    TFile* rtfile = new TFile(outputrootfilename.c_str(),"recreate");
    // directories - environmental
    TDirectory* dirEnv = rtfile->mkdir("Environment");
    TDirectory* dirPedestal = rtfile->mkdir("Pedestals");
    TDirectory* dirEvtRate = rtfile->mkdir("EvtRate");
    TDirectory* dirRebuildPulse = rtfile->mkdir("RebuiltPulse");
    // directories - by trigger type
    TDirectory* dir_TrigType[nbTrig];
    TDirectory* dir_QDC[nbTrig];
    TDirectory* dir_TDC[nbTrig];
    TDirectory* dir_Digitizer[nbTrig];
    TDirectory* dir_evtTimeDiff[nbTrig];
    char dirnames[100];
    for(int i=0; i<nbTrig; i++){
        sprintf(dirnames,"dir_%s",sTriggerType[i].c_str());
        dir_TrigType[i] = rtfile->mkdir(dirnames);
        //dir_TrigType[i]->cd();
        // and sub-directories
        //dir_QDC[i] = dir[i]->mkdir("QDC");
        //dir_TDC[i] = dir[i]->mkdir("TDC");
       // dir_Digitizer[i] = dir[i]->mkdir("Digitizer");
       // dir_evtTimeDiff[i] = rtfile->mkdir("TimeDiff");
    }

    // -------------------------------------
    // define the histograms
    // -------------------------------------
    //TH1F* hEvtTimeDiff[nbTrig]; // trigger timestampe difference, unit in second. I confirmed it is expoentially distributed, not very useful in my analysis
    // QDC
    //TH1F* hQDC1[8][nbTrig]; // nbTrig is defined in the beginning, for different triggers: all, multi, hodo, led, m_or_h, m_and_h 
    //TH1F* hQDC2[8][nbTrig];
    // digitizer - pedestals and threshold for every PMT, get their statistics. Not very useful, but as a reference
    TH1F* hPedestalMean[8];
    TH1F* hPedestalWidth[8];
    //TH1F* hWaveformThreshold[8];
    // digitizer
    // -- amplitude: it is global, every waveform will have only one minimum (if no pulse the first appreance of minimum is taken)
    //TH1F* hDigiAmp[8][nbTrig]; // the amplitude information of the pulses 
    //TH1F* hDigiAmpTimeBin[8][nbTrig]; // the time bin of the amplitude
    //TH2F* hScat_DigiAmp_vs_TimeBin[8][nbTrig]; // amplitude vs time bin
    //TH2F* hScat_DigiAmpTS0_vs_DigiAmpTS1[nbTrig]; // time bin of the amplitude correlation between S0 and S1
    // -- pulses: the first pulse is the most important
    TH1F* hDigiST[pmtnb][nbTrig];  // start time bin of the first pulse
    //TH1F* hDigiST_all[8][nbTrig]; // start time bin of EVERY pulse
    //TH1F* hDigiPW[8][nbTrig];  // pulse width of the first pulse
    //TH1F* hNumberOfPulses[8][nbTrig]; // the number of pulses 
    //TH1F* hDigitizerPulseCharge_npe[8][nbTrig]; // total charge of all found pulses in npe
    //TH2F* hDigiPulseChargeNpe_vs_PulseWidth[8][nbTrig]; // charge of first pulse in npe vs. width of first pulse
    TH1F* hDigiPulseAmp[pmtnb][nbTrig]; // amplitude of every found pulse in each PMT
    //TH1F* hDigitizerPulseCharge1[8][nbTrig];// first pulse's charge in npe
    //TH1F* hDigitizerChargeRatio[8][nbTrig];// charge ratio of the first pulse to all pulses
    //TH2F* hDigitizerPulseCharge_vs_PulseTime[8][nbTrig]; // pulse charge in npe vs. pulse start time bin, for every pulse
    //TH2F* hAfterPulseChargeCorrelation[8][nbTrig]; // the correlation between after charge pulse and the initial muon charge
    //TH2F* hDigitizerPulseWidth_vs_PulseTime[8][nbTrig];  // pulse width of every found pulse vs. its start time bin
    //TH1F* hDigitizerPulse_TimeDifference[8][nbTrig];  // the time difference between the first and last pulses (require at least two pulses are found)
    //TH1F* hStopMuonTime[8][nbTrig]; // the time where the stop muon candidate can give a pulse in each pmt.
    //TH2F* hDigiChargeNpe_vs_DigiAmp[8][nbTrig];
    //TH1F* hDigitizerPulseChargeSum[nbTrig]; // the total charge measured by all the PMTs, for different trigger types
    //TH1F* hDigitizerPulseChargeSumTop[nbTrig]; // the total charge measured by the top PMTs (S4 and S5), for different trigger types
    //TH1F* hDigitizerPulseChargeSumBottom[nbTrig]; // thee total charge measured by the bottom PMTs (S0->S3,S6,S7), for different trigger types
    //TH2F* hDigitizerPulseChargeSumCorrelation[nbTrig]; // the correlation between summed charge of top PMTs and summed charge of bottom PMTs
    //TH1F* hDigitizerPulseChargeSum1[nbTrig]; // the total charge measured by all the PMTs, for different trigger types, the first pulse
    //TH1F* hDigitizerPulseChargeSumTop1[nbTrig]; // the total charge measured by the top PMTs (S4 and S5), for different trigger types, the first pulse
    //TH1F* hDigitizerPulseChargeSumBottom1[nbTrig]; // the total charge measured by the bottom PMTs (S0->S3,S6,S7), for different trigger types, the first pulse
    //TH1F* hDigitizerPulseChargeSumS2S3[nbTrig]; // the charge sum of S2 and S3, to be compared with sum of S4 and S5, because they can be considered having similar geometry view.
    //TH1F* hNHit[nbTrig];
    //TH1F* hTimeDiff[nbTrig]; // the histogram of difference between pulse time and muon start time for every pulse in every PMT. This can tell the coincidence of decay muons
    //TH1F* hTimeDiff1[nbTrig]; // the histogram of difference between pulse time and muon start time for every pulse in every PMT.
    //TH1F* hTimeDiffSelect[nbTrig]; // the histogram of difference between pulse time and muon start time for every pulse in every PMT.
    //TH1F* hTimeDiffSelect1[nbTrig]; // the histogram of difference between pulse time and muon start time for every pulse in every PMT.
    //TH1F* hDigitizerPulseChargeDecay[nbTrig]; // the total charge (in npe) of the electron caused pulses from the identified muon decay candidates
    //TH1F* hDigitizerPulseChargeDecay1[nbTrig]; // the total charge (in npe) of the electron caused pulses from the identified muon decay candidates, with further cut
    //TH2F* hDecayChargeVsTime[nbTrig]; // the 2D histogram of total charge from decay electron pulses vs. the time
    //TH2F* hDecayChargeVsTime1[nbTrig]; // the 2D histogram of total charge from decay electron pulses vs. the time
    //TH2F* hChargeRatioVsTimeInDecay[nbTrig]; // In the decayed candidates: the ratio of charge in top 2 PMTs over the charge in bottom 6 PMTs, as a function of the time difference between time in top PMTs and time in bottom PMTs.
    //TH2F* hChargeRatioVsTimeInDecay1[nbTrig]; // with further cut
    //some useful scattering plots
    //TH2F* hScat_TDC_vs_PMTT[8][nbTrig]; // TDC entry of hodoTrig vs. PMT pulse start bin, this is for individual PMTs
    //TH2F* hScat_TDC_vs_PMTT1[8][nbTrig]; // TDC entry of PMT vs. PMT pulse start bin, this is for individual PMTs
    //TH1F* hTDC_minus_PMTT[8][nbTrig]; // TDC entry of PMT vs. PMT pulse start bin, this is for individual PMTs  
    //TH2F* hScat_hodoTDC_vs_PMTT[8]; // TDC entry of hodoscope trigger vs. PMT pulse start bin
    //TH2F* hScat_TDCS0_vs_TDCS1[nbTrig]; // TDC entresi in PMT S0 and PMT S1
    //TH2F* hScat_NPulses_vs_S0[7][nbTrig]; // number of pulses in other pmts vs. number of pulses in S0
    //TH1F* hTotalNumberOfPulses[nbTrig]; // total number of pulses in all the pmts, 
    //TH2F* hScat_PulseTime_vs_S0[7][nbTrig]; // pulse time of all pulses;
    //TH2F* hScat_PulseTime_vs_S1[6][nbTrig]; // pulse time of all pulses;
    //TH2F* hScat_PulseTime_vs_S2[5][nbTrig]; // pulse time of all pulses;
    //TH2F* hScat_PulseTime_vs_S3[4][nbTrig]; // pulse time of all pulses;
    //TH1F* hPulseTimeDiffToS0[7][nbTrig]; // the difference in pulse times between other PMTs and S0
    //TH1F* hPulseTimeDiffToS1[6][nbTrig]; // the difference in pulse times between other PMTs and S1
    //TH1F* hPulseTimeDiffToS2[5][nbTrig]; // the difference in pulse times between other PMTs and S2
    //TH1F* hPulseTimeDiffToS3[4][nbTrig]; // the difference in pulse times between other PMTs and S3
    //TH1F* hPulseTimeDiffToS4[3][nbTrig]; // the difference in pulse times between other PMTs and S4
    //TH1F* hPulseTimeDiffToS5[2][nbTrig]; // the difference in pulse times between other PMTs and S5
    //TH1F* hPulseTimeDiffToS6[1][nbTrig]; // the difference in pulse times between other PMTs and S6
    //TH1F* hCoincidence[nbTrig]; // histogram to indicate how many pulses are in coincidencce in time
    //TH2F* hCoincidenceVsDecayTime[nbTrig]; // 
    //TH2F* hCoincidenceNumbers[nbTrig];

    int selectHodo = 0, selectHs = 0;
    int testcnt1=0, testcnt2=0;
    int selectHodoTrig=0, selectHsTrig=0;
    //TH2F* hTestHodoSelection = new TH2F("hTestHodoSelection","",4,0,2,4,0,2);
    //TH2F* hTDChodo_vs_H0 = new TH2F("TDChodo_vs_TDCH0","",2500,-2,4998,2500,-2,4998);
    //TH2F* hTDChodo_vs_H1 = new TH2F("TDChodo_vs_TDCH1","",2500,-2,4998,2500,-2,4998);
    //TH2F* hTDChodo_vs_H2 = new TH2F("TDChodo_vs_TDCH2","",2500,-2,4998,2500,-2,4998);
    //TH2F* hTDChodo_vs_H3 = new TH2F("TDChodo_vs_TDCH3","",2500,-2,4998,2500,-2,4998);
    //TH2F* hTDChodo_vs_H4 = new TH2F("TDChodo_vs_TDCH4","",2500,-2,4998,2500,-2,4998);
    //TH2F* hTDChodo_vs_H5 = new TH2F("TDChodo_vs_TDCH5","",2500,-2,4998,2500,-2,4998);
    //TH2F* hTDChodo_vs_multi = new TH2F("TDChodo_vs_TDCMulti","",2500,-2,4998,2500,-2,4998);
    //TH2F* hTDChodo_vs_Led = new TH2F("TDChodo_vs_TDCLed","",2500,-2,4998,2500,-2,4998);	

    TGraph* gMeanNPE[8]; // the mean npe vs. time for every pmts, for hodotriggers, used to monitor wbls stability

    char tempname[150], tempname1[150];
    // initialze the histograms
    for(int i=0;i<pmtnb;i++){
        //sprintf(tempname,"hTDCentryHodo_vs_waveformStartBin_%s", sDigiChMap[i].c_str());
        //hScat_hodoTDC_vs_PMTT[i] = new TH2F(tempname,"",1000,0,1000,500,-0.5,5000-0.5);
        //hScat_hodoTDC_vs_PMTT[i]->SetXTitle("Start bin of the first pulse");
        //hScat_hodoTDC_vs_PMTT[i]->SetYTitle("TDC hodo trigger ch. entry");

        sprintf(tempname,"hPedMean_%s",sDigiChMap[i].c_str());
        hPedestalMean[i] = new TH1F(tempname,"",100,-2,2);
        hPedestalMean[i]->SetXTitle("Pedestal mean (ADC channel)");
        hPedestalMean[i]->SetYTitle("Counts");

        sprintf(tempname,"hPedWidth_%s",sDigiChMap[i].c_str());
        hPedestalWidth[i] = new TH1F(tempname,"",150,0,3);
        hPedestalWidth[i]->SetXTitle("Pedestal width (ADC channel)");
        hPedestalWidth[i]->SetYTitle("Counts");

        //sprintf(tempname,"hWaveformThreshold_%s",sDigiChMap[i].c_str());
        //hWaveformThreshold[i] = new TH1F(tempname,"",300,0,15);
        //hWaveformThreshold[i]->SetXTitle("Waveform threshold (ADC channel)");
        //hWaveformThreshold[i]->SetYTitle("Counts");

        for(int j=0;j<nbTrig;j++){
            //QDC1 & QDC2 histogram
            //if(strcmp(sQDC1ChMap[i].c_str(),"")!=0){
            //    sprintf(tempname,"hQDC1_CH%d_%s_%s",i, sQDC1ChMap[i].c_str(),sTriggerType[j].c_str());
            //    hQDC1[i][j] = new TH1F(tempname,"",100,0,1000);
            //    hQDC1[i][j]->SetXTitle("QDC channel");
            //    hQDC1[i][j]->SetYTitle("Counts");
            //cout<<hQDC1[i][j]->GetName()<<endl;
            //}
            //if(strcmp(sQDC2ChMap[i].c_str(),"")!=0){
            //    sprintf(tempname,"hQDC2_CH%d_%s_%s",i,sQDC2ChMap[i].c_str(),sTriggerType[j].c_str());
            //    hQDC2[i][j] = new TH1F(tempname,"",100,0,1000);
            //    hQDC2[i][j]->SetXTitle("QDC channel");
            //    hQDC2[i][j]->SetYTitle("Counts");
            //}

            //histogram of event time difference,
            //if(i==0){
            //  sprintf(tempname,"hEvtTimeDiff_%s",sTriggerType[j].c_str());
            //  hEvtTimeDiff[j] = new TH1F(tempname,"",1000,0,300);
            //  hEvtTimeDiff[j]->SetXTitle("Time difference (s)");
            //  hEvtTimeDiff[j]->SetYTitle("Counts");
            //}

            // total charge of all PMTs
            //if(i==0){
              //sprintf(tempname,"hPulseChargeSum_npe_%s",sTriggerType[j].c_str());
              //hDigitizerPulseChargeSum[j] = new TH1F(tempname,"",1000,0,1000);
              //hDigitizerPulseChargeSum[j]->SetXTitle("Total measured n.p.e");
              //hDigitizerPulseChargeSum[j]->SetYTitle("Counts");
              //sprintf(tempname,"hPulseChargeSum_topPMTs_npe_%s",sTriggerType[j].c_str());
              //hDigitizerPulseChargeSumTop[j] = new TH1F(tempname,"",1000,0,1000);
              //hDigitizerPulseChargeSumTop[j]->SetXTitle("Total measured n.p.e");
              //hDigitizerPulseChargeSumTop[j]->SetYTitle("Counts");
              //sprintf(tempname,"hPulseChargeSum_botPMTs_npe_%s",sTriggerType[j].c_str());
              //hDigitizerPulseChargeSumBottom[j] = new TH1F(tempname,"",2000,0,500);
              //hDigitizerPulseChargeSumBottom[j]->SetXTitle("Total measured n.p.e");
              //hDigitizerPulseChargeSumBottom[j]->SetYTitle("Counts");
              //sprintf(tempname,"hPulseChargeSum_corrTB_npe_%s",sTriggerType[j].c_str());
              //hDigitizerPulseChargeSumCorrelation[j] = new TH2F(tempname,"",1000,0,1000,1000,0,1000);
              //hDigitizerPulseChargeSumCorrelation[j]->SetXTitle("Total measured n.p.e in bottom PMTs");
              //hDigitizerPulseChargeSumCorrelation[j]->SetYTitle("Total measured n.p.e in top PMTs");
              //sprintf(tempname,"hPulseChargeSum_firstPulse_npe_%s",sTriggerType[j].c_str());
              //hDigitizerPulseChargeSum1[j] = new TH1F(tempname,"",2000,0,500);
              //hDigitizerPulseChargeSum1[j]->SetXTitle("Total measured n.p.e");
              //hDigitizerPulseChargeSum1[j]->SetYTitle("Counts");
              //sprintf(tempname,"hPulseChargeSum_firstPulse_topPMTs_npe_%s",sTriggerType[j].c_str());
              //hDigitizerPulseChargeSumTop1[j] = new TH1F(tempname,"",2000,0,500);
              //hDigitizerPulseChargeSumTop1[j]->SetXTitle("Total measured n.p.e");
              //hDigitizerPulseChargeSumTop1[j]->SetYTitle("Counts");
              //sprintf(tempname,"hPulseChargeSum_firstPulse_botPMTs_npe_%s",sTriggerType[j].c_str());
              //hDigitizerPulseChargeSumBottom1[j] = new TH1F(tempname,"",2000,0,500);
              //hDigitizerPulseChargeSumBottom1[j]->SetXTitle("Total measured n.p.e");
              //hDigitizerPulseChargeSumBottom1[j]->SetYTitle("Counts");
              //sprintf(tempname,"hPulseChargeSum_firstPulse_S2S3_npe_%s",sTriggerType[j].c_str());
              //hDigitizerPulseChargeSumS2S3[j] = new TH1F(tempname,"",1000,0,1000);
              //hDigitizerPulseChargeSumS2S3[j]->SetXTitle("Total measured n.p.e");
              //hDigitizerPulseChargeSumS2S3[j]->SetYTitle("Counts");
              //sprintf(tempname,"hNHit_%s",sTriggerType[j].c_str());
              //hNHit[j] = new TH1F(tempname,"",9,0,9);
              //hNHit[j]->SetXTitle("Total measured n.p.e");
              //hNHit[j]->SetYTitle("Counts");
              //sprintf(tempname,"hTimeDiff_%s",sTriggerType[j].c_str());
              //hTimeDiff[j] = new TH1F(tempname,"",266,-100,2560);
              //hTimeDiff[j]->SetXTitle("Time difference between pulse time and muon start time");
              //hTimeDiff[j]->SetYTitle("Counts per 10 ns bins");
              //sprintf(tempname,"hTimeDiff_SomeSelection_%s",sTriggerType[j].c_str());
              //hTimeDiff1[j] = new TH1F(tempname,"",266,-100,2560);
              //hTimeDiff1[j]->SetXTitle("Time difference between pulse time and muon start time");
              //hTimeDiff1[j]->SetYTitle("Counts per 10 ns bins");
              //sprintf(tempname,"hTimeDiffSelect_%s",sTriggerType[j].c_str());
              //hTimeDiffSelect[j] = new TH1F(tempname,"",266,-100,2560);
              //hTimeDiffSelect[j]->SetXTitle("Time difference between pulse time and muon start time");
              //hTimeDiffSelect[j]->SetYTitle("Counts per 10 ns bins");
              //sprintf(tempname,"hTimeDiffSelect1_%s",sTriggerType[j].c_str());
              //hTimeDiffSelect1[j] = new TH1F(tempname,"",266,-100,2560);
              //hTimeDiffSelect1[j]->SetXTitle("Time difference between pulse time and muon start time");
              //hTimeDiffSelect1[j]->SetYTitle("Counts per 10 ns bins");
              //sprintf(tempname,"hCharge_decay_%s",sTriggerType[j].c_str());
              //hDigitizerPulseChargeDecay[j] = new TH1F(tempname,"",200,0,200);
              //hDigitizerPulseChargeDecay[j]->SetXTitle("Total charge in all PMTs (npe)");
              //hDigitizerPulseChargeDecay[j]->SetYTitle("Counts per 10 ns bins");
              //sprintf(tempname,"hCharge1_decay_%s",sTriggerType[j].c_str());
              //hDigitizerPulseChargeDecay1[j] = new TH1F(tempname,"",200,0,200);
              //hDigitizerPulseChargeDecay1[j]->SetXTitle("Total charge in all PMTs (npe)");
              //hDigitizerPulseChargeDecay1[j]->SetYTitle("Counts per 10 ns bins");
              //sprintf(tempname,"hChargeVsTime_decay_%s",sTriggerType[j].c_str());
              //hDecayChargeVsTime[j] = new TH2F(tempname,"",266,-100,2560,200,0,200);
              //hDecayChargeVsTime[j]->SetXTitle("Decay time (ns)");
              //hDecayChargeVsTime[j]->SetYTitle("Total charge in all PMTs (npe)");
              //sprintf(tempname,"hChargeVsTime1_decay_%s",sTriggerType[j].c_str());
              //hDecayChargeVsTime1[j] = new TH2F(tempname,"",266,-100,2560,200,0,200);
              //hDecayChargeVsTime1[j]->SetXTitle("Decay time (ns)");
              //hDecayChargeVsTime1[j]->SetYTitle("Total charge in all PMTs (npe)");
              //sprintf(tempname,"hChargeRatioVsTimeInDecay_%s",sTriggerType[j].c_str());
              //hChargeRatioVsTimeInDecay[j] = new TH2F(tempname,"",80,-20,20,240,-0.2,2.2);
              //hChargeRatioVsTimeInDecay[j]->SetXTitle("t_{top}-t_{bottom} (ns)");
              //hChargeRatioVsTimeInDecay[j]->SetYTitle("Ratio");
              //hChargeRatioVsTimeInDecay[i]->SetTitle("Ratio of charge between top and bottom PMTs, weighted by number of hitting PMTs");
              //sprintf(tempname,"hChargeRatioVsTimeInDecay1_%s",sTriggerType[j].c_str());
              //hChargeRatioVsTimeInDecay1[j] = new TH2F(tempname,"",80,-20,20,240,-0.2,2.2);
              //hChargeRatioVsTimeInDecay1[j]->SetXTitle("t_{top}-t_{bottom} (ns)");
              //hChargeRatioVsTimeInDecay1[j]->SetYTitle("Ratio");
              //hChargeRatioVsTimeInDecay1[i]->SetTitle("Ratio of charge between top and bottom PMTs, weighted by number of hitting PMTs");

              //sprintf(tempname,"hScat_AmpTimeBin_S0_vs_S1_%s_%s",sDigiChMap[i].c_str(),sTriggerType[j].c_str());
              //hScat_DigiAmpTS0_vs_DigiAmpTS1[j] = new TH2F(tempname,"",2560/2,0,2560,2560/2,0,2560);
              //hScat_DigiAmpTS0_vs_DigiAmpTS1[j]->SetXTitle("Time bin at minimum for S0");
              //hScat_DigiAmpTS0_vs_DigiAmpTS1[j]->SetYTitle("Time bin at minimum for S1");

              //sprintf(tempname,"hScat_TDCEntry_S0_vs_S1_%s_%s",sDigiChMap[i].c_str(),sTriggerType[j].c_str());
              //hScat_TDCS0_vs_TDCS1[j] = new TH2F(tempname,"",501,-10,5000,501,-10,5000);
              //hScat_TDCS0_vs_TDCS1[j]->SetXTitle("TDC for S0");
              //hScat_TDCS0_vs_TDCS1[j]->SetYTitle("TDC for S1");

              //sprintf(tempname,"hTotalNumberOfPulses_%s",sTriggerType[j].c_str());
              //hTotalNumberOfPulses[j] = new TH1F(tempname,"",20,0,20);
              //hTotalNumberOfPulses[j]->SetXTitle("Total number of pulses in all PMTs");
              //hTotalNumberOfPulses[j]->SetYTitle("Counts");
              //sprintf(tempname,"hCoincidence",sTriggerType[j].c_str());
              //hCoincidence[j] = new TH1F(tempname,"",100,0,100);
              //hCoincidence[j]->SetXTitle("Number of coincident pulses in each event");
              //hCoincidence[j]->SetYTitle("Number of events");
              //sprintf(tempname,"hCoincidenceVsDecayTime",sTriggerType[j].c_str());
              //hCoincidenceVsDecayTime[j] = new TH2F(tempname,"",256,0,2560,100,0,100);
              //hCoincidenceVsDecayTime[j]->SetXTitle("Decay time (ns)");
              //hCoincidenceVsDecayTime[j]->SetYTitle("Number of coincident pulses in each event");
              //sprintf(tempname,"hCoincidenceNumbers",sTriggerType[j].c_str());
              //hCoincidenceNumbers[j] = new TH2F(tempname,"",10,0,10,10,0,10);
              //hCoincidenceNumbers[j]->SetXTitle("Nb of pulses in the 1st coincidence");
              //hCoincidenceNumbers[j]->SetYTitle("Nb of pulses in the 2nd coincidence");
            //}
            /*
            if(i>0 && i<8){
              sprintf(tempname,"hPulseTimeDiffToS0_%s_%s",sDigiChMap[i].c_str(),sTriggerType[j].c_str());
              hPulseTimeDiffToS0[i-1][j] = new TH1F(tempname,"",2560*2,-2560,2560);
              sprintf(tempname,"Pulse time differencce between %s and S0", sDigiChMap[i].c_str());
              hPulseTimeDiffToS0[i-1][j]->SetXTitle(tempname);
              hPulseTimeDiffToS0[i-1][j]->SetYTitle("Counts / 1 ns");
              if(i>1){
                sprintf(tempname,"hPulseTimeDiffToS1_%s_%s",sDigiChMap[i].c_str(),sTriggerType[j].c_str());
                hPulseTimeDiffToS1[i-2][j] = new TH1F(tempname,"",2560*2,-2560,2560);
                sprintf(tempname,"Pulse time differencce between %s and S1", sDigiChMap[i].c_str());
                hPulseTimeDiffToS1[i-2][j]->SetXTitle(tempname);
                hPulseTimeDiffToS1[i-2][j]->SetYTitle("Counts / 1 ns");
              }
              if(i>2){
                sprintf(tempname,"hPulseTimeDiffToS2_%s_%s",sDigiChMap[i].c_str(),sTriggerType[j].c_str());
                hPulseTimeDiffToS2[i-3][j] = new TH1F(tempname,"",2560*2,-2560,2560);
                sprintf(tempname,"Pulse time differencce between %s and S2", sDigiChMap[i].c_str());
                hPulseTimeDiffToS2[i-3][j]->SetXTitle(tempname);
                hPulseTimeDiffToS2[i-3][j]->SetYTitle("Counts / 1 ns");
              }
              if(i>3){
                sprintf(tempname,"hPulseTimeDiffToS3_%s_%s",sDigiChMap[i].c_str(),sTriggerType[j].c_str());
                hPulseTimeDiffToS3[i-4][j] = new TH1F(tempname,"",2560*2,-2560,2560);
                sprintf(tempname,"Pulse time differencce between %s and S3", sDigiChMap[i].c_str());
                hPulseTimeDiffToS3[i-4][j]->SetXTitle(tempname);
                hPulseTimeDiffToS3[i-4][j]->SetYTitle("Counts / 1 ns");
              }
              if(i>4){
                sprintf(tempname,"hPulseTimeDiffToS4_%s_%s",sDigiChMap[i].c_str(),sTriggerType[j].c_str());
                hPulseTimeDiffToS4[i-5][j] = new TH1F(tempname,"",2560*2,-2560,2560);
                sprintf(tempname,"Pulse time differencce between %s and S4", sDigiChMap[i].c_str());
                hPulseTimeDiffToS4[i-5][j]->SetXTitle(tempname);
                hPulseTimeDiffToS4[i-5][j]->SetYTitle("Counts / 1 ns");
              }
              if(i>5){
                sprintf(tempname,"hPulseTimeDiffToS5_%s_%s",sDigiChMap[i].c_str(),sTriggerType[j].c_str());
                hPulseTimeDiffToS5[i-6][j] = new TH1F(tempname,"",2560*2,-2560,2560);
                sprintf(tempname,"Pulse time differencce between %s and S5", sDigiChMap[i].c_str());
                hPulseTimeDiffToS5[i-6][j]->SetXTitle(tempname);
                hPulseTimeDiffToS5[i-6][j]->SetYTitle("Counts / 1 ns");
              }
              if(i>6){
                sprintf(tempname,"hPulseTimeDiffToS6_%s_%s",sDigiChMap[i].c_str(),sTriggerType[j].c_str());
                hPulseTimeDiffToS6[i-7][j] = new TH1F(tempname,"",2560*2,-2560,2560);
                sprintf(tempname,"Pulse time differencce between %s and S6", sDigiChMap[i].c_str());
                hPulseTimeDiffToS6[i-7][j]->SetXTitle(tempname);
                hPulseTimeDiffToS6[i-7][j]->SetYTitle("Counts / 1 ns");
              }
            }
            */

            if(i<pmtnb){
                //for TDC, the 8 channels.
                //sprintf(tempname,"hTDCentryHodoTrig_vs_waveformStartBin_%s_%s",sDigiChMap[i].c_str(),sTriggerType[j].c_str());
                //hScat_TDC_vs_PMTT[i][j] = new TH2F(tempname,"",1000,0,1000,500,-0.5,5000-0.5);
                //hScat_TDC_vs_PMTT[i][j]->SetXTitle("Start bin of the first pulse");
                //hScat_TDC_vs_PMTT[i][j]->SetYTitle("TDC entry");
                //sprintf(tempname,"hTDCentryPMT_vs_waveformStartBin_%s_%s",sDigiChMap[i].c_str(),sTriggerType[j].c_str());
                //hScat_TDC_vs_PMTT1[i][j] = new TH2F(tempname,"",1000,0,1000,400,0,200);
                //hScat_TDC_vs_PMTT1[i][j]->SetXTitle("Start bin of the first pulse");
                //hScat_TDC_vs_PMTT1[i][j]->SetYTitle("TDC entry (in ns)");

                //sprintf(tempname,"hTimeDiff_%s_%s",sDigiChMap[i].c_str(),sTriggerType[j].c_str());
                //hTDC_minus_PMTT[i][j] = new TH1F(tempname,"",2700,-100,2600);
                //hTDC_minus_PMTT[i][j]->SetXTitle("PMT first pulse's time - TDC entry (in ns)");
                //hTDC_minus_PMTT[i][j]->SetYTitle("Counts");

                sprintf(tempname,"hStartTime_%s_%s",sDigiChMap[i].c_str(),sTriggerType[j].c_str());
                hDigiST[i][j] = new TH1F(tempname,"",2560,0,2560);
                hDigiST[i][j]->SetXTitle("Time bin number");
                hDigiST[i][j]->SetYTitle("Counts");

                //sprintf(tempname,"hStartTimeAllPulses_%s_%s",sDigiChMap[i].c_str(),sTriggerType[j].c_str());
                //hDigiST_all[i][j] = new TH1F(tempname,"",2560,0,2560);
                //hDigiST_all[i][j]->SetXTitle("Time bin number");
                //hDigiST_all[i][j]->SetYTitle("Counts");

                //sprintf(tempname,"hDigiPulseCharge_npe_%s_%s",sDigiChMap[i].c_str(),sTriggerType[j].c_str());
                //hDigitizerPulseCharge_npe[i][j] = new TH1F(tempname,"",50,0,50);
                //hDigitizerPulseCharge_npe[i][j]->SetXTitle("npe");
                //hDigitizerPulseCharge_npe[i][j]->SetYTitle("Counts");

                sprintf(tempname,"hDigiPulseAmplitude_%s_%s",sDigiChMap[i].c_str(),sTriggerType[j].c_str());
                hDigiPulseAmp[i][j] = new TH1F(tempname,"",5000,0,10000);
                hDigiPulseAmp[i][j]->SetXTitle("Amplitude of every fould pulse (ADC count)");
                hDigiPulseAmp[i][j]->SetYTitle("Counts");

                //sprintf(tempname,"hDigitizerPulseCharge_1stPulse_%s_%s",sDigiChMap[i].c_str(),sTriggerType[j].c_str());
                //hDigitizerPulseCharge1[i][j] = new TH1F(tempname,"",50,0,50);
                //hDigitizerPulseCharge1[i][j]->SetXTitle("charge in npe of first pulse");
                //hDigitizerPulseCharge1[i][j]->SetYTitle("Counts");

                //sprintf(tempname,"hDigitizerPulseCharge_vs_StartBin_%s_%s",sDigiChMap[i].c_str(),sTriggerType[j].c_str());
                //hDigitizerPulseCharge_vs_PulseTime[i][j] = new TH2F(tempname,"",520,-2,50,128,0,2560);
                //hDigitizerPulseCharge_vs_PulseTime[i][j]->SetXTitle("Charge (npe)");
                //hDigitizerPulseCharge_vs_PulseTime[i][j]->SetYTitle("Pulse start time bin number");

                //sprintf(tempname,"hAfterPulseChargeCorrelation_%s_%s",sDigiChMap[i].c_str(),sTriggerType[j].c_str());
                //hAfterPulseChargeCorrelation[i][j] = new TH2F(tempname,"",520,-2,50,520,-2,52);
                //hAfterPulseChargeCorrelation[i][j]->SetXTitle("Initial pulse charge (npe)");
                //hAfterPulseChargeCorrelation[i][j]->SetYTitle("Afterpulse charge (npe)");

                //sprintf(tempname,"hDigitizerPulse_TimeDifference_%s_%s",sDigiChMap[i].c_str(),sTriggerType[j].c_str());
                //hDigitizerPulse_TimeDifference[i][j] = new TH1F(tempname,"",128,0,2560);
                //hDigitizerPulse_TimeDifference[i][j]->SetXTitle("Pulse time diff. between two pulses (only for 2-pulse waveforms)");
                //hDigitizerPulse_TimeDifference[i][j]->SetYTitle("Counts");

                sprintf(tempname,"gMeanNPE_vs_elapsedTime_%s_%s",sDigiChMap[i].c_str(),sTriggerType[2].c_str());
                gMeanNPE[i]=new TGraph(); gMeanNPE[i]->SetName(tempname);
            }
        }
    }//end initialzing histogramgs

    //Histograms for plotting TDC
    TH1F* hTDC[16][nbTrig]; //histograms holding raw TDC data for each channel
    //TH1F* hTDC1[16]; //histograms holding raw TDC - TDC in ch9
    TH1F* hTDCnew[16][nbTrig]; // new histograms holding raw TDC data for each channel, the tirgger types are separated
    //TH1F* hTDCnewDelta[16][nbTrig]; // new histograms holding (TDC data - raw TDC data of trigger) for each channel, the tirgger types are separated
    for(int i=0;i<16;i++){
        char name[30];

        //sprintf(name,"hTDC_wrtCh9_%s",sTDC_ch_Map[i].c_str());
        //hTDC1[i]=new TH1F(name,name,4000,-400,400);
        //hTDC1[i]->SetXTitle("Time (ns)");        hTDC1[i]->SetYTitle("Counts");
        for(int j=0; j<nbTrig; j++){
            sprintf(name, "hTDCnew_%s_trigType_%s",sTDC_ch_Map[i].c_str(),sTriggerType[j].c_str());
            hTDCnew[i][j] = new TH1F(name,name,2500,-2,4998);
            hTDCnew[i][j]->SetXTitle("TDC");  hTDCnew[i][j]->SetYTitle("Counts"); 

	    sprintf(name,"hTDC_%s_%s",sTDC_ch_Map[i].c_str(),sTriggerType[j].c_str());
            hTDC[i][j]=new TH1F(name,name,210,-10,200);
            hTDC[i][j]->SetXTitle("TDC time in ns (1 TDC = 0.035 ns)");
	    hTDC[i][j]->SetYTitle("Counts");
            //sprintf(name, "hTDCRelative_%s_trigType_%s",sTDC_ch_Map[i].c_str(),sTriggerType[j].c_str());
            //hTDCnewDelta[i][j] = new TH1F(name,name,600,-1000.5,5000-0.5);
            //hTDCnewDelta[i][j]->SetXTitle("TDC");  hTDCnewDelta[i][j]->SetYTitle("Counts"); 
        }
    }
    //TH1F* hTDC_H0H2 = new TH1F("hTDC_H0_H2","TDC_H0 - TDC_H2", 2000,-2000,2000);
    //TH1F* hTDC_H1H2 = new TH1F("hTDC_H1_H2","TDC_H1 - TDC_H2", 2000,-2000,2000);
    //TH1F* hTDC_H3H2 = new TH1F("hTDC_H3_H2","TDC_H3 - TDC_H2", 2000,-2000,2000);
    //TH1F* hTDC_H1H0 = new TH1F("hTDC_H1_H0","TDC_H1 - TDC_H0", 2000,-2000,2000);
    //TH1F* hTDC_H3H0 = new TH1F("hTDC_H3_H0","TDC_H3 - TDC_H0", 2000,-2000,2000);
    //TH1F* hTDC_H5H2 = new TH1F("hTDC_H5_H2","TDC_H5 - TDC_H2", 2000,-2000,2000);
    //TH1F* hTDC_H4H2 = new TH1F("hTDC_H4_H2","TDC_H4 - TDC_H2", 2000,-2000,2000);
    //
    //TH1F* hTimeDiff_hodoDet_hodoTrig[6]; // TDC time difference betweeen hodoscope detectors and the hodoTrig
    //TH1F* hTimeDiff_PMTs_hodoTrig[8]; // TDC time converted to ns, PMT pulse start bin is in ns, look at their difference
    /*
    for(int i=0;i<6;i++){
        char name[100];
        sprintf(name,"hTimeDiff_H%d_vs_hodoTrig",i);
        hTimeDiff_hodoDet_hodoTrig[i]=new TH1F(name,name,400,-200,200);
        hTimeDiff_hodoDet_hodoTrig[i]->SetXTitle("Time diff. b/t hodoDet and hodoTrig in ns");
        hTimeDiff_hodoDet_hodoTrig[i]->SetYTitle("Counts");
    }		
    for(int i=0; i<8; i++){
        char name[100];
        sprintf(name,"hTimeDiff_S%d_vs_hodoTrig",i);
        hTimeDiff_PMTs_hodoTrig[i]=new TH1F(name,name,600,-300,300);
        hTimeDiff_PMTs_hodoTrig[i]->SetXTitle("Time diff b/t PMT and hodoTrig in ns");
        hTimeDiff_PMTs_hodoTrig[i]->SetYTitle("Counts");
    }
    */

    // histograms for plotting event rate (scaler)
    TH1F* hEvtRate[17]; // scaler data, in ROOT Tree it is represented as event rate
    // id = 7 -> hodo trigger, id = 8 -> multi trigger; id = 9 -> led trigger
    for(int i=0; i<17; i++){
      char name[30]; sprintf(name, "hEvtRate_ScalerChannel_%d",i);
      hEvtRate[i] = new TH1F(name,"",2500,0-0.5,2500-0.5);
      hEvtRate[i]->SetXTitle("Event rate (Hz)");
      hEvtRate[i]->SetYTitle("Counts");
    }
    TGraph* gEvtRateVsTime[3]; // event rate vs. time, I'll take the average event rate in each run and plot it vs. run time (time stamp of the last event in each run
    gEvtRateVsTime[0] = new TGraph(); gEvtRateVsTime[0]->SetName("gEvtVsTime_LedTrig"); 
    gEvtRateVsTime[1] = new TGraph(); gEvtRateVsTime[1]->SetName("gEvtVsTime_HodoTrig"); 
    gEvtRateVsTime[2] = new TGraph(); gEvtRateVsTime[2]->SetName("gEvtVsTime_MultiTrig"); 

    //Graphs for ploting temperature and resistivity
    TGraph* gTempInBox = new TGraph();
    TGraph* gTempOnRack = new TGraph();
    TGraph* gResistivity = new TGraph();
    // histograms for correlations between level, temperature, event rate
    TH2F* hLevelVsTemperature=new TH2F("hLevel1_vs_DarkboxT","hLevel1_vs_DarkboxT",200,15,35,500,0,5); // only use one level sensor data (from level sensor 1, inlet)
    TH2F* hLevelVsMultiEvtRate = new TH2F("hLevel1_vs_MultiEvtRate","hLevel1_vs_MultiEvtRate",100,0,50,500,0,5); // level vs. multiplicity evt rate
    //TH2F* hHodoRateVsMultiRate = new TH2F("hHodoRate_vs_MultiRate","hodoEvtRate vs multiEvtRate", 100,0,50,100,0,2);// hodo evt rate vs. multi evt rate

    // in the LED triggers, PMTs sometimes don't get a photon, for each PMT, the number of led triggers that give >0 photons can be computed
    // a ratio can be made between top PMTs and bottom PMTs. If material attenuation is changed, this ratio could change.
    // I made an 2*6 array of graphs, to compare the two top PMTs with the bottom 6 PMTs
    //TGraph* gLedEffRatioVsTime[2][6];// = new TGraph(); 
    int temp1index[2]={4,5};
    int temp2index[6]={0,1,2,3,6,7};
    //for(int i=0;i<2;i++){
    //  for(int j=0; j<6; j++){
    //    char name[100]; sprintf(name, "gLedEffRatioVsTime_S%d_vs_S%d",temp1index[i],temp2index[j]);
    //    gLedEffRatioVsTime[i][j]=new TGraph();
    //    gLedEffRatioVsTime[i][j]->SetName(name);
    //  }
    //}


    TH1F* hCharge[pmtnb][nbTrig]; // the total charge in every pmt within a time window of 2.5us (from event start time, 0 ns, in the MC: muon starts at Z=1200 mm).
    TH1F* hCharge1[pmtnb][nbTrig]; // the total charge in every pmt within a time window of 50ns. This might be useful
    for(int i=0;i<pmtnb;i++){
        for(int j=0;j<nbTrig;j++){
            sprintf(tempname,"hCharge_2500ns_S%d_%sTrig",i,sTriggerType[j].c_str());
            hCharge[i][j] = new TH1F(tempname,"",50,0,50);
            hCharge[i][j]->SetXTitle("Number of detected photons");
            hCharge[i][j]->SetYTitle("Counts");
            sprintf(tempname,"hCharge_50ns_S%d_%sTrig",i,sTriggerType[j].c_str());
            hCharge1[i][j] = new TH1F(tempname,"",50,0,50);
            hCharge1[i][j]->SetXTitle("npe in 50 ns window after the first photon");
            hCharge1[i][j]->SetYTitle("Counts");
        }
    }
    
    TGraphErrors* gLedEffRatioVsTimeTot = new TGraphErrors();  // total effect, 
    double evtLedTrigs_topPMT=0, evtLedTrigs_botPMT=0;
    double evtLedTrigs_total = 0;
    int Nmonitor=0;
    TDatime tmonitor;

    int totalEvtNb = 0;
    int totalEnvEvtNb = 0; // the counter for enviornment data, used for filling graphs of temperature, resistivity, etc.
    int trigcnts_all = 0;
    int trigcnts_multi=0;
    int trigcnts_hodo=0;
    int trigcnts_hodo_H5_H4=0;
    int trigcnts_hodo_H5_notH4=0;
    int trigcnts_hodo_notH5_H4=0;
    int trigcnts_hodo_notH5_notH4=0;
    int trigcnts_stricthodo = 0;
    int trigcnts_led=0;
    int trigcnts_stoppedmuon = 0;
    int trigcnts_multi_and_hodo = 0;
    int trigcnts_multi_or_hodo = 0;
    int trigcnts_crossmuonHodo = 0;
    int trigcnts_crossmuonStrictHodo = 0;
    int test_cnt1 = 0, test_cnt2 = 0;
    // starting reading ROOT tree files and fill histograms, one loop is for one root tree (one run)
    TDatime t00; // the time of the first event in the first run
 
    //try to put all definition outside loop
    TFile* rootfile;
    TTree* tree;
    // the branchses in that tree
    TBranch* bEvtNumber;// = tree->GetBranch("EvtNumber");
    TBranch* bEvtTime;// = tree->GetBranch("EvtTime");
    TBranch* bTemperatureBox;// = tree->GetBranch("TemperatureBox");
    TBranch* bTemperatureRack;// = tree->GetBranch("TemperatureRack");
    TBranch* bResistivity;// = tree->GetBranch("Resistivity");
    TBranch* bLevel1;// = tree->GetBranch("LiquidLevel1"); // liquid levels, needed later
    TBranch* bLevel2;// = tree->GetBranch("LiquidLevel2");
    TBranch* bQDC1;// = tree->GetBranch("QDC1");
    TBranch* bQDC2;// = tree->GetBranch("QDC2");

    TBranch* bTDC;// = tree->GetBranch("TDC");
    TBranch* bScaler;// = tree->GetBranch("Scaler");
    TBranch* bDigitizerCharge;// = tree->GetBranch("DigitizerCharge");
    TBranch* bDigitizerAmplitude;// = tree->GetBranch("DigitizerAmplitude");
    TBranch* bDigitizerAmplitudeTimeBin;// = tree->GetBranch("DigitizerAmplitudeTimeBin");
    TBranch* bDigitizerStartTimeBin;// = tree->GetBranch("DigitizerStartTimeBin");
    TBranch* bDigitizerRiseTime;// = tree->GetBranch("DigitizerRiseTime");
    TBranch* bDigitizerPulseWidth;// = tree->GetBranch("DigitizerPulseWidth");
    TBranch* bDigitizerPedMean;// = tree->GetBranch("DigitizerPedMean");
    TBranch* bDigitizerPedWidth;// = tree->GetBranch("DigitizerPedWidth");
    TBranch* bTrigType;// = tree->GetBranch("TrigType")   ;
    //new added branches
    TBranch* bTrigFlag_Multi;// = tree->GetBranch("TrigTypeFlag_Multi");
    TBranch* bTrigFlag_Hodo;// = tree->GetBranch("TrigTypeFlag_Hodo");
    TBranch* bTrigFlag_Led;// = tree->GetBranch("TrigTypeFlag_Led");
    TBranch* bDigitizerNumberOfPulses;// = tree->GetBranch("DigitizerNumberOfPulses");
    TBranch* bDigitizerPulseCharge;// = tree->GetBranch("DigitizerPulseCharge");
    TBranch* bDigitizerPulseStartBin;// = tree->GetBranch("DigitizerPulseStartBin");
    TBranch* bDigitizerPulseEndBin;// = tree->GetBranch("DigitizerPulseEndBin");
    TBranch* bDigitizerPulseAmplitude;// = tree->GetBranch("DigitizerPulseAmplitude");
    for(int aa = 0; aa<NbOfFilesToRead ; aa++){ //NbOfFilesToRead, begin to read one ROOT Treee
        std::string thisfilename = filepath + rootfilename[aa];// + postfix;
        cout<<thisfilename<<endl;
        int thisRunNumber = atoi(thisfilename.substr(thisfilename.length()-18,5).c_str());
        // apply the PMT calibrations
	for(int i=0; i<input_runNumberBegin.size();i++){
    	    if(thisRunNumber>=input_runNumberBegin[i] && thisRunNumber<=input_runNumberEnd[i]){
	        for(int j=0; j<input_PMTcalibration[i].size(); j++){
		    converterfactor[j]=input_PMTcalibration[i][j];
	        }
		cout<<"Calibration applied for run "<<thisRunNumber<<endl;
		for(int j=0; j<8; j++)
		    cout<<converterfactor[j]<<", ";
		cout<<endl;
		break;
	    }
	}

        rootfile = new TFile(thisfilename.c_str(),"read");
        tree = (TTree*)rootfile->Get("OneTonEvent");

        // the branchses in that tree
        bEvtNumber = tree->GetBranch("EvtNumber");
        bEvtTime = tree->GetBranch("EvtTime");
        bTemperatureBox = tree->GetBranch("TemperatureBox");
        bTemperatureRack = tree->GetBranch("TemperatureRack");
        bResistivity = tree->GetBranch("Resistivity");
        bLevel1 = tree->GetBranch("LiquidLevel1"); // liquid levels, needed later
        bLevel2 = tree->GetBranch("LiquidLevel2");
        bQDC1 = tree->GetBranch("QDC1");
        bQDC2 = tree->GetBranch("QDC2");
        bTDC = tree->GetBranch("TDC");
        bScaler = tree->GetBranch("Scaler");
        bDigitizerCharge = tree->GetBranch("DigitizerCharge");
        bDigitizerAmplitude = tree->GetBranch("DigitizerAmplitude");
        bDigitizerAmplitudeTimeBin = tree->GetBranch("DigitizerAmplitudeTimeBin");
        bDigitizerStartTimeBin = tree->GetBranch("DigitizerStartTimeBin");
        bDigitizerRiseTime = tree->GetBranch("DigitizerRiseTime");
        bDigitizerPulseWidth = tree->GetBranch("DigitizerPulseWidth");
        bDigitizerPedMean = tree->GetBranch("DigitizerPedMean");
        bDigitizerPedWidth = tree->GetBranch("DigitizerPedWidth");
        bTrigType = tree->GetBranch("TrigType")   ;
        //new added branches
        bTrigFlag_Multi = tree->GetBranch("TrigTypeFlag_Multi");
        bTrigFlag_Hodo = tree->GetBranch("TrigTypeFlag_Hodo");
        bTrigFlag_Led = tree->GetBranch("TrigTypeFlag_Led");
        bDigitizerNumberOfPulses = tree->GetBranch("DigitizerNumberOfPulses");
        bDigitizerPulseCharge = tree->GetBranch("DigitizerPulseCharge");
        bDigitizerPulseStartBin = tree->GetBranch("DigitizerPulseStartBin");
        bDigitizerPulseEndBin = tree->GetBranch("DigitizerPulseEndBin");
        bDigitizerPulseAmplitude = tree->GetBranch("DigitizerPulseAmplitude");
        // variables of the tree branches
        int evtNb;
        double evtT;
        float evtTempBox;
        float evtTempRack;
        float evtResist;
        float evtLevel1;
        float evtLevel2;
        float evtQDC1[8]     ;
        float evtQDC2[8];
        float evtTDC[16];
        float evtScaler[17];
        float evtDigiCharge[8];
        float evtDigiAmplitude[8];
        float evtDigiAmplitudeTimeBin[8]; 
        float evtDigiST[8]; //start time bin
        float evtDigiRT[8]; //rise time
        float evtDigiPW[8]; // pulse width
        float evtDigiPM[8]; // pedestal mean
        float evtDigiPS[8];//  pedestal width
        int evtTrigType;
        //new added branches
        int evtTrigFlag_Multi;
        int evtTrigFlag_Hodo;
        int evtTrigFlag_Led;
        int evtDigitizerNumberOfPulses[8];
        double evtDigitizerPulseCharge[8][20]; // evtDigitizerNumberOfPulses tells how many pulses, maximum 20 is allowed (defined)
        double evtDigitizerPulseStartBin[8][20];
        double evtDigitizerPulseEndBin[8][20];
        double evtDigitizerPulseAmplitude[8][20];

        bEvtNumber->SetAddress(&evtNb);
        bEvtTime->SetAddress(&evtT);
        bTemperatureBox->SetAddress(&evtTempBox);
        bTemperatureRack->SetAddress(&evtTempRack);
        bResistivity->SetAddress(&evtResist);
        bLevel1->SetAddress(&evtLevel1);
        bLevel2->SetAddress(&evtLevel2);
        bQDC1->SetAddress(&evtQDC1);
        bQDC2->SetAddress(&evtQDC2);
        bTDC->SetAddress(&evtTDC);
        bScaler->SetAddress(&evtScaler);
        bDigitizerCharge->SetAddress(&evtDigiCharge);
        bDigitizerAmplitude->SetAddress(&evtDigiAmplitude);
        bDigitizerAmplitudeTimeBin->SetAddress(&evtDigiAmplitudeTimeBin);
        bDigitizerStartTimeBin->SetAddress(&evtDigiST);
        bDigitizerRiseTime->SetAddress(&evtDigiRT);
        bDigitizerPulseWidth->SetAddress(&evtDigiPW);
        bDigitizerPedMean->SetAddress(&evtDigiPM);
        bDigitizerPedWidth->SetAddress(&evtDigiPS);
        bTrigType->SetAddress(&evtTrigType);
        // set newly added branches
        bTrigFlag_Multi->SetAddress(&evtTrigFlag_Multi);
        bTrigFlag_Hodo->SetAddress(&evtTrigFlag_Hodo);
        bTrigFlag_Led->SetAddress(&evtTrigFlag_Led);
        bDigitizerNumberOfPulses->SetAddress(&evtDigitizerNumberOfPulses);
        bDigitizerPulseCharge->SetAddress(&evtDigitizerPulseCharge);
        bDigitizerPulseStartBin->SetAddress(&evtDigitizerPulseStartBin);
        bDigitizerPulseEndBin->SetAddress(&evtDigitizerPulseEndBin);
        bDigitizerPulseAmplitude->SetAddress(&evtDigitizerPulseAmplitude);

        int nentries = (int)tree->GetEntries();
        double firstEvtTime = 0;
        double fEvtTime[nbTrig]; for(int itime=0;itime<nbTrig;itime++)  fEvtTime[itime]= 0;
        int date0, time0; //the date and time of the first event
        TDatime t0; // to store the date and time of the first event
        TDatime t1; // to store the date and time of the 5000th event in a run
        int evtLedTrigs = 0;
        int evtLedTrigs_PMT[8]; for(int mm=0;mm<8;mm++){evtLedTrigs_PMT[mm]=0;}// event numbers in LedTrigs, in which PMTs get >0 photons.
        for(int i=0;i<nentries;i++){ // begin to read events in one ROOT tree
            totalEvtNb++;
            //if(i%2000 ==0) cout<<i<<endl;
            //if(i==1000) break;
            bEvtNumber->GetEntry(i);
            bEvtTime->GetEntry(i);
            bTemperatureBox->GetEntry(i);
            bTemperatureRack->GetEntry(i);
            bLevel1->GetEntry(i);
            bLevel2->GetEntry(i);

            bResistivity->GetEntry(i);
            bQDC1->GetEntry(i);
            bQDC2->GetEntry(i);
            bTDC->GetEntry(i);
            bScaler->GetEntry(i);
            bDigitizerCharge->GetEntry(i);
            bDigitizerAmplitude->GetEntry(i);
            bDigitizerAmplitudeTimeBin->GetEntry(i);
            bDigitizerStartTimeBin->GetEntry(i);
            bDigitizerRiseTime->GetEntry(i);

            bDigitizerPulseWidth->GetEntry(i);
            bDigitizerPedMean->GetEntry(i);
            bDigitizerPedWidth->GetEntry(i);
            bTrigType->GetEntry(i);
            bEvtTime->GetEntry(i);

            bTrigFlag_Multi->GetEntry(i);
            bTrigFlag_Hodo->GetEntry(i);
            bTrigFlag_Led->GetEntry(i);
            bDigitizerNumberOfPulses->GetEntry(i);

            bDigitizerPulseCharge->GetEntry(i);
            bDigitizerPulseStartBin->GetEntry(i);
            bDigitizerPulseEndBin->GetEntry(i);
            bDigitizerPulseAmplitude->GetEntry(i);

            //define trigger types -- begins
            bool isledtrig = false;
            trigcnts_all++;
            int newtrigtype = 0;
			if(evtTDC[8]!=-1){ // LED triggers
              trigcnts_led++;
              evtLedTrigs++;
              isledtrig = true;
            }
            if(evtTDC[6]!=-1 && evtTDC[2]==-1 && evtTDC[0]!=-1 && evtTDC[3]!=-1 && evtTDC[1]!=-1){ // hodoscope triggers
              trigcnts_hodo++;
            }
            if( evtTDC[7]!=-1){ // multi triggers only!
              trigcnts_multi++;
            }
            if(evtTDC[7]!=-1 && evtTDC[6]!=-1 && evtTDC[2]==-1 && evtTDC[0]!=-1 && evtTDC[3]!=-1 && evtTDC[1]!=-1){  // multi and hodo triggers
              trigcnts_multi_and_hodo++;
            }
            if( (evtTDC[7]!=-1 || evtTDC[6]!=-1) && evtTDC[2]==-1 && evtTDC[0]!=-1 && evtTDC[3]!=-1 && evtTDC[1]!=-1 ){ // multi or hodo triggers, not led trigger
              trigcnts_multi_or_hodo++;
            }
	    if(evtTDC[6]!=-1 && evtTDC[4]==-1 && evtTDC[5]!=-1 && evtTDC[2]==-1 && evtTDC[0]!=-1 && evtTDC[3]!=-1 && evtTDC[1]!=-1){ // hodoTrig + H5 + !H4
                trigcnts_hodo_H5_notH4++;
	    }
	    if(evtTDC[6]!=-1 && evtTDC[4]==-1 && evtTDC[5]==-1 && evtTDC[2]==-1 && evtTDC[0]!=-1 && evtTDC[3]!=-1 && evtTDC[1]!=-1){ // hodoTrig + !H5 + !H4
	        trigcnts_hodo_notH5_notH4++;
	    }
	    if(evtTDC[6]!=-1 && evtTDC[4]!=-1 && evtTDC[5]==-1 && evtTDC[2]==-1 && evtTDC[0]!=-1 && evtTDC[3]!=-1 && evtTDC[1]!=-1){ // hodoTrig + !H5 + H4
		trigcnts_hodo_notH5_H4++;
	    }
	    if(evtTDC[6]!=-1 && evtTDC[4]!=-1 && evtTDC[5]!=-1 && evtTDC[2]==-1 && evtTDC[0]!=-1 && evtTDC[3]!=-1 && evtTDC[1]!=-1){ // hodoTrig + H5 + H4, this should be 0
		trigcnts_hodo_H5_H4++;
	    }

			
	    if(evtTDC[6]!=-1 && evtTDC[2]==-1 && evtTDC[0]!=-1 && evtTDC[3]!=-1 && evtTDC[1]!=-1){
		selectHodo++; 
		selectHodoTrig=1;
		if(evtTDC[7]!=-1) testcnt1++;
		if(evtTDC[8]!=-1) testcnt2++;
	    } // the TDC_hodo
	    if( (evtTDC[0]!=-1 || evtTDC[2]!=-1) && (evtTDC[1]!=-1 || evtTDC[3]!=-1) ){selectHs++; selectHsTrig=1;}
		//hTestHodoSelection->Fill(selectHodoTrig, selectHsTrig);
                bool testTrig = false;
            /*	
            if(evtTDC[6]!=-1){
                newtrigtype = 0;
                istrigger[newtrigtype]=true;
                triggercnt[newtrigtype]++;
			  
	        if(evtTDC[2]!=-1 && evtTDC[0]!=-1) hTDC_H0H2->Fill(evtTDC[0]-evtTDC[2]);
		if(evtTDC[2]!=-1 && evtTDC[1]!=-1) hTDC_H1H2->Fill(evtTDC[1]-evtTDC[2]);
		if(evtTDC[2]!=-1 && evtTDC[3]!=-1) hTDC_H3H2->Fill(evtTDC[3]-evtTDC[2]);
		if(evtTDC[0]!=-1 && evtTDC[1]!=-1) hTDC_H1H0->Fill(evtTDC[1]-evtTDC[0]);
		if(evtTDC[0]!=-1 && evtTDC[3]!=-1) hTDC_H3H0->Fill(evtTDC[3]-evtTDC[0]);
		if(evtTDC[2]!=-1 && evtTDC[5]!=-1) hTDC_H5H2->Fill(evtTDC[5]-evtTDC[2]);
		if(evtTDC[2]!=-1 && evtTDC[4]!=-1) hTDC_H4H2->Fill(evtTDC[4]-evtTDC[2]);
                for(int hi=0; hi<6;hi++){
		    if(evtTDC[hi]!=-1) hTimeDiff_hodoDet_hodoTrig[hi]->Fill((evtTDC[hi]-evtTDC[6])*0.035);
		}
            }
            */
            //selectHodoTrig=selectHsTrig=0;
            //continue;

            //now make my trigger conbinations based on the normal hodo trigger
            int PsTrigFlag[6]={0,0,0,0,0,0};
            for(int pscnt=0; pscnt<6; pscnt++){
                if(evtTDC[pscnt]!=-1){
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
            for(int trigindex=0; trigindex<nbTrig; trigindex++){
                int testtrignum = std::strtol(sTriggerType[trigindex].c_str(), 0, 2);//astring.Atoi();
                //cout<<strig[trigindex]<<"\t"<<testtrignum<<"\t"<<thistriggerflag<<endl;
                if(thistriggerflag == testtrignum){
                    istrigger[trigindex]=true;
                    triggercnt[trigindex]++;
                    //eventStartTime[trigindex] = myTriggerHitTime[indexMap_forEvtStartTime[trigindex]];
                    break;
                }
            }

            // define trigger types -- ends

            //get event rate using the information in Scaler, and take average from each processed run
            // t00 is the timestamp of the first event in the first run,
            // t1 is the timestamp of the event in the middle of each run
            if(i==0){ 
              firstEvtTime = double(evtT);
              //t0.Set(date0,time0);
              if(aa==0){ 
                TTimeStamp timestamp0(firstEvtTime);
                date0 = timestamp0.GetDate(0); // get date and time, not in UTC zone
                time0 = timestamp0.GetTime(0);
                t00.Set(date0,time0);
              }
            }
            else{
              if(i==(int)(nentries/2-1)){
                TTimeStamp timestamp1(evtT); // get date and time, for the 5000th event
                date0=timestamp1.GetDate(0);
                time0 = timestamp1.GetTime(0);
                t1.Set(date0,time0);
              }
            }
            // event rate from scaler data, fill histograms
            if(i!=0){
              for(int scalerid=0; scalerid<17; scalerid++){
                hEvtRate[scalerid]->Fill(evtScaler[scalerid]);
              }
            }
           // event rate vs. run time, at the last event, get average event rate for the three trigger types
            if(i==nentries-1){
              gEvtRateVsTime[2]->SetPoint(aa,t1.Convert(0)-t00.Convert(0),hEvtRate[8]->GetMean()); // multi triggers
              gEvtRateVsTime[1]->SetPoint(aa,t1.Convert(0)-t00.Convert(0),hEvtRate[7]->GetMean()); // hodo triggers
              gEvtRateVsTime[0]->SetPoint(aa,t1.Convert(0)-t00.Convert(0),hEvtRate[9]->GetMean()); // led triggers
              //because I need to get mean values for each run, after each run I need to clear the histogram
              for(int ch=0; ch<17; ch++){
                  hEvtRate[ch]->Reset("ICESM");
              }
            }

            if(i%100==0){//fill level vs. temperature, level vs. rate, and hodoRate vs multiRate every 100 events
              hLevelVsTemperature->Fill(evtTempBox-15, evtLevel1);
              hLevelVsMultiEvtRate->Fill(evtScaler[8], evtLevel1); // multiTrigRate vs level
              //hHodoRateVsMultiRate->Fill(evtScaler[8], evtScaler[7]);
            }
            // for temperature and resistivity
            if(i%2000==0){
              gTempInBox->SetPoint(totalEnvEvtNb, evtT, evtTempBox-15);
              gTempOnRack->SetPoint(totalEnvEvtNb, evtT, evtTempRack);
              gResistivity->SetPoint(totalEnvEvtNb, evtT, evtResist);
              //cout<<"Time "<<evtT<<endl;
              totalEnvEvtNb++;
            }

            //save to TDC histograms
            for(int mm=0;mm<16;mm++){
                //hTDC[mm]->Fill(evtTDC[mm]); // the raw TDC entries
                //hTDC1[mm]->Fill(0.1*(evtTDC[mm]-evtTDC[9]));
                for(int trigindex=0; trigindex<nbTrig; trigindex++){
                    if(istrigger[trigindex]==true){
                        hTDCnew[mm][trigindex]->Fill(evtTDC[mm]);
			hTDC[mm][trigindex]->Fill(evtTDC[mm]*0.035);
                    }
                }
            }
            // test #1: I notice for S2 sometimes pedestal mean is higher, want to find those events and check their waveforms
            // if(evtDigiPM[2]>0.7){cout<<"evt with large pedestal mean: evt id "<<i<<endl;}

            double totalpmtcharge = 0; // total pmt charge in npe. in the full 2500 ns window
            double totalpmtcharge_top = 0; // total pmt charge in top PMTs in npe in the full time window
            double totalpmtcharge_bottom = 0; // total pmt charge in bottom PMTs in npe in the full time window
            double totalpmtcharge1 = 0; // total pmt charge in npe. from the first pulse
            double totalpmtcharge_top1 = 0; // total pmt charge in top PMTs in npe from the first pulse
            double totalpmtcharge_bottom1 = 0; // total pmt charge in bottom PMTs from the first pulse
            double totalpmtchargeS2S3 = 0;
            int totalNbOfPulses = 0;
            double thispmt_charge[pmtnb], thispmt_charge_50ns[pmtnb];
            for(int t=0; t<pmtnb; t++){ thispmt_charge[t]=0; thispmt_charge_50ns[t]=0;}
            // my vectors to hold pulses in the PMTs:
            vector< vector<double> > pulseTimeInPMTs; // pulse time
            pulseTimeInPMTs.resize(pmtnb,vector<double>(0,0)); // 8 columns for the 8 pmts, each contains 0 pulses in the initialization
            vector< vector<double> > pulseChargeInPMTs; // pulse charge
            pulseChargeInPMTs.resize(pmtnb,vector<double>(0,0)); // 8 columns for the 8 pmts, each contains 0 pulses in the initialization
            for(int k=0;k<pmtnb;k++){ // begin looping through pmts
                /*
                //save to QDC1 histograms
                if(strcmp(sQDC1ChMap[k].c_str(),"")!=0){
                    //cout<<k<<"\t"<<evtQDC1[k]<<endl;
                   for(int trigindex=0; trigindex<nbTrig; trigindex++){
                     if(istrigger[trigindex]==true){
                       hQDC1[k][trigindex]->Fill(evtQDC1[k]);
                     }
                   }
                }
                //save to QDC2 histograms
                if(strcmp(sQDC2ChMap[k].c_str(),"")!=0){
                   for(int trigindex=0; trigindex<nbTrig; trigindex++){
                     if(istrigger[trigindex]==true){
                       hQDC2[k][trigindex]->Fill(evtQDC2[k]);
                     }
                   }
                }
                */
              // now look at digitizers
              // save pedestal mean, width and threshold, no need to seperate into different trigger types
              hPedestalMean[k]->Fill(evtDigiPM[k]);
              hPedestalWidth[k]->Fill(evtDigiPS[k]);
              //hWaveformThreshold[k]->Fill(5.0*evtDigiPS[k]-evtDigiPM[k]);
              // pmt pulses charge, time etc.
              double evtTotalCharge = 0.0; // total charge in the full time window, note this is not converted to npe yet
              double evtTotalCharge1 = 0.0; // total charge of the first pulse,
              // I want to remove those pulses with <=20 amplitude
              // also for each pulse I require its charge to be >=0.5 pe.
              vector<double> newEvtDigitizerPulseCharge;
              vector<double> newEvtDigitizerPulseStartBin;
              vector<double> newEvtDigitizerPulseEndBin;
              vector<double> newEvtDigitizerPulseAmplitude;
              double ledtrig_pmtcharge = 0;
              //cout<<"number of pulses before "<<evtDigitizerNumberOfPulses[k]<<endl;

             for(int trig=0; trig<nbTrig; trig++){
                 if(istrigger[trig]==true){
                     for(int pulsecnt=0 ; pulsecnt<evtDigitizerNumberOfPulses[k]; pulsecnt++){
                         hDigiPulseAmp[k][trig]->Fill(-1.0*evtDigitizerPulseAmplitude[k][pulsecnt]);
                     }
                 }
             }

              for(int pulsecnt=0; pulsecnt<evtDigitizerNumberOfPulses[k]; pulsecnt++){
                  //cout<< k << "\t" << pulsecnt << "\t" << evtDigitizerPulseCharge[k][pulsecnt] << "\t" << evtDigitizerPulseStartBin[k][pulsecnt] <<endl;
                  if(isledtrig==true){ledtrig_pmtcharge += -1.0*evtDigitizerPulseCharge[k][pulsecnt]/converterfactor[k];}
                  if(evtDigitizerPulseAmplitude[k][pulsecnt]>=-15.0  // the found pulse has amplitude <=15 ADC counts is not good because spe amplitude ~50 ADC counts
                       || (-1.0*evtDigitizerPulseCharge[k][pulsecnt]/converterfactor[k])<0.5 // the charge for the pulse is not good if < 0.5 p.e.,
                  ) continue;
                  // otherwise I keep this as a good pulse
                  newEvtDigitizerPulseCharge.push_back(evtDigitizerPulseCharge[k][pulsecnt]);
                  newEvtDigitizerPulseStartBin.push_back(evtDigitizerPulseStartBin[k][pulsecnt]);
                  newEvtDigitizerPulseEndBin.push_back(evtDigitizerPulseEndBin[k][pulsecnt]);
                  newEvtDigitizerPulseAmplitude.push_back(evtDigitizerPulseAmplitude[k][pulsecnt]);
                  pulseTimeInPMTs[k].push_back(evtDigitizerPulseStartBin[k][pulsecnt]);
                  pulseChargeInPMTs[k].push_back(-1.0*evtDigitizerPulseCharge[k][pulsecnt]/converterfactor[k]); // pulse charge in npe
             }
             evtDigitizerNumberOfPulses[k] = newEvtDigitizerPulseCharge.size();
             //cout<<"number of pulses after "<<evtDigitizerNumberOfPulses[k]<<endl;
             double t_firstpulse = 0;
             for(int pulsecnt=0; pulsecnt<evtDigitizerNumberOfPulses[k]; pulsecnt++){
                 //cout<< k <<"\t" << pulsecnt <<"\t" << pulseChargeInPMTs[k][pulsecnt] << "\t" << pulseTimeInPMTs[k][pulsecnt] <<endl;
                 thispmt_charge[k] += pulseChargeInPMTs[k][pulsecnt];
                 //cout<<k<<"\t"<<pulsecnt<<"\t"<<pulseTimeInPMTs[k][pulsecnt]<<"\t"<<meanPulseStartBin[k]<<endl;
                 if(pulsecnt==0){
                     t_firstpulse = pulseTimeInPMTs[k][pulsecnt];
                     //cout << "test output " << k <<"\t"<< t_firstpulse << endl;
                     for(int trig=0; trig<nbTrig; trig++) {
                         if(istrigger[trig]==true)
                             //cout<<k<<"\t"<<trig<<"\t"<<sTriggerType[trig]<<"\t"<<t_firstpulse<<endl;
                             hDigiST[k][trig]->Fill(t_firstpulse);
                     }
                     
                 }
                 if(evtDigitizerNumberOfPulses[k]==1){ thispmt_charge_50ns[k]=pulseChargeInPMTs[k][pulsecnt];}
                 else{
                     if( (pulseTimeInPMTs[k][pulsecnt]-t_firstpulse) <=50 ){
                         thispmt_charge_50ns[k]+=pulseChargeInPMTs[k][pulsecnt];
                     }
                 }
             }
             
             totalNbOfPulses += evtDigitizerNumberOfPulses[k];
             if(isledtrig==true) evtLedTrigs_total++;
             if(isledtrig==true && ledtrig_pmtcharge>0){  // record the number of led trigger in PMT with >0 npe charge
                 evtLedTrigs_PMT[k]++;
                 if(k==4 || k==5) evtLedTrigs_topPMT++;
                 else evtLedTrigs_botPMT++;
             }
             // now begin my process
             if(evtDigitizerNumberOfPulses[k]>0) nhit++; 
             /*
             for(int pulsecnt=0; pulsecnt<evtDigitizerNumberOfPulses[k]; pulsecnt++){
                 evtTotalCharge += newEvtDigitizerPulseCharge[pulsecnt];
                 if(pulsecnt==0) evtTotalCharge1 += newEvtDigitizerPulseCharge[pulsecnt];
             }
             totalpmtcharge += -1.0*evtTotalCharge/converterfactor[k]; // sum up the pmt charge in the full window
             totalpmtcharge1 += -1.0*evtTotalCharge1/converterfactor[k]; // sum up the pmt charge of the first pulse
             if(k==4 || k==5){
               totalpmtcharge_top += -1.0*evtTotalCharge/converterfactor[k]; // sum up the top pmt charge
               totalpmtcharge_top1 += -1.0*evtTotalCharge1/converterfactor[k]; // sum up the top pmt charge
             }
             if(k!=4 && k!=5){
               totalpmtcharge_bottom += -1.0*evtTotalCharge/converterfactor[k]; // sum up the top pmt charge
               totalpmtcharge_bottom1 += -1.0*evtTotalCharge1/converterfactor[k]; // sum up the top pmt charge
             }
             if(k==2 || k==3){
               totalpmtchargeS2S3 += -1.0*evtTotalCharge1/converterfactor[k]; // sum up charge in pmt S2 and S3
             }
             */
             //for(int trigindex=0; trigindex<nbTrig; trigindex++){
               //if(istrigger[trigindex]==true){
                 //hNumberOfPulses[k][trigindex]->Fill(evtDigitizerNumberOfPulses[k]);
                 //hDigiAmp[k][trigindex]->Fill(evtDigiAmplitude[k]); //pulse's amplitude, global
                 //hDigiAmpTimeBin[k][trigindex]->Fill(evtDigiAmplitudeTimeBin[k]+80); //pulse's amplitude time bin
                 //hScat_DigiAmp_vs_TimeBin[k][trigindex]->Fill(evtDigiAmplitudeTimeBin[k]+80,evtDigiAmplitude[k]);
                 //if(evtDigitizerNumberOfPulses[k]>0){//pulse start time vs. TDC trigger time, remeber no TDC entries for S6 and S7
                   //hScat_TDC_vs_PMTT[k][trigindex]->Fill(newEvtDigitizerPulseStartBin[0],evtTDC[6]); //evtTDC[TDCIndexMap[k]]
                   //if(k<6) hScat_TDC_vs_PMTT1[k][trigindex]->Fill(newEvtDigitizerPulseStartBin[0],0.035*evtTDC[TDCIndexMap[k]]); //
                   //hTDC_minus_PMTT[k][trigindex]->Fill(newEvtDigitizerPulseStartBin[0]-0.035*evtTDC[6]);
                 //}
                 //for(int mm=0; mm<evtDigitizerNumberOfPulses[k];mm++){
                   //if(mm==0){
                     //hDigitizerPulseCharge1[k][trigindex]->Fill(-1.0*newEvtDigitizerPulseCharge[mm]/converterfactor[k]);
                     //hDigiPulseChargeNpe_vs_PulseWidth[k][trigindex]->Fill(-1.0*newEvtDigitizerPulseCharge[mm]/converterfactor[k],newEvtDigitizerPulseEndBin[mm]-newEvtDigitizerPulseStartBin[mm]);
					 //hTimeDiff_PMTs_hodoTrig[k]->Fill(newEvtDigitizerPulseStartBin[mm]-evtTDC[6]*0.035);
                   //}
                   //hDigitizerPulseCharge_vs_PulseTime[k][trigindex]->Fill(-1.0*newEvtDigitizerPulseCharge[mm]/converterfactor[k],newEvtDigitizerPulseStartBin[mm]);
                   //hDigitizerPulseWidth_vs_PulseTime[k][trigindex]->Fill(newEvtDigitizerPulseEndBin[mm]-newEvtDigitizerPulseStartBin[mm],newEvtDigitizerPulseStartBin[mm]);
                   //hDigiPulseAmp[k][trigindex]->Fill(-1.0*newEvtDigitizerPulseAmplitude[k]);
                   //hDigiST_all[k][trigindex]->Fill(newEvtDigitizerPulseStartBin[mm]);
                 //}
                 /*
                 if(evtDigitizerNumberOfPulses[k]>1){
                     double dT2 = newEvtDigitizerPulseStartBin[1]-newEvtDigitizerPulseStartBin[0];
                     hDigitizerPulse_TimeDifference[k][trigindex]->Fill(dT2);
                     //check correlation between after pulse charge and initial pulse charge
                     // this has to be done for events with >=2 pulses. S4 has fewest after pulse, here I don't use S4. The numbers here I read from the hDigitizerPulse_TimeDifference histograms
                     if( k==0 && ( (dT2>60 && dT2<140) || (dT2>340 && dT2<440) ) ){
                       hAfterPulseChargeCorrelation[k][trigindex]->Fill(-1.0*newEvtDigitizerPulseCharge[0]/converterfactor[k], -1.0*newEvtDigitizerPulseCharge[1]/converterfactor[k]);
                     }
                     if( (k==1 || k==3) && ( (dT2>100 && dT2<160) || (dT2>340 && dT2<440) ) ){
                       hAfterPulseChargeCorrelation[k][trigindex]->Fill(-1.0*newEvtDigitizerPulseCharge[0]/converterfactor[k], -1.0*newEvtDigitizerPulseCharge[1]/converterfactor[k]);
                     }
                     if( (k==5) && ( (dT2>360 && dT2<460) ) ){
                       hAfterPulseChargeCorrelation[k][trigindex]->Fill(-1.0*newEvtDigitizerPulseCharge[0]/converterfactor[k], -1.0*newEvtDigitizerPulseCharge[1]/converterfactor[k]);
                     }
                     if( (k==2) && ( (dT2>100 && dT2<160)  ) ){
                       hAfterPulseChargeCorrelation[k][trigindex]->Fill(-1.0*newEvtDigitizerPulseCharge[0]/converterfactor[k], -1.0*newEvtDigitizerPulseCharge[1]/converterfactor[k]);
                     }
                     if( (k==4) && ( (dT2>380 && dT2<460)  ) ){
                       hAfterPulseChargeCorrelation[k][trigindex]->Fill(-1.0*newEvtDigitizerPulseCharge[0]/converterfactor[k], -1.0*newEvtDigitizerPulseCharge[1]/converterfactor[k]);
                     }
                     if( (k==6 || k==7) && ( (dT2>140 && dT2<160) || (dT2>220 && dT2<320) ) ){
                       hAfterPulseChargeCorrelation[k][trigindex]->Fill(-1.0*newEvtDigitizerPulseCharge[0]/converterfactor[k], -1.0*newEvtDigitizerPulseCharge[1]/converterfactor[k]);
                     }
                 }
                 */
                 //hDigitizerPulseCharge_npe[k][trigindex]->Fill(-1.0*evtTotalCharge/converterfactor[k]);
                 //hDigiChargeNpe_vs_DigiAmp[k][trigindex]->Fill(evtDigiAmplitude[k], -1.0*evtTotalCharge/converterfactor[k]);

                 //if(evtDigitizerNumberOfPulses[k]>0){
                 //    hDigiST[k][trigindex]->Fill(newEvtDigitizerPulseStartBin[0]);
                     //hDigiPW[k][trigindex]->Fill(newEvtDigitizerPulseEndBin[0]-newEvtDigitizerPulseStartBin[0]);
                 //}
               //}
             //}//end of filling histograms for one pmt
        } // end of looping through pmts


        for(int trigindex=0; trigindex<nbTrig; trigindex++){
            for(int pmtID=0; pmtID<pmtnb; pmtID++){
                if(istrigger[trigindex]==true){
                    hCharge[pmtID][trigindex]->Fill(thispmt_charge[pmtID]);
                    hCharge1[pmtID][trigindex]->Fill(thispmt_charge_50ns[pmtID]);
                }
            }
	}

        // now for summed pmt charge
        for(int index=0; index<nbTrig; index++){
          if(istrigger[index]==true){
            //hDigitizerPulseChargeSum[index]->Fill(totalpmtcharge); // all triggers
            //hDigitizerPulseChargeSumTop[index]->Fill(totalpmtcharge_top); // all triggers
            //hDigitizerPulseChargeSumBottom[index]->Fill(totalpmtcharge_bottom); // all triggers
            //hDigitizerPulseChargeSumCorrelation[index]->Fill(totalpmtcharge_bottom, totalpmtcharge_top); // all triggers
            //hDigitizerPulseChargeSum1[index]->Fill(totalpmtcharge1); // all triggers
            //hDigitizerPulseChargeSumTop1[index]->Fill(totalpmtcharge_top1); // all triggers
            //hDigitizerPulseChargeSumBottom1[index]->Fill(totalpmtcharge_bottom1); // all triggers
            int numberOfCoincidence = 0; //
            int numberOfCoincidenceNew = 0;
            double nbsigma = 3.0;
            // get the time difference between pmts, for all photons
            /*
            for(int ngroup = 0; ngroup<7; ngroup++){
              for(int nbpulseS0=0; nbpulseS0<pulseTimeInPMTs[ngroup].size(); nbpulseS0++){
                for(int pmtid=ngroup+1;pmtid<8;pmtid++){
                  for(int nbpulseS1=0; nbpulseS1<pulseTimeInPMTs[pmtid].size(); nbpulseS1++){
                    double dT = (pulseTimeInPMTs[pmtid][nbpulseS1]-meanPulseStartBin[pmtid]) - (pulseTimeInPMTs[ngroup][nbpulseS0]-meanPulseStartBin[ngroup]);
                    if(ngroup==0) hPulseTimeDiffToS0[pmtid-1-ngroup][index]->Fill( dT );
                    if(ngroup==1) hPulseTimeDiffToS1[pmtid-1-ngroup][index]->Fill( dT );
                    if(ngroup==2) hPulseTimeDiffToS2[pmtid-1-ngroup][index]->Fill( dT );
                    if(ngroup==3) hPulseTimeDiffToS3[pmtid-1-ngroup][index]->Fill( dT );
                    if(ngroup==4) hPulseTimeDiffToS4[pmtid-1-ngroup][index]->Fill( dT );
                    if(ngroup==5) hPulseTimeDiffToS5[pmtid-1-ngroup][index]->Fill( dT );
                    if(ngroup==6) hPulseTimeDiffToS6[pmtid-1-ngroup][index]->Fill( dT );
                    if( dT > coincidenceTimeMeanMatrix[pmtid-1-ngroup][ngroup]-nbsigma*coincidenceTimeWidthMatrix[pmtid-1-ngroup][ngroup] && dT < coincidenceTimeMeanMatrix[pmtid-1-ngroup][ngroup]+nbsigma*coincidenceTimeWidthMatrix[pmtid-1-ngroup][ngroup]){
                      numberOfCoincidence++;
                    }
                  }
                }
              }
            }
            */
            /*
            vector< pair <double, int> > mytimevector; // I use this vector to hold pair of pulse time and pmtid
            vector< pair <double, double> > mychargevector; // I use this vector to hold pair of pulse time and pulse charge. I will sort by time
            for(int pmtid=0;pmtid<8;pmtid++){ // no selection, fill all pulse times (subtract muon start time in each pmt) in all pulses in all pmts
              for(int nbpulse=0; nbpulse<pulseTimeInPMTs[pmtid].size(); nbpulse++){
                // if the first pulse is not in the muon start region, then I use the average muon start time as the start time
                // otherwise, that means the first pulse is the true muon start time, so I use it directly as the start time
                double deltaT=0;// = pulseTimeInPMTs[pmtid][nbpulse]-evtTDC[6]*0.035-meanPulseStartBin[pmtid];
                if(nbpulse==0 && pulseTimeInPMTs[pmtid][0]>=upperBoundFirstPulseStartBin[pmtid])
                  deltaT = pulseTimeInPMTs[pmtid][nbpulse]-evtTDC[6]*0.035-meanPulseStartBin[pmtid];
                else deltaT = pulseTimeInPMTs[pmtid][nbpulse]-pulseTimeInPMTs[pmtid][0];
                hTimeDiff[index]->Fill(deltaT);
                mytimevector.push_back( make_pair(deltaT, pmtid) );
                mychargevector.push_back( make_pair(deltaT, pulseChargeInPMTs[pmtid][nbpulse]) );
              }
            }
            
            // sort my vector by pulse time
            sort(mytimevector.begin(), mytimevector.end());
            sort(mychargevector.begin(), mychargevector.end());
            */
            // search for time coincidences in this vector. 
            // if the length of the vector is <=1, there must be no coincidenced time
            /*
            int cntindex=0; 
            if(mytimevector.size()>1){
              double thistime = mytimevector[0].first;
              double thischarge = mychargevector[0].second;
              int thispmt = mytimevector[0].second;
              int indexR, indexC; 
              int mycoincidencecounting[20]; for(int yy=0;yy<20;yy++) mycoincidencecounting[yy]=0;
              vector< vector<int> > thecoincidencepmt; thecoincidencepmt.resize(1); //thecoincidencepmt[0].resize(1);
              vector< vector<double> > thecoincidencetime; thecoincidencetime.resize(1); //thecoincidencetime[0].resize(1);
              vector< vector<double> > thecoincidencecharge; thecoincidencecharge.resize(1);
              int hidenflag = 0;
              //mycoincidencetiming.push_back(thistime);

              for(int vi=1; vi<mytimevector.size(); vi++){
                double atime = mytimevector[vi].first;
                double acharge = mychargevector[vi].second;
                int apmt = mytimevector[vi].second;
                if(apmt==thispmt){thistime = atime; thispmt = apmt; continue; }
                if(thispmt>apmt){ indexR = thispmt; indexC = apmt; }
                else{             indexC = thispmt; indexR = apmt; }
                double rangeleft  = coincidenceTimeMeanMatrix[indexR][indexC]-5.0*coincidenceTimeWidthMatrix[indexR][indexC];
                double rangeright = coincidenceTimeMeanMatrix[indexR][indexC]+5.0*coincidenceTimeWidthMatrix[indexR][indexC];
                if( (atime-thistime) >= rangeleft && (atime-thistime) <= rangeright ){ // this time is in conincidence with the previous one
                  if(hidenflag==0){ // I use a "hiden" flag to record the 0th pulse in mytimevector, which I put outside this for loop
                    thecoincidencepmt[cntindex].push_back( thispmt );
                    thecoincidencetime[cntindex].push_back( thistime );
                    thecoincidencecharge[cntindex].push_back( thischarge );
                    mycoincidencecounting[cntindex]+=1;
                    hidenflag=1;  // I only need to save that pair once, so here I must change the flag value to by pass this if-structure
                  }
                  mycoincidencecounting[cntindex]+=1;
                  thecoincidencepmt[cntindex].push_back( apmt );
                  thecoincidencetime[cntindex].push_back( atime );
                  thecoincidencecharge[cntindex].push_back( acharge );
                }
                else{ // this time must belong to the next group
                  cntindex++;
                  thecoincidencepmt.resize(cntindex+1); //thecoincidencepmt[cntindex].resize(1);
                  thecoincidencetime.resize(cntindex+1); //thecoincidencetime[cntindex].resize(1);
                  thecoincidencecharge.resize(cntindex+1);
                  hidenflag=0;
                  //mycoincidencetiming.push_back(atime);
                }
                thistime = atime; thispmt = apmt; thischarge = acharge;
              }

              cntindex=0;
              for(int yy=0;yy<20;yy++){
                if(mycoincidencecounting[yy]>1){
                  cntindex++;
                }
              }
              if(cntindex>=2){test_cnt1++; if(cntindex>=3) test_cnt2++;}
              if(cntindex==2){ // there are very few events: in 50 water runs, >=2 events: 39 , >=3 events: 1, the ratio is 1/39=2.6%.
                //cout<<aa<<"\t"<<i<<"\t\n";
                double electroncharge = 0;
                double avgtime = 0;
                int decaypulsecnt = 0; // to count how many pmts have hits by the decayed electrons
                int decaypulsecnt_top = 0;  // to count how many top PMTs have hits by the decayed electrons, should be <=2
                int decaypulsecnt_bot = 0;  // to count how many bottom PMTs have hits by the decayed electrons, should be <=6
                double electroncharge_top = 0;
                double electroncharge_bot = 0;
                double avgtime_top=0, avgtime_bot=0;
                int decayflag = 0;
                int coincidencecnt[2];
                for(int yy=0;yy<20;yy++){
                  if(mycoincidencecounting[yy]>1){
                    decayflag++;
                    //for(int pp=0; pp<thecoincidencetime[yy].size(); pp++){
                    //  cout<<aa<<"\t"<<i<<"\t"<<yy<<"\t"<<decayflag<<"\t"<<thecoincidencepmt[yy][pp]<<"\t"<<thecoincidencetime[yy][pp]<<"\t"<<thecoincidencecharge[yy][pp]<<endl;
                    //}
                    if(decayflag==1) coincidencecnt[decayflag-1] = thecoincidencetime[yy].size();
                    if(decayflag==2){
                      coincidencecnt[decayflag-1] = thecoincidencetime[yy].size();
                      for(int pp=0; pp<thecoincidencetime[yy].size(); pp++){
                        avgtime += thecoincidencetime[yy][pp];
                        electroncharge += thecoincidencecharge[yy][pp];
                        if(thecoincidencepmt[yy][pp]==4 || thecoincidencepmt[yy][pp]==5){
                          decaypulsecnt_top++;
                          avgtime_top += thecoincidencetime[yy][pp];
                          electroncharge_top += thecoincidencecharge[yy][pp];
                        }
                        else{
                          decaypulsecnt_bot++;
                          avgtime_bot += thecoincidencetime[yy][pp];
                          electroncharge_bot += thecoincidencecharge[yy][pp];
                        }
                      }
                      avgtime /= thecoincidencetime[yy].size();
                    }
                  }
                }
                hCoincidence[index]->Fill(numberOfCoincidence);
                hTimeDiffSelect[index]->Fill(avgtime);
                hDigitizerPulseChargeDecay[index]->Fill(electroncharge); // now this must be the total charge of the decay electron caused pulses
                hDecayChargeVsTime[index]->Fill(avgtime,electroncharge);
                hCoincidenceVsDecayTime[index]->Fill(avgtime,numberOfCoincidence);
                hCoincidenceNumbers[index]->Fill(coincidencecnt[0],coincidencecnt[1]);
                if(decaypulsecnt_top!=0 && decaypulsecnt_bot!=0){
                  hChargeRatioVsTimeInDecay1[index]->Fill(avgtime_top/decaypulsecnt_top-avgtime_bot/decaypulsecnt_bot,electroncharge_top/(electroncharge));
                }
              }
            }
            */
            /*
            bool store_decay_pulses = false;
            if(cntindex==2 && store_decay_pulses ==true){
              rtfile->cd(); // the root file to save information
              sprintf(tempname,"hTimeDiff_FileID_%d_Evt_%d",aa,i);
              TH1F* hTimeDiffOneEvt = new TH1F(tempname,tempname,266,-100,2560); // 10 ns bin 
              for(int pmtid=0;pmtid<8;pmtid++){ // no selection, fill all pulse times (subtract muon start time in each pmt) in all pulses in all pmts
                for(int nbpulse=0; nbpulse<pulseTimeInPMTs[pmtid].size(); nbpulse++){
                  hTimeDiffOneEvt->Fill(pulseTimeInPMTs[pmtid][nbpulse]-evtTDC[6]*0.035-meanPulseStartBin[pmtid]);
                }
              }
              hTimeDiffOneEvt->Write();
              dirRebuildPulse->cd();
              char pulsename[100]; sprintf(pulsename, "File_%d_Evt_%d",aa,i);
              TCanvas* cpulserebuilt = new TCanvas(pulsename,pulsename,1600,800);
              cpulserebuilt->Divide(1,8);
              for(int pmtid=0;pmtid<8;pmtid++){
                cpulserebuilt->cd(pmtid+1);
                sprintf(pulsename, "File_%d_Evt_%d_PMT_%d",aa,i,pmtid);
                TH1F* pulsehist = new TH1F(pulsename,pulsename,2560,0,2560);
                for(int nbpulse=0; nbpulse<pulseTimeInPMTs[pmtid].size(); nbpulse++){
                  pulsehist->SetBinContent((int)(pulseTimeInPMTs[pmtid][nbpulse]-evtTDC[6]*0.035-meanPulseStartBin[pmtid]+1), pulseChargeInPMTs[pmtid][nbpulse]);
                }
                pulsehist->Draw();
                pulsehist->GetXaxis()->SetLabelSize(0.1);
                pulsehist->GetYaxis()->SetLabelSize(0.1);
              }
              cpulserebuilt->Write();
              delete hTimeDiffOneEvt;
            }
            */
          }
          istrigger[index]=false;
        }

        //end of one event
        nhit=0; // reset nhit

     } // end of reading all events in one ROOT tree

     // fill the led event number ratios
     cout<<"Number of Led triggers >0 pe in PMTs: ";
     for(int mm=0;mm<pmtnb;mm++){cout<<evtLedTrigs_PMT[mm]<<", ";}
     cout<<"\t/\t total led triggers "<<evtLedTrigs<<endl;
//     cout<<endl;
     //for(int mm=0;mm<2;mm++){
     //  for(int nn=0;nn<6;nn++){
     //    if(evtLedTrigs_PMT[temp2index[nn]]==0){cout<<"no led trigger in PMT S"<<temp2index[nn]<<" in run "<<aa<<endl; continue;}
     //    gLedEffRatioVsTime[mm][nn]->SetPoint(aa,t1.Convert(0)-t00.Convert(0), (1.0*evtLedTrigs_PMT[temp1index[mm]])/(1.0*evtLedTrigs_PMT[temp2index[nn]]));
     //  }
     //}
     if( aa==NbOfFilesToRead-1 ){
//       if( (t1.Convert(0)-t00.Convert(0))%86400 -  )
       //calculate deviation
       cout<<"Time "<< t1.Convert(0)-t00.Convert(0) <<"\t"<< (t1.Convert(0)-t00.Convert(0))%86400 <<endl;
       double p1 = evtLedTrigs_topPMT/evtLedTrigs_total;
       double dv1 = TMath::Sqrt(evtLedTrigs_topPMT*(1-p1));
       double p2 = evtLedTrigs_botPMT/evtLedTrigs_total;
       double dv2 = TMath::Sqrt(evtLedTrigs_botPMT*(1-p2));
       double value = 3.0*evtLedTrigs_topPMT/evtLedTrigs_botPMT;
       double error = value*TMath::Sqrt((dv1/evtLedTrigs_topPMT)*(dv1/evtLedTrigs_topPMT)+(dv2/evtLedTrigs_botPMT)*(dv2/evtLedTrigs_botPMT));
       gLedEffRatioVsTimeTot->SetPoint(Nmonitor, t1.Convert(0)-t00.Convert(0), value);
       gLedEffRatioVsTimeTot->SetPointError(Nmonitor, 0, error);
       //fstream fout_ratio("LedTrig_TopBot_EvtNbRatio.txt",ios::out|ios::app);
       //fout_ratio<<t1.GetYear()<<"\t"<<t1.GetMonth()<<"\t"<<t1.GetDay()
       //    <<"\t"<<t1.GetHour()<<"\t"<<t1.GetMinute()<<"\t"<<t1.GetSecond()
       //    <<"\t"<<value<<"\t"<<error
	//	   <<"\t"<<hDigitizerPulseChargeSumBottom[0]->GetMean()
	//	   <<"\t"<<hDigitizerPulseChargeSumBottom[0]->GetRMS()<<endl;
       //fout_ratio.close();
       Nmonitor++;
       evtLedTrigs_topPMT=0; evtLedTrigs_botPMT=0; evtLedTrigs_total=0;
     }

     //fill the meanPMTnpe for stability monitornig purposes
     for(int pmtid=0; pmtid<pmtnb; pmtid++){
       //gMeanNPE[pmtid]->SetPoint(aa,t1.Convert(0)-t00.Convert(0), hDigitizerPulseCharge1[pmtid][2]->GetMean());
       gMeanNPE[pmtid]->SetPoint(aa,t1.Convert(0)-t00.Convert(0), hCharge1[pmtid][0]->GetMean());
     }

     rootfile->Close(); //delete rootfile;
     //cout<<" -- done"<<endl;
  }//end reading ROOT tree

  cout<< "finished reading all input root files..." <<endl;
// for temperature and resistivity, write to ROOT file
    /*
    dirEnv->cd();
    TCanvas* cEnv1 = new TCanvas();
    gTempInBox->Draw("AL"); gTempInBox->SetLineColor(kBlack);
    gTempOnRack->Draw("L"); gTempOnRack->SetLineColor(kBlue);
    gTempInBox->GetXaxis()->SetTimeDisplay(1);
    gTempInBox->GetXaxis()->SetTimeFormat("#splitline{%H:%M}{%m/%d}%F1970-01-01 00:00:00");
    gTempInBox->GetYaxis()->SetTitle("Temperature (#circC)");
    TLegend* leg = new TLegend(0.2,0.7,0.5,0.9);
    leg->AddEntry(gTempInBox,"Temperature in dark room","l");
    leg->AddEntry(gTempOnRack,"Temperature at electronics rack","l");
    leg->Draw();
    cEnv1->SetName("cTemperatures"); cEnv1->Write(); cEnv1->Close();
    TCanvas* cEnv2 = new TCanvas();
    gResistivity->Draw("AL");
    gResistivity->GetXaxis()->SetTimeDisplay(1);
    gResistivity->GetXaxis()->SetTimeFormat("#splitline{%H:%M}{%m/%d}%F1970-01-01 00:00:00");
    gResistivity->GetYaxis()->SetTitle("Resistivity (M#Omega#timescm)");
    cEnv2->SetName("cResistivity"); cEnv2->Write(); cEnv2->Close();
    TCanvas* cEnv3 = new TCanvas();
    TLegend* leg3 = new TLegend(0.75,0.5,0.9,0.9);
    for(int pmtid=0;pmtid<8;pmtid++){
      if(pmtid==0) gMeanNPE[pmtid]->Draw("AP");
      else gMeanNPE[pmtid]->Draw("P");
      gMeanNPE[pmtid]->SetMarkerStyle(8);
      gMeanNPE[pmtid]->SetMarkerColor(pmtid+1);
      leg3->AddEntry(gMeanNPE[pmtid],sDigiChMap[pmtid].c_str(),"lp");
    }
    leg3->Draw();
    gMeanNPE[0]->GetXaxis()->SetTimeDisplay(1);
    gMeanNPE[0]->GetXaxis()->SetTimeOffset(t00.Convert(0));
    gMeanNPE[0]->GetXaxis()->SetTimeFormat("#splitline{%H:%M}{%m/%d}");
//    gMeanNPE[0]->GetXaxis()->SetTimeFormat("#splitline{%H:%M}{%m-%d-%y}%F1970-01-01 00:00:00");
    gMeanNPE[0]->GetYaxis()->SetTitle("Mean npe");
    cEnv3->SetName("cMeanNPE"); cEnv3->Write(); cEnv3->Close();
    */
/*
    TCanvas* cLedEvtRatio[2];
    TLegend* legLedEvtRatio[2];
    for(int mm=0; mm<2; mm++){
      char name[100]; sprintf(name,"cLedEventRate_S%d",temp1index[mm]);
      cLedEvtRatio[mm] = new TCanvas(name,name,1200,600);
      legLedEvtRatio[mm] = new TLegend(0.7,0.6,0.9,0.9);
      for(int nn=0;nn<6;nn++){
        if(nn==0){
          gLedEffRatioVsTime[mm][nn]->Draw("AL");
          gLedEffRatioVsTime[mm][nn]->GetXaxis()->SetTimeDisplay(1);
          gLedEffRatioVsTime[mm][nn]->GetXaxis()->SetTimeOffset(t00.Convert(0));
          gLedEffRatioVsTime[mm][nn]->GetXaxis()->SetTimeFormat("%m/%d");
          gLedEffRatioVsTime[mm][nn]->GetYaxis()->SetTitle("Event ratio");
        }
        else{
          gLedEffRatioVsTime[mm][nn]->Draw("L");
        }
        //gLedEffRatioVsTime[mm][nn]->SetMarkerStyle(8);
        //gLedEffRatioVsTime[mm][nn]->SetMarkerColor(nn+1);
        gLedEffRatioVsTime[mm][nn]->SetLineColor(nn+1);
        sprintf(name,"PMTS%d_vs_PMTS%d",temp1index[mm],temp2index[nn]);
        legLedEvtRatio[mm]->AddEntry(gLedEffRatioVsTime[mm][nn],name,"l");
      }
      legLedEvtRatio[mm]->Draw();
      cLedEvtRatio[mm]->Write(); cLedEvtRatio[mm]->Close();
    }
*/
    /*
    TCanvas* cLedEvtRatioTot = new TCanvas(); 
    cLedEvtRatioTot->SetName("cLedEventRate_top_over_bot");
    gLedEffRatioVsTimeTot->Draw("AP");
    gLedEffRatioVsTimeTot->SetMarkerStyle(8);
    gLedEffRatioVsTimeTot->SetMarkerSize(0.7);
    gLedEffRatioVsTimeTot->GetXaxis()->SetTimeDisplay(1);
    gLedEffRatioVsTimeTot->GetXaxis()->SetTimeOffset(t00.Convert(0));
    gLedEffRatioVsTimeTot->GetXaxis()->SetTimeFormat("%m/%d");
    gLedEffRatioVsTimeTot->GetYaxis()->SetTitle("Event ratio");    
    cLedEvtRatioTot->Write(); cLedEvtRatioTot->Close();

// write TDC
//    dirTDC->cd();
    for(int i=0;i<16;i++){
        for(int j=0;j<nbTrig;j++){
            //dir[j]->cd();
            dir_TDC[j]->cd();
            hTDCnew[i][j]->Write();
	    hTDC[i][j]->Write();
            //hTDCnewDelta[i][j]->Write();
        }
    }
    */
    //for(int i=0; i<6;i++) hTimeDiff_hodoDet_hodoTrig[i]->Write();
    //for(int i=0; i<8;i++) hTimeDiff_PMTs_hodoTrig[i]->Write();
    //hTDC_H0H2->Write();
    //hTDC_H1H2->Write();
    //hTDC_H3H2->Write();
    //hTDC_H1H0->Write();
    //hTDC_H3H0->Write();
    //hTDC_H5H2->Write();
    //hTDC_H4H2->Write();
	
     //cout<<" -- 1"<<endl;
/*
    TCanvas* cTDC = new TCanvas(); cTDC->SetName("cAllTDCChannels");
    TLegend* legTDC = new TLegend(0.6,0.7,0.9,0.9);
    for(int i=0;i<16; i++) {
        hTDC[i]->Write();
        hTDC[i]->SetLineWidth(2);
        if(i<=7){hTDC[i]->SetLineColor(i+1); hTDC[i]->SetLineStyle(2);}
        else {hTDC[i]->SetLineColor(i+1-8); hTDC[i]->SetLineStyle(1);}
        cTDC->cd();
        if(i==0) hTDC[i]->Draw("");
        else  hTDC[i]->Draw("same");
        legTDC->AddEntry(hTDC[i],hTDC[i]->GetName(),"l");
    }
    legTDC->Draw();
    cTDC->Update(); cTDC->Write();  cTDC->Close();
    TCanvas* cTDC1=new TCanvas(); cTDC1->SetName("cChannelOfInterest");
    hTDC[6]->Draw(); hTDC[7]->Draw("same"); hTDC[8]->Draw("same");
    TLegend* legTDC1 = new TLegend(0.6,0.7,0.9,0.9);
    legTDC1->AddEntry(hTDC[6],hTDC[6]->GetName(),"l");
    legTDC1->AddEntry(hTDC[7],hTDC[7]->GetName(),"l");
    legTDC1->AddEntry(hTDC[8],hTDC[8]->GetName(),"l");
    legTDC1->Draw();
    cTDC1->Write(); cTDC1->Close();
*/
/*
    TCanvas* cTDC2 = new TCanvas(); cTDC2->SetName("cAllTDCChannels_Relative");
    TLegend* legTDC2 = new TLegend(0.6,0.7,0.9,0.9);
    for(int i=0;i<16; i++) {
        hTDC1[i]->Write();
        hTDC1[i]->SetLineWidth(2);
        if(i<=7){hTDC1[i]->SetLineColor(i+1); hTDC1[i]->SetLineStyle(2);}
        else {hTDC1[i]->SetLineColor(i+1-8); hTDC1[i]->SetLineStyle(1);}
        cTDC2->cd();
        if(i==0) hTDC1[i]->Draw("");
        else  hTDC1[i]->Draw("same");
        legTDC2->AddEntry(hTDC1[i],hTDC1[i]->GetName(),"l");
    }
    legTDC2->Draw();
    cTDC2->Update(); cTDC2->Write();  cTDC2->Close();
    TCanvas* cTDC3=new TCanvas(); cTDC3->SetName("cChannelOfInterest_Relative");
    hTDC1[6]->Draw(); hTDC1[7]->Draw("same"); hTDC1[8]->Draw("same");
    TLegend* legTDC3 = new TLegend(0.6,0.7,0.9,0.9);
    legTDC3->AddEntry(hTDC1[6],hTDC1[6]->GetName(),"l");
    legTDC3->AddEntry(hTDC1[7],hTDC1[7]->GetName(),"l");
    legTDC3->AddEntry(hTDC1[8],hTDC1[8]->GetName(),"l");
    legTDC3->Draw();
    cTDC3->Write(); cTDC3->Close();
*/

// write ADC, QDC ..
    /*
    for(int i=0;i<8;i++){
        for(int j=0;j<nbTrig;j++){
            dir[j]->cd();
            dir_QDC[j]->cd();
            if(strcmp(sQDC1ChMap[i].c_str(),"")!=0) hQDC1[i][j]->Write();
            if(strcmp(sQDC2ChMap[i].c_str(),"")!=0) hQDC2[i][j]->Write();
        }
        //if(i<6) hDigiSTVsTDC_CosmicTrig[i]->Write();
    }
    */

// write new digitizer plots
//   dirDigitizerNew->cd();

   for(int i=0;i<pmtnb;i++){
         //dir_Digitizer[2]->cd(); 
         //hScat_hodoTDC_vs_PMTT[i]->Write();
         /*
      for(int j=0;j<nbTrig;j++){
         dir[j]->cd();
         dir_Digitizer[j]->cd();
         //hNumberOfPulses[i][j]->Write();
         //hDigiAmp[i][j]->Write();
         //hDigiAmpTimeBin[i][j]->Write();
         //hScat_DigiAmp_vs_TimeBin[i][j]->Write();
         hDigitizerPulseCharge_npe[i][j]->Write();
         hDigitizerPulseCharge1[i][j]->Write();
         hDigiPulseAmp[i][j]->Write();
         //hDigitizerChargeRatio[i][j]->Write();
         hDigitizerPulseCharge_vs_PulseTime[i][j]->Write();
         //hDigiPulseChargeNpe_vs_PulseWidth[i][j]->Write();
         //hDigitizerPulseWidth_vs_PulseTime[i][j]->Write();
         hDigitizerPulse_TimeDifference[i][j]->Write();
         hAfterPulseChargeCorrelation[i][j]->Write();
         //hStopMuonTime[i][j]->Write();
         //hDigiChargeNpe_vs_DigiAmp[i][j]->Write();
         hDigiST[i][j]->Write();
         hDigiST_all[i][j]->Write();
         //hDigiPW[i][j]->Write();

         if(i==0){
           hDigitizerPulseChargeSum[j]->Write();
           hDigitizerPulseChargeSumTop[j]->Write();
           hDigitizerPulseChargeSumBottom[j]->Write();
           hDigitizerPulseChargeSumCorrelation[j]->Write(); // all triggers
           hDigitizerPulseChargeSum1[j]->Write();
           hDigitizerPulseChargeSumTop1[j]->Write();
           hDigitizerPulseChargeSumBottom1[j]->Write();
           //hDigitizerPulseChargeSumS2S3[j]->Write();
           //hNHit[j]->Write();
           hTimeDiff[j]->Write();
           //hTimeDiff1[j]->Write();
           hTimeDiffSelect[j]->Write();
           hTimeDiffSelect1[j]->Write();
           hDigitizerPulseChargeDecay[j]->Write();
           hDigitizerPulseChargeDecay1[j]->Write();
           hDecayChargeVsTime[j]->Write();
           hDecayChargeVsTime1[j]->Write();
           hChargeRatioVsTimeInDecay[j]->Write();
           hChargeRatioVsTimeInDecay1[j]->Write();
           //hScat_DigiAmpTS0_vs_DigiAmpTS1[j]->Write();
           //hScat_TDCS0_vs_TDCS1[j]->Write();
           hTotalNumberOfPulses[j]->Write();
           hCoincidence[j]->Write();
           hCoincidenceVsDecayTime[j]->Write();
           hCoincidenceNumbers[j]->Write();
         }
         if(i<8) {
           hScat_TDC_vs_PMTT[i][j]->Write();
           if(i<6) hScat_TDC_vs_PMTT1[i][j]->Write();
           hTDC_minus_PMTT[i][j]->Write();
         }
         if(i>0) {
           //hScat_NPulses_vs_S0[i-1][j]->Write();
           //hScat_PulseTime_vs_S0[i-1][j]->Write();
           //if(i>1) hScat_PulseTime_vs_S1[i-2][j]->Write();
           //if(i>2) hScat_PulseTime_vs_S2[i-3][j]->Write();
           //if(i>3) hScat_PulseTime_vs_S3[i-4][j]->Write();
           hPulseTimeDiffToS0[i-1][j]->Write();
           if(i>1) hPulseTimeDiffToS1[i-2][j]->Write();
           if(i>2) hPulseTimeDiffToS2[i-3][j]->Write();
           if(i>3) hPulseTimeDiffToS3[i-4][j]->Write();
           if(i>4) hPulseTimeDiffToS4[i-5][j]->Write();
           if(i>5) hPulseTimeDiffToS5[i-6][j]->Write();
           if(i>6) hPulseTimeDiffToS6[i-7][j]->Write();
         }
      }
      */
      //dirPedestal->cd();
      //hPedestalMean[i]->Write();
      //hPedestalWidth[i]->Write();
      //hWaveformThreshold[i]->Write();
   }

// write evtTimeDiff histograms
//   dirEvtTimeDiff->cd();
   //for(int i=0;i<nbTrig;i++){
   //  dir_evtTimeDiff[i]->cd();
   //  hEvtTimeDiff[i]->Write();
   //}
    /*
    // write event  rate (scaler ) 
    dirEvtRate->cd();
    //for(int i=0;i<17;i++) hEvtRate[i]->Write();
    TCanvas* cEvtRate = new TCanvas("cEvtRate","",800,600);
    cEvtRate->Divide(1,2); 
    cEvtRate->cd(1);
    gEvtRateVsTime[2]->Draw("AL"); gEvtRateVsTime[2]->SetLineColor(2); 
    gEvtRateVsTime[2]->GetXaxis()->SetTitle("Time");
    gEvtRateVsTime[2]->GetYaxis()->SetTitle("Average rate (Hz)");
    gEvtRateVsTime[2]->GetXaxis()->SetTimeDisplay(1);
    gEvtRateVsTime[2]->GetXaxis()->SetTimeOffset(t00.Convert(0));
    gEvtRateVsTime[2]->GetXaxis()->SetTimeFormat("#splitline{%H:%M}{%m/%d}");
    gEvtRateVsTime[2]->GetXaxis()->SetLabelOffset(0.02);
    gEvtRateVsTime[2]->GetXaxis()->SetTitleOffset(1.5);
    cEvtRate->cd(2);
    gEvtRateVsTime[1]->Draw("AL"); gEvtRateVsTime[1]->SetLineColor(1);
    gEvtRateVsTime[0]->Draw("L"); gEvtRateVsTime[0]->SetLineColor(4);
    gEvtRateVsTime[1]->GetXaxis()->SetTitle("Time");
    gEvtRateVsTime[1]->GetYaxis()->SetTitle("Average rate (Hz)");
    gEvtRateVsTime[1]->GetXaxis()->SetTimeDisplay(1);
    gEvtRateVsTime[1]->GetXaxis()->SetTimeOffset(t00.Convert(0));
    gEvtRateVsTime[1]->GetXaxis()->SetTimeFormat("#splitline{%H:%M}{%m/%d}");
    gEvtRateVsTime[1]->GetXaxis()->SetLabelOffset(0.02);
    gEvtRateVsTime[1]->GetXaxis()->SetTitleOffset(1.5);
    cEvtRate->cd(1);
    TLegend* legEvtRate = new TLegend(0.6,0.7,0.9,0.9);
    legEvtRate->AddEntry(gEvtRateVsTime[0],"Led triggers","l");
    legEvtRate->AddEntry(gEvtRateVsTime[1],"Hodo triggers","l");
    legEvtRate->AddEntry(gEvtRateVsTime[2],"Multi triggers","l");
    legEvtRate->Draw();
    cEvtRate->Write(); cEvtRate->Close();

    hLevelVsMultiEvtRate->Write();
    hLevelVsTemperature->Write();
    //hHodoRateVsMultiRate->Write();
    //hTestHodoSelection->Write();
	//hTDChodo_vs_H0->Write();
	//hTDChodo_vs_H1->Write();
	//hTDChodo_vs_H2->Write();
	//hTDChodo_vs_H3->Write();
	//hTDChodo_vs_H4->Write();
	//hTDChodo_vs_H5->Write();
	//hTDChodo_vs_multi->Write();
	//hTDChodo_vs_Led->Write();
    */
    for(int i=0; i<pmtnb; i++){
        for(int j=0; j<nbTrig; j++){
            dir_TrigType[j]->cd();
            hCharge[i][j]->Write();
            hCharge1[i][j]->Write();
        }
    }
    for(int i=0; i<pmtnb; i++){
        for(int j=0; j<nbTrig; j++){
            dir_TrigType[j]->cd();
            hDigiST[i][j]->Write();
        }
    }
    for(int i=0; i<pmtnb; i++){
        for(int j=0; j<nbTrig; j++){
            dir_TrigType[j]->cd();
            hDigiPulseAmp[i][j]->Write();
        }
    }        

    rtfile->Close();
    /*
    cout<<"trigger counts (all, multi, hodo, led, multi_or_hodo, multi_and_hodo, stopped muon candidate, strict_hodo):"<<endl;
    cout<<trigcnts_all<<"\t"<<trigcnts_multi<<"\t"<<trigcnts_hodo<<"\t"
        <<trigcnts_led<<"\t"<<trigcnts_multi_or_hodo<<"\t"<<trigcnts_multi_and_hodo<<"\t"
        <<trigcnts_stoppedmuon<<"\t"<<trigcnts_stricthodo
        <<endl;
	cout<<"trigger: hodo + H5 + not H4 "<<trigcnts_hodo_H5_notH4<<endl;
	cout<<"trigger: hodo + not H5 + not H4 "<<trigcnts_hodo_notH5_notH4<<endl;
	cout<<"trigger: hodo + not H5 + H4 "<<trigcnts_hodo_notH5_H4<<endl;
	cout<<"trigger: hodo + H5 + H4 "<<trigcnts_hodo_H5_H4<<" (should be 0) "<<endl;
    */
    for(int index=0; index<nbTrig; index++){
        cout<<"NbOfEvents for "<<sTriggerType[index]<<": \t"<<triggercnt[index]<<endl;
    }

    //cout<<"test cnts 2 coincidences "<<test_cnt1<<"\t 3 coincidences "<<test_cnt2<<endl;
	
    //	cout<<"selectHodo "<<selectHodo<<"\t"<<"selectHs "<<selectHs<<endl;
	
    //	cout<<"test number of events in TDC_hodo!=-1: <<"<<selectHodo<<"\n"
    //   <<"                         TDC_hodo!=-1 and TDC_multi!=-1: "<<testcnt1<<"\n"
    //   <<"                         TDC_hodo!=-1 and TDC_led!=-1:   "<<testcnt2<<"\n";

  return 0;
}


