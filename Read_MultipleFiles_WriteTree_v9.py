import h5py
import numpy as np
#import matplotlib.pyplot as plt
#import graphUtils
import ROOT
from ROOT import TH1F, TTree, TFile, TSpectrum, TCanvas, TLine, TDirectory
import datetime
import time
import timeit
from array import array
import zipfile
import os
import sys

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
#
#  Version v7 (20180719)
#this programs is to be run as the following
#     python *.py *.zip
# each zip file contains the .h5 data file
#it is designed to be run as jobs at the racf cluster
####

#
#some initial numbers
#two QDCs (CAEN V965A), each has 8 channels
NbQDC1CH = 8  #QDC1 CHs 0 1 2 3 4 5 are used for H0, 1, 2, 3, 4, 5. now just include unused channels in tree
NbQDC2CH = 8  #QDC2 CHs 0 6 2 3 4 5 are used for S0, 1, 2, 3, 4, 5.
#TDC (CAEN V775) has 32 channels but we use 16 channels
NbTDCCH = 16
#Scaler (CAEN V830) has 32 channels but we use 17 channels
NbScalerCH = 17
#two digitizers (CAEN V1729A), each has 4 channels
NbDigiCH = 8
NSamples = 2560 # Number of samples for each digitizer channel, runs at 1GS/s so sampling point interval 1 ns.
SaveWaveforms = 0 # a bool value to indicate whether or not save waveforms: 0 for not save and 1 for save

#prepare file names to be read
rawdatapath = "./"
hdf5_file_zip_name = []
hdf5_file_zip_name.append(sys.argv[1])
directory_to_extract_to = "./"

hdf5_file_name = []
for aZipfile in hdf5_file_zip_name:
    aname = rawdatapath + aZipfile
    zip_ref = zipfile.ZipFile(aname, 'r')
    zip_ref.extractall(directory_to_extract_to)
    zip_ref.close()
    h5filename = str.replace(aZipfile,".zip",".h5")
    hdf5_file_name.append(h5filename)    

start_time = timeit.default_timer()

# start reading the files
totalEvtNumber = 0
for filename in hdf5_file_name:
    #prepare histograms and ROOT file
    #for each .h5 file create a ROOT Tree
    rootfilename = ""
    if SaveWaveforms==0:
        rootfilename = str.replace(filename,".h5","_Tree_v8.root")
    else:
        rootfilename = str.replace(filename,".h5","_Tree_v8_waveformsSaved.root")
    rootfile = TFile(rootfilename,"recreate")
    print rootfilename
    print filename
    if SaveWaveforms == 1: # save waveforms in this ROOT file too
        waveformsDir = rootfile.mkdir("waveformsDir")
    waveformsAvgDir = rootfile.mkdir("waveforms_Avg_Dir")
    mytree = TTree( 'OneTonEvent', 'My Event Tree' )
    #the variables to be stored in the ROOT Tree
    EvtNumber = array('i', [0]) #event number, a signle integer
    EvtTime = array('d',[0])    #event time, it is a timestamp, I use a double to hold it
    TemperatureBox = array('f',[0.])  #event temperature in the dark room, a float (this temperature sensor, NI-TC01, is not calibrated, unit in celsius degrees)
    TemperatureRack = array('f',[0.]) #event temperature at the daq rack, a float, unit in celsius degrees
    Resistivity = array('f',[0.])     #event resistivity, a float, unit MOhm*cm
    LiquidLevel1 = array('f',[0.])    #event liquid level 1 (inlet) from the old float sensor
    LiquidLevel2 = array('f',[0.])    #event liquid level 2 (outlet) from the old float sensor
#    usLiquidLevel1 = array('f',[0.])  #event liquid level 1 (measured by ultrasonic level sensor), need it later
#    usLiquidLevel2 = array('f',[0.])  #event liquid level 2 (measured by ultrasonic level sensor), need it later
    QDC1 = array('f',NbQDC1CH*[0.])   #event QDC1, array, size 8, float
    QDC2 = array('f',NbQDC2CH*[0.])   #event QDC2, array, size 8,float
    TDC = array('f',NbTDCCH*[0.])     #event TDC, array, size 16, float
    Scaler = array('f',NbScalerCH*[0.]) #event Scaler, array, size 17, float
    TrigType = array('i',[0]) # I define 1 -> multiplicity, 2 -> Cosmic and 3->LED triggers
    TrigTypeFlag_Multi = array('i',[0]) #trigger type multipliticy, integer numbers, another way of determining trigger types
    TrigTypeFlag_Hodo = array('i',[0])  #trigger type hodoscope
    TrigTypeFlag_Led = array('i',[0])   #trigger type led
    #digitizer variables                         #in my old method, I treat every waveform either one or zero pulse, now I want to count number of pulses in each waveform!
    DigitizerCharge = array('f',NbDigiCH*[0.])   #digitizer charge, array, size 8, float
    DigitizerAmplitude = array('f',NbDigiCH*[0.])#digitizer amplitude, array, size 8, float
    DigitizerAmplitudeTimeBin = array('f',NbDigiCH*[0.])#digitizer amplitude time bin, array, size 8, float
    DigitizerStartTimeBin = array('f',NbDigiCH*[0.])#digitizer start time bin
    DigitizerRiseTime =  array('f',NbDigiCH*[0.]) #digitizer rise time
    DigitizerPulseWidth = array('f',NbDigiCH*[0.])#digitizer pulse width
    DigitizerFullPulseWidth = array('f',NbDigiCH*[0.])
    #digitizer pedestal
    DigitizerPedMean = array('f',NbDigiCH*[0.]) #array of 8, float
    DigitizerPedWidth = array('f',NbDigiCH*[0.])
    #digitizer, new method, counting number of pulses
    DigitizerNumberOfPulses = array('i',NbDigiCH*[0]) #in new method, I need to count the number of pulses in each waveform.
    #digitizer charge, start bin, end bin
    # below the three variables are 2D arrays. 
    DigitizerPulseCharge = np.array([ [0.0]*20 for i in range(8)], np.dtype('d')) # 8 PMT channels, each assume maximum 20 pulses, this is a 8 by 20 array, for more than 20 pulses, only save the first 20 pulses the event. The final number of pulses is stored for later easier process.
    DigitizerPulseStartBin = np.array([ [0.0]*20 for i in range(8)], np.dtype('d')) # 8 PMT channel, each assume maximum 20 pulses, this is a 8 by 20 array. stores start bin of each found pulse
    DigitizerPulseEndBin = np.array([ [0.0]*20 for i in range(8)], np.dtype('d')) # 8 PMT channel, each assume maximum 20 pulses, this is a 8 by 20 array. stores end bin of each found pulse
    DigitizerPulseAmplitude = np.array([ [0.0]*20 for i in range(8)], np.dtype('d')) # 8 PMT channel, each assume maximum 20 pulses, this is a 8 by 20 array. stores amplitude of each found pulse
    DigitizerPulseAmplitudeTimeBin = np.array([ [0]*20 for i in range(8)], np.dtype('i')) # 8 PMT channel, each assume maximum 20 pulses, this is a 8 by 20 array. stores time bin at amplitude of each found pulse
    MuonDecayTime = array('i',NbDigiCH*[0]) #this info is not checked, could be wrong. But leave it here for now.
    #
    mytree.Branch("EvtNumber",EvtNumber,"EvtNumber/I")
    mytree.Branch("EvtTime",EvtTime,"EvtTime/D")
    mytree.Branch("TemperatureBox",TemperatureBox,"TemperatureBox/F")
    mytree.Branch("TemperatureRack",TemperatureRack,"TemperatureRack/F")
    mytree.Branch("Resistivity",Resistivity,"Resistivity/F")
    mytree.Branch("LiquidLevel1",LiquidLevel1,"LiquidLevel1/F")
    mytree.Branch("LiquidLevel2",LiquidLevel2,"LiquidLevel2/F")
#    mytree.Branch("usLiquidLevel1",usLiquidLevel1,"usLiquidLevel1/F") # to be used later
#    mytree.Branch("usLiquidLevel2",usLiquidLevel2,"usLiquidLevel2/F")
    mytree.Branch("QDC1",QDC1,"QDC1[8]/F")
    mytree.Branch("QDC2",QDC2,"QDC2[8]/F")
    mytree.Branch("TDC",TDC,"TDC[16]/F")
    mytree.Branch("Scaler",Scaler,"Scaler[17]/F")
    mytree.Branch("DigitizerCharge",DigitizerCharge,"DigitizerCharge[8]/F")
    mytree.Branch("DigitizerAmplitude",DigitizerAmplitude,"DigitizerAmplitude[8]/F")
    mytree.Branch("DigitizerAmplitudeTimeBin",DigitizerAmplitudeTimeBin,"DigitizerAmplitudeTimeBin[8]/F")
    mytree.Branch("DigitizerStartTimeBin",DigitizerStartTimeBin,"DigitizerStartTimeBin[8]/F")
    mytree.Branch("DigitizerRiseTime",DigitizerRiseTime,"DigitizerRiseTime[8]/F")
    mytree.Branch("DigitizerPulseWidth",DigitizerPulseWidth,"DigitizerPulseWidth[8]/F")
    mytree.Branch("DigitizerFullPulseWidth",DigitizerFullPulseWidth,"DigitizerFullPulseWidth[8]/F")
    mytree.Branch("DigitizerPedMean",DigitizerPedMean,"DigitizerPedMean[8]/F")
    mytree.Branch("DigitizerPedWidth",DigitizerPedWidth,"DigitizerPedWidth[8]/F")
    mytree.Branch("TrigType",TrigType,"TrigType/I")
    mytree.Branch("TrigTypeFlag_Multi",TrigTypeFlag_Multi,"TrigTypeFlag_Multi/I")
    mytree.Branch("TrigTypeFlag_Hodo",TrigTypeFlag_Hodo,"TrigTypeFlag_Hodo/I")
    mytree.Branch("TrigTypeFlag_Led",TrigTypeFlag_Led,"TrigTypeFlag_Led/I")
    # new method, to store digitizer number of pulses and charges etc.
    mytree.Branch("DigitizerNumberOfPulses",DigitizerNumberOfPulses,"DigitizerNumberOfPulses[8]/I")
    mytree.Branch('DigitizerPulseCharge', DigitizerPulseCharge, 'DigitizerPulseCharge[8][20]/D') # 8 PMT channels, each assume maximum 20 pulses
    mytree.Branch('DigitizerPulseStartBin', DigitizerPulseStartBin, 'DigitizerPulseStartBin[8][20]/D') # 8 PMT channel, each assume maximum 20 pulses
    mytree.Branch('DigitizerPulseEndBin', DigitizerPulseEndBin, 'DigitizerPulseEndBin[8][20]/D') # 8 PMT channel, each assume maximum 20 pulses
    mytree.Branch('DigitizerPulseAmplitude', DigitizerPulseAmplitude, 'DigitizerPulseAmplitude[8][20]/D') # 8 PMT channels, each assume maximum 20 pulses
    mytree.Branch('DigitizerPulseAmplitudeTimeBin', DigitizerPulseAmplitudeTimeBin, 'DigitizerPulseAmplitudeTimeBin[8][20]/I') # 8 PMT channels, each assume maximum 20 pulses
    mytree.Branch("MuonDecayTime",MuonDecayTime,"MuonDecayTime[8]/I")
    #
    #okay, now start reading the data file
    file = h5py.File(filename,'r') # open one h5 file
    #CalData1 = file['Calibration']['Digitizer_1']['Pedestal']
    #CalData2 = file['Calibration']['Digitizer_2']['Pedestal']
    EvtDir = file['Events'] # get the Events card
    #timeCnt = 0;
    initCnts=[] #initial counts (the first event) in scaler, all active channels
    initTime=0 #initial time (the first event)
    # average waveforms for led, hodo and multiplicity triggers
    # I will use those signal pulses for average (skip if no pulse in the signal)
    cnt_evt_ledTrig = cnt_evt_hodoTrig = cnt_evt_multiTrig = 0
    #if 1: #SaveWaveforms == 1: record average waveforms for each PMT no matter what, this doesnot take too much time
    hWaveforms_avg = []
    totalwave = []
    totalwave_cnts = [0 for i in range(0,24)]
    channelname = ["S0","S1","S2","S3","S4","S5","S6","S7"]
    trigFlagname = ["LedTrig","HodoTrig","MultiTrig"]
    for channel in range(0,8): 
        for trigFlag in range(0,3):
            histname = "AvgWaveform_"+channelname[channel]+"_"+trigFlagname[trigFlag]
            onehist = TH1F(histname,"",2560,-0.5,2560-0.5)
            hWaveforms_avg.append(onehist)
            totalwave.append([0.0 for i in range(0,2560)])
    # start reading in events
    for evtkey in EvtDir.keys():
        if not(evtkey=='evt_table'):
            EvtNumber[0] = int(evtkey)
#            if EvtNumber[0]>20:# and EvtNumber[0]<1900:
#                break;
            totalEvtNumber += 1
            thisDigi1 = EvtDir[evtkey]["Digitizer_1"]
            thisDigi2 = EvtDir[evtkey]["Digitizer_2"]
            thisQDC1 = EvtDir[evtkey]["QDC_1"]
            thisQDC2 = EvtDir[evtkey]["QDC_2"]
            thisTDC  = EvtDir[evtkey]["TDC"]
            thisScaler = EvtDir[evtkey]["Scaler"]
            thisTime = EvtDir[evtkey]['Event_Time']
            thisResistivity = EvtDir[evtkey]['Event_Resistivity']
            thisTemperature = EvtDir[evtkey]['Event_Temp']
            thisTempeAtRack = EvtDir[evtkey]['Event_TempRack']
	    LiquidLevel1[0] = EvtDir[evtkey]['Event_LevelSensor'][0]
	    LiquidLevel2[0] = EvtDir[evtkey]['Event_LevelSensor'][1]
	    #print thisLevel1, thisLevel2
            EvtTime[0]=thisTime[()]
            TemperatureBox[0]=thisTemperature[()]
            TemperatureRack[0]=thisTempeAtRack[()]
            Resistivity[0]=thisResistivity[()]
            # determine the trigger type using TDC information
            # I will use trigFlag for MultiTrig (1), CosmicTrig (2) and LedTrig (3)
            trigFlag = 0
            trig_temp1, trig_temp2, trig_temp3, thistrig = 0, 0, 0, 1 # this is the old method
            TrigFlag_Multi, TrigFlag_Hodo, TrigFlag_Led = 0, 0, 0 # this is the new method, it is (tested) same as the old method
            for ch, num in thisTDC:
                if int(ch) == 7: # multiplicity
                    trig_temp1 = 1
                    TrigFlag_Multi = 1
                if int(ch) == 6: # cosmic
                    trig_temp2 = 2
                    TrigFlag_Hodo = 1
                if int(ch) == 8: # led trigger
                    trig_temp3 = 3
                    TrigFlag_Led = 1
                if int(ch) == 9: # all trigger, should be always 1, otherwise there is a problem
                    thistrig = 1
            #well, there could be other trigger types, for example, both multiTrig and cosmic trigger
            #but I will first focus on the first three types: multi, cosmic, and led.
            if thistrig==1 and trig_temp1 == 1 and trig_temp2 != 2 and trig_temp3 != 3:
                trigFlag = 1 #multiplicity trigger only
                cnt_evt_multiTrig += 1
            if thistrig==1 and trig_temp1 != 1 and trig_temp2 == 2 and trig_temp3 != 3:
                trigFlag = 2 #cosmic trigger only
                cnt_evt_hodoTrig += 1
            if thistrig==1 and trig_temp1 != 1 and trig_temp2 != 2 and trig_temp3 == 3:
                trigFlag = 3 #LED trigger only
                cnt_evt_ledTrig += 1
            if thistrig==1 and trig_temp1 == 1 and trig_temp2 == 2 and trig_temp3 !=3:
                trigFlag = 4 #both multi and cosmic triggers, not led trigger
            #if cnt_evt_multiTrig>50 and cnt_evt_hodoTrig>50 and cnt_evt_ledTrig>50:
            #    break
            TrigType[0] = trigFlag
            TrigTypeFlag_Multi[0] = TrigFlag_Multi
            TrigTypeFlag_Hodo[0] = TrigFlag_Hodo
            TrigTypeFlag_Led[0] = TrigFlag_Led
            #process QDC data
            for ch, num, c, d in thisQDC1:
                QDC1[int(ch)] = num
            for ch, num, c, d in thisQDC2:
                QDC2[int(ch)] = num
            #process Scaler data
            i=0
            if int(evtkey) == 1: # for the first event, set rate to zero
                for counts in thisScaler:
                    initTime = thisTime[()]
                    initCnts.append(counts)
                    Scaler[i] = 0
                    i += 1
            else: # for all other events, use the counts between this event and last event to determine a rate
                for counts in thisScaler:
                    rate = (counts - initCnts[i])/(thisTime[()] - initTime)
                    Scaler[i] = rate  
                    i += 1
            #process TDC data
            # Before assigning values, clear the TDC array.
            # intentatively, I set the array to have -1. Later when processing results, if I see -1, that means no entry for the channel in the event.
            for ch in range(0,16,1):
                TDC[ch] = -1.0
            for ch, num in thisTDC: # cont. in the data, only TDC channels with entries are stored.
                TDC[int(ch)] = num
            #process digitizer data
            digi1_ped_mean = [] # pedestal for the 4 channels in this digitizer
            digi1_threshold = []
            digi2_ped_mean = []
            digi2_threshold = []
            mean, mean2 = 0, 0
            sigma, sigma2 = 0, 0
            if np.mod(int(evtkey),1) == 0:
                if np.mod(int(evtkey),100) ==0:    print evtkey, totalEvtNumber
                #prepare for digitizer pedestals
                #get pedestals digitizer 1
                for ch in range(0,4,1):
                    # use bins 1500-2560, and reduce bias by removing some possible signals in that range,
                    # I assume baseline should be less than +-10 ADC units, this is decided after looking at the processed pedestals
                    subwavelist = [] # make a sub waveform list to hold those bins for pedestals
                    for bin in range(1000,1800,1): # I was using a range 1500-2560, but for S2 some waveforms have a upwarp at the tail, to avoid that I change to 1000-1800
                        if abs(thisDigi1[ch][bin])<10:
                            subwavelist.append(thisDigi1[ch][bin])
                    mean = np.mean(subwavelist)
                    sigma = np.std(subwavelist)
                    del subwavelist
                    digi1_ped_mean.append(mean)
                    digi1_threshold.append(mean-5.0*sigma)
                    DigitizerPedMean[ch] = mean
                    DigitizerPedWidth[ch] = sigma
                #get pedestals digitizer 2
                for ch in range(0,4,1):
                    subwavelist = []
                    for bin in range(1000,1800,1):
                        if abs(thisDigi2[ch][bin])<10:
                            subwavelist.append(thisDigi2[ch][bin])
                    mean2 = np.mean(subwavelist)
                    sigma2 = np.std(subwavelist)
                    del subwavelist
                    digi2_ped_mean.append(mean2)
                    digi2_threshold.append(mean2-5.0*sigma2)
                    DigitizerPedMean[ch+4] = mean2
                    DigitizerPedWidth[ch+4] = sigma2
                #prepare for digitizer waveform pulse finding and charge integration
                #initialize the 2d arrays that hold pmt pulse charge, startbin,...
                for pmtNNN in range(8): # because I have 8 PMTs
                    for pulseN in range(20): # each waveform I assume maximum 20 pulses
                        DigitizerPulseCharge[pmtNNN][pulseN] = 0.0
                        DigitizerPulseStartBin[pmtNNN][pulseN] = 0.0
                        DigitizerPulseEndBin[pmtNNN][pulseN] = 0.0
                waveformname_head = "Evt_" # used for waveforms
                # for digitizer 1
                for ch in range(0,4,1):
                    #begin my own algorithm for pulse finding
                    #first make pulse islands by zero suppression -- assign zeros to all those bins that are below my threshold
                    thiswaveformcopy = [0 if x>digi1_threshold[ch] else x for x in thisDigi1[ch]]
                    #maximum = thisDigi1[ch][np.argmin(thisDigi1[ch][80:])+80] # get the amplitude of the pulse. note: don't consider the first 80 bins because of some weird small pulse in the beginning of the waveform (between bin 20 and 50)!!
                    DigitizerAmplitudeTimeBin[ch] = np.argmin(thisDigi1[ch][80:])
                    maximum = thisDigi1[ch][int(80+DigitizerAmplitudeTimeBin[ch])]
                    DigitizerAmplitude[ch] = -1.0*(maximum-digi1_ped_mean[ch])
                    #here do avearage waveform, sum waveforms first
                    if SaveWaveforms == 0 or SaveWaveforms == 1:
                        if TrigFlag_Led == 1 and TrigFlag_Hodo == 0 and TrigFlag_Multi ==0:# and totalwave_cnts[ch*3+0]<=50:
                            totalwave_cnts[ch*3+0] += 1
                            for bin in range(0,2560,1):
                                totalwave[ch*3+0][bin] += thisDigi1[ch][bin]
                        if TrigFlag_Led == 0 and TrigFlag_Hodo == 1 and TrigFlag_Multi ==0:# and totalwave_cnts[ch*3+1]<=50:
                            totalwave_cnts[ch*3+1] += 1
                            for bin in range(0,2560,1):
                                totalwave[ch*3+1][bin] += thisDigi1[ch][bin]
                        if TrigFlag_Led == 0 and TrigFlag_Hodo == 0 and TrigFlag_Multi ==1:# and totalwave_cnts[ch*3+2]<=50:
                            totalwave_cnts[ch*3+2] += 1
                            for bin in range(0,2560,1):
                                totalwave[ch*3+2][bin] += thisDigi1[ch][bin]
                    nPulses = 0
                    flagStartFound = 0
                    #if no signal (maximum is below threshold), no need to do anything, but fill nPulses=0 is needed
                    #after looking futher, I realize some pulses have amplitude surpass threshold, but the amplitude is really low (smaller than 1 p.e.), I choose to cut on amplitude at 10 ADC counts
                    if DigitizerAmplitude[ch] <= 10: # if maximum >= digi1_threshold[ch]:
                        DigitizerNumberOfPulses[ch] = nPulses
                        #save waveforms if needed
                        if SaveWaveforms == 1:
                            waveformname = waveformname_head + str(int(evtkey)) + "_trigFlag_" + str(TrigFlag_Multi) + str(TrigFlag_Hodo) + str(TrigFlag_Led) + "_ch_" + str(ch) + "_nPulses_" + str(nPulses)
                            waveformtitle = waveformname + "_mean_" + str(digi1_ped_mean[ch]) + "_sigma_" + str(sigma)
                            hist = TH1F(waveformname,waveformtitle,2560,-0.5,2560-0.5)
                            for bin in range(0,2560,1):
                                hist.SetBinContent(bin+1,thisDigi1[ch][bin])
                            waveformsDir.cd()
                            hist.Write()
                            del hist
                        continue
                    #otherwise at least there is a pulse 
                    bin=80
                    while bin>=80 and bin<(2560-4):
                        if nPulses==20: #maximum number of pulses, if 20 pulses are found then stop the searching
                            break
                        #find the start bin of a pulse in a waveform
                        #To define a cross over event: two consecutive bins below threshold and following two bins above threshold
                        #Okay, begin my searching. Find the first signal rising cross-threshold
                        if flagStartFound==0:
                            if thiswaveformcopy[bin]<0:
                                #if the condition is true, then it crosses threshold between bin+1 and bin+2
                                DigitizerPulseStartBin[ch][nPulses] = bin - 3
                                flagStartFound = 1
                                bin+=3
                                continue
                        #then, find the signal falling cross-threshold
                        if flagStartFound==1:
                            if thiswaveformcopy[bin]==0:
                                # if the above condition is ture, the cross is between bin+1 and bin+2
                                DigitizerPulseEndBin[ch][nPulses] = bin+3
                                nPulses=nPulses+1
                                flagStartFound = 0
                        bin+=1
                    #Combine neighbor pulses as one if they are separated <= 25 ns.
                    peakxL = [] # stores start bin of pulses
                    peakxR = [] # stores end bin of pulses
                    if nPulses>0:
                        realN = 1
                        peakxL.append(DigitizerPulseStartBin[ch][0])
                        peakxR.append(DigitizerPulseEndBin[ch][0])
                        for oldN in range(1,nPulses,1):
                            if ( DigitizerPulseStartBin[ch][oldN] - peakxR[realN-1] ) <= 25:
                                peakxR[realN-1] = DigitizerPulseEndBin[ch][oldN]
                            else:
                                realN = realN + 1
                                peakxL.append(DigitizerPulseStartBin[ch][oldN])
                                peakxR.append(DigitizerPulseEndBin[ch][oldN])
                        nPulses = realN
                    #remove the bad pulses whose width is <=8 ns.
                    #print evtkey, ch, nPulses
                    #for nn in range(nPulses):
                    #    print "        ", peakxL[nn], peakxR[nn]
                    nremove = 0
                    for cnt in range(0,nPulses,1):
                        if (peakxR[cnt-nremove]-peakxL[cnt-nremove]) <= 8:
                            peakxR.remove(peakxR[cnt-nremove])
                            peakxL.remove(peakxL[cnt-nremove])
                            nPulses = nPulses-1
                            nremove = nremove + 1
                    #this is specific: remove the pulses with overshoot
                    #nremove = 0
                    #for cnt in range(0,nPulses,1):
                    #    #localmaximum = thisDigi1[ch][int(peakxR[cnt-nremove])+np.argmax(thisDigi1[ch][int(peakxR[cnt-nremove]):int(peakxR[cnt-nremove])+20])]
                    #    localmaximum = max(thisDigi1[ch][int(peakxR[cnt-nremove])-4:(int(peakxR[cnt-nremove])+16)])
                    #    if localmaximum >= digi1_ped_mean[ch] + 3.0*sigma :
                    #        peakxR.remove(peakxR[cnt-nremove])
                    #        peakxL.remove(peakxL[cnt-nremove])
                    #        nPulses = nPulses-1
                    #        nremove = nremove+1
                    #finally peakxL and peakxR can tell how many pulses and in which bins
                    #print "evt " + str(int(evtkey)) + ", ch " + str(ch) + ", pulses " + str(nPulses)
                    #for ncnt in range(0,nPulses,1):
                    #    print "start at ", peakxL[ncnt], ", end at ", peakxR[ncnt]
                    #this following condition should NOT happen with above strategy, otherwise there is a problem
                    if nPulses>20:
                        print "evt ", int(evtkey), "ch ", ch, "nPulses ", nPulses
                        nPulses=20
                    #Now compute charge integration
                    for nn in range(nPulses):
                        DigitizerPulseCharge[ch][nn] = np.sum(thisDigi1[ch][int(peakxL[nn]):int(peakxR[nn])+1])-digi1_ped_mean[ch]*(peakxR[nn]-peakxL[nn]+1)
                        DigitizerPulseStartBin[ch][nn] = peakxL[nn]
                        DigitizerPulseEndBin[ch][nn] = peakxR[nn]+1
                        DigitizerPulseAmplitudeTimeBin[ch][nn] = np.argmin(thisDigi1[ch][int(DigitizerPulseStartBin[ch][nn]):int(DigitizerPulseEndBin[ch][nn])]) + int(DigitizerPulseStartBin[ch][nn])
                        DigitizerPulseAmplitude[ch][nn] = thisDigi1[ch][DigitizerPulseAmplitudeTimeBin[ch][nn]]
                    DigitizerNumberOfPulses[ch] = nPulses
                    #save waveforms if needed
                    if SaveWaveforms == 1:
                        waveformname = waveformname_head + str(int(evtkey)) + "_trigFlag_" + str(TrigFlag_Multi) + str(TrigFlag_Hodo) + str(TrigFlag_Led) + "_ch_" + str(ch) + "_nPulses_" + str(nPulses)
                        waveformtitle = waveformname + "_mean_" + str(digi1_ped_mean[ch]) + "_sigma_" + str(sigma)
                        for nn in range(nPulses):
                            waveformtitle = waveformtitle + "_" + str(peakxL[nn]) + "_" + str(peakxR[nn])
                        hist = TH1F(waveformname,waveformtitle,2560,-0.5,2560-0.5)
                        for bin in range(0,2560,1):
                            hist.SetBinContent(bin+1,thisDigi1[ch][bin])
                        waveformsDir.cd()
                        hist.Write()
                        del hist
                #for digitizer 2
                for ch in range(0,4,1):
                    #begin my own algorithm for pulse finding
                    #maximum = thisDigi2[ch][np.argmin(thisDigi2[ch][80:])+80]
                    thiswaveformcopy = [0 if x>digi2_threshold[ch] else x for x in thisDigi2[ch]]
                    DigitizerAmplitudeTimeBin[ch+4] = np.argmin(thisDigi2[ch][80:])
                    maximum = thisDigi2[ch][int(80+DigitizerAmplitudeTimeBin[ch+4])]
                    DigitizerAmplitude[ch+4] = -1.0*(maximum-digi2_ped_mean[ch])
                    #here do avearage waveform
                    if SaveWaveforms == 0 or SaveWaveforms == 1:
                        if TrigFlag_Led == 1 and TrigFlag_Hodo == 0 and TrigFlag_Multi ==0: # and totalwave_cnts[(ch+4)*3+0] <= 50:
                            totalwave_cnts[(ch+4)*3+0] += 1
                            for bin in range(0,2560,1):
                                totalwave[(ch+4)*3+0][bin] += thisDigi2[ch][bin]
                        if TrigFlag_Led == 0 and TrigFlag_Hodo == 1 and TrigFlag_Multi ==0: # and totalwave_cnts[(ch+4)*3+1] <= 50:
                            totalwave_cnts[(ch+4)*3+1] += 1
                            for bin in range(0,2560,1):
                                totalwave[(ch+4)*3+1][bin] += thisDigi2[ch][bin]
                        if TrigFlag_Led == 0 and TrigFlag_Hodo == 0 and TrigFlag_Multi ==1: # and totalwave_cnts[(ch+4)*3+2] <= 50:
                            totalwave_cnts[(ch+4)*3+2] += 1
                            for bin in range(0,2560,1):
                                totalwave[(ch+4)*3+2][bin] += thisDigi2[ch][bin]
                    nPulses = 0
                    flagStartFound = 0
                    #flagEndFound = 0
                    if DigitizerAmplitude[ch+4] <= 10: # if maximum >= digi2_threshold[ch]:
                        DigitizerNumberOfPulses[ch+4] = nPulses
                        #save waveforms if needed
                        if SaveWaveforms == 1:
                            waveformname = waveformname_head + str(int(evtkey)) + "_trigFlag_" + str(TrigFlag_Multi) + str(TrigFlag_Hodo) + str(TrigFlag_Led) + "_ch_" + str(ch+4) + "_nPulses_" + str(nPulses)
                            waveformtitle = waveformname + "_mean_" + str(digi1_ped_mean[ch]) + "_sigma_" + str(sigma)
                            hist = TH1F(waveformname,waveformtitle,2560,-0.5,2560-0.5)
                            for bin in range(0,2560,1):
                                hist.SetBinContent(bin+1,thisDigi2[ch][bin])
                            waveformsDir.cd()
                            hist.Write()
                            del hist
                        continue
                    bin=80
                    while bin>=80 and bin<(2560-4):
                        #find the start bin of a pulse in a waveform
                        #To define a cross over threshold: two bins below threshold and two bins above threshold
                        if nPulses==20: #maximum number of pulses 
                            break
                        if flagStartFound==0:
                            if thiswaveformcopy[bin]<0:
                                #if the condition is true, then it crosses threshold between bin+1 and bin+2
                                DigitizerPulseStartBin[ch+4][nPulses] = bin - 3
                                flagStartFound = 1
                                bin+=3
                                continue
                                #print "        finds a start", bin-1
                        if flagStartFound==1:
                            if thiswaveformcopy[bin] == 0:
                                DigitizerPulseEndBin[ch+4][nPulses] = bin+3
                                nPulses = nPulses+1
                                flagStartFound = 0
                                #print "            and it ends", bin+3
                        bin+=1
                    peakxL = []
                    peakxR = []
                    if nPulses>0:
                        realN = 1
                        peakxL.append(DigitizerPulseStartBin[ch+4][0])
                        peakxR.append(DigitizerPulseEndBin[ch+4][0])
                        for oldN in range(1,nPulses,1):
                            if ( DigitizerPulseStartBin[ch+4][oldN] - peakxR[realN-1] ) <= 25:
                                peakxR[realN-1] = DigitizerPulseEndBin[ch+4][oldN]
                            else:
                                realN = realN + 1
                                peakxL.append(DigitizerPulseStartBin[ch+4][oldN])
                                peakxR.append(DigitizerPulseEndBin[ch+4][oldN])
                        nPulses = realN
                    #remove the bad peaks in the spectrum whose width is <=8 ns.
                    #print evtkey, ch+4, nPulses
                    #for nn in range(nPulses):
                    #    print "        ", peakxL[nn], peakxR[nn]
                    nremove = 0
                    for cnt in range(0,nPulses,1):
                        if (peakxR[cnt-nremove]-peakxL[cnt-nremove]) <= 8:
                            peakxR.remove(peakxR[cnt-nremove])
                            peakxL.remove(peakxL[cnt-nremove])
                            nPulses = nPulses-1
                            nremove = nremove + 1
                    #this is specific: remove the pulses with overshoot
                    #nremove = 0
                    #for cnt in range(0,nPulses,1):
                    #    #localmaximum = thisDigi2[ch][int(peakxR[cnt-nremove])+np.argmax(thisDigi2[ch][int(peakxR[cnt-nremove]):int(peakxR[cnt-nremove])+20])]
                    #    localmaximum = max(thisDigi2[ch][int(peakxR[cnt-nremove])-4:(int(peakxR[cnt-nremove])+16)])
                    #    if localmaximum >= digi2_ped_mean[ch] + 3.0*sigma2 :
                    #        peakxR.remove(peakxR[cnt-nremove])
                    #        peakxL.remove(peakxL[cnt-nremove])
                    #        nPulses = nPulses-1
                    #        nremove = nremove+1
                    #print "evt " + str(int(evtkey)) + ", ch " + str(ch+4) + ", pulses " + str(nPulses)
                    #for ncnt in range(nPulses):
                    #    print "start at ", peakxL[ncnt], ", end at ", peakxR[ncnt]
                    if nPulses>20:
                        print "evt ", int(evtkey), "ch ", ch+4, "nPulses ", nPulses
                        nPulses = 20
                    #Now compute charge integration
                    for nn in range(nPulses):
                        DigitizerPulseCharge[ch+4][nn] = np.sum(thisDigi2[ch][int(peakxL[nn]):int(peakxR[nn])+1])-digi2_ped_mean[ch]*(peakxR[nn]-peakxL[nn]+1)
                        DigitizerPulseStartBin[ch+4][nn] = peakxL[nn]
                        DigitizerPulseEndBin[ch+4][nn] = peakxR[nn]+1
                        DigitizerPulseAmplitudeTimeBin[ch+4][nn] = np.argmin(thisDigi2[ch][int(DigitizerPulseStartBin[ch+4][nn]):int(DigitizerPulseEndBin[ch+4][nn])])+int(DigitizerPulseStartBin[ch+4][nn])
                        DigitizerPulseAmplitude[ch+4][nn] = thisDigi2[ch][DigitizerPulseAmplitudeTimeBin[ch+4][nn]]
                    DigitizerNumberOfPulses[ch+4] = nPulses
                    #save waveforms if needed
                    if SaveWaveforms == 1:
                        waveformname = waveformname_head + str(int(evtkey)) + "_trigFlag_" + str(TrigFlag_Multi) + str(TrigFlag_Hodo) + str(TrigFlag_Led) + "_ch_" + str(ch+4) + "_nPulses_" + str(nPulses)
                        waveformtitle = waveformname + "_mean_" + str(digi2_ped_mean[ch]) + "_sigma_" + str(sigma2)
                        for nn in range(nPulses):
                            waveformtitle = waveformtitle + "_" + str(peakxL[nn]) + "_" + str(peakxR[nn])
                        hist = TH1F(waveformname,waveformtitle,2560,-0.5,2560-0.5)
                        for bin in range(0,2560,1):
                            hist.SetBinContent(bin+1,thisDigi2[ch][bin])
                        waveformsDir.cd()
                        hist.Write()
                        del hist
            #done
            mytree.Fill()
            #if int(evtkey) >100:    break
    # average waveforms
    if SaveWaveforms == 0 or SaveWaveforms == 1:
        waveformsAvgDir.cd()
        for channel in range(0,8):
            for trigFlag in range(0,3):
                #print "AvgWaveforms: ", channelname[channel], trigFlagname[trigFlag], totalwave_cnts[channel*3+trigFlag]
                name = hWaveforms_avg[channel*3+trigFlag].GetName()
                name = name + "_NbOfTriggers_" + str(totalwave_cnts[channel*3+trigFlag])
                hWaveforms_avg[channel*3+trigFlag].SetTitle(name)
                if totalwave_cnts[channel*3+trigFlag] == 0:
                    print "caution 0 trigger: ch ", channelname[channel], "trigType ", trigFlagname[trigFlag]
                    continue
                for bin in range(0,2560):
                    hWaveforms_avg[channel*3+trigFlag].SetBinContent(bin+1,totalwave[channel*3+trigFlag][bin]/totalwave_cnts[channel*3+trigFlag])
                hWaveforms_avg[channel*3+trigFlag].Write()
    rootfile.cd()
    rootfile.Write()
    rootfile.Close()
    file.close()
    os.remove(filename)

elaspsed = timeit.default_timer() - start_time

print "Time: ", elaspsed, " seconds."

#draw histograms of QDC channels



