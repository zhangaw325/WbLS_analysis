# WbLS_analysis
Script collection for analyzing the one ton detector

1. Analyzing the raw data files
-> Raw data files are in .h5 format, and are zipped in .zip files
-> Raw data files include:
    waveforms from 8 digitizer channels for the 8 PMTs, 
    16 TDC channels (if no hit in the detector channel, no number will be recorded, in other words, the length of the array is <=16),
    17 scaler channels recording the accumalated number of triggers in each detector
    some other environmental variables: level, temperature, time, etc
