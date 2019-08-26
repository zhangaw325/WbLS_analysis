# WbLS_analysis
Script collection for analyzing the one ton detector

1. Read_MultipleFiles_WriteTree_v9.py [Analyzing the raw data files]
-> Raw data files are in .h5 format, and are zipped in .zip files
-> Raw data files include:
    waveforms from 8 digitizer channels for the 8 PMTs, 
    16 TDC channels (if no hit in the detector channel, no number will be recorded, in other words, the length of the array is <=16),
    17 scaler channels recording the accumalated number of triggers in each detector
    some other environmental variables: level, temperature, time, etc
-> the python code takes one parameter of the raw data file in zip format
-> the ProcessRawDataJob_batch.job is used to process data files in a batch mode by submitting jobs
-> the batch_rawdata.sh is needed, 
-> also a text file listing the raw data file names is necessary
