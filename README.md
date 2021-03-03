# Linelength-spike-detector
Detects abnormal events in brain wave data (namely, interictal spikes in EEG data) using a linelength transform algorithm.

Transforms data into linelength then detects events (spikes) surpassing
the designated percentile threshold. Note that this function assumes any 
detections in any channel occurring simultaneously are involved in the 
same spike event. 
Based on Estellar et al 2001, DOI 10.1109/IEMBS.2001.1020545
INPUTS
  d: vector or matrix of ICEEG data and 
  sfx: sampling frequency
  llw: linelength window (in seconds) over which to calculate transform
  prc: percentile to use as a threshold for detections
  badch: logical index of bad channels (1=bad, 0=ok)
OUTPUTS
  ets: matrix of events (rows) and their on/off times (2 columns) in samples
  ech: logical index of which channels are involved in each detection, 
       thus having the same number of rows (spikes) as ets
  
Example: [ets,ech]=LLspikedetector(d,512,.04,99.99)

jon.kleen@ucsf.edu 2016-2021
