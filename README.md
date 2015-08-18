# pythiaNPEh
Code for Pythia NPE-h Template Construction


This code is run on RCF using condor_submit. The .job file submits all of the jobs individually, acting as a distributor of jobs over the various nodes. It grabs the script file, which actually runs the job by grabbing the card of interest. The card has all the Pythia specific inputs, such as available processes, beam type and energy, etc. The script file also determines where the output file goes and the name of the histograms. 

The .cpp file determines how to handle each generated event after Pythia generates it. This looks at tracks in the event and decides which ones are of interest, then histograms the output. 
