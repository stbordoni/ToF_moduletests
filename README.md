# ToF_moduletests
preliminary data analysis for ND280upgrade ToF modules


Instruction from Alexander: 


I send you 3 files:

  MakeRootTree.C  -  convert the *.bin file to *.root

  DrawTrackOnly.C -  read *.root tree and plot online the position of hits (also fitted track)

  read_data_WC.h  -  decoding information

So, you should do the following. Copy all three files and also *.bin files in one directory.
Update the *.bin name at the beginning of MakeRootTree.C and run it:

  root -b -q MakeRootTree.C+

As an output you will have a root file with the same name (  .bin -> .root ). It will contain
the root tree with a decoded information, like time measured at each side of every bar.

Update the root file name at the beginning of DrawTrackOnly.C and run it

   root DrawTrackOnly.C

A canvas will popup and you will see hits and tracks. As I remember, you should have 
~500 good straight tracks per 10 kev raw data file.

The goal of the alignment procedure will be to turn the signal velocity along the bars
and find out their relative position. The steps should be the following:

1. Calculate the hit location along the bar assuming some default values of velocity and bar’s position.

2. Fit the track with a straight line and remove all bad chi2 tracks.

3. Plot residuals,  (x_track - x_hit) vs x_hit,  for every bar. If you have a nicely horizontal 
     distribution that means all is perfect (no need to for additional alignment). If you see 
     a non-zero slope in the distribution, you should adjust the velocity. If you see a constant
     shift from zero, you should adjust the bar’s position.

4. the procedure should be interactive: once you adjust the velocity and position, you should
     recalculate the hit position and refit the track. 
