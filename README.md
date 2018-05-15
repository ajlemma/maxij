# maxi j1820+070 image analysis code
for ~1s cadence optical data images taken with RHO 14" optical telescope


SSE, AJT, SKJ 


* maxijdefs.py: module with functions referenced by other programs (need this to run other programs)
* maxij_imshifts.py: calculates shifts with cross-correlation & multiprocessing
* maxij_imalign.py: shifts images & rebins them, saves aligned images
