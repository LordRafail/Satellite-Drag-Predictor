These are the files that calculate the drag coefficient using historical and forecasted space weather data. The drag ceofficient is calculated using GSI models and environment functions from ADBSAT. 
The user needs to download and install ADBSAT in the same directory as the files in this repository for everything to work. An active internet connection is also required. 
ADBSAT Github repository download link: https://github.com/nhcrisp/ADBSat
MATLAB addons needed to run the files are: Aerospace Blockset, Aerospace Toolbox, Econometrics toolbox, optimization toolbox, Simulink, statistics and machine learning toolbox and UAV toolbox.  

The rest of the files are my own work. 
You also need to go on to https://swe_bgs_ac_uk.content.swe.s2p.esa.int/indices.shtml?index=apforecast72hr, create an account with ESA, 
and download the .Json file with the 72hr AP index forecast and put it in the working directory. If you don't do this predicting the drag for a future date will not work. 
I hope everything works fine, if not just send me an email at rafail.panagiotidis@postgrad.manchester.ac.uk

Rafail Panagiotidis 2025.
