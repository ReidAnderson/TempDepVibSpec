T = [250, 275, 300, 325, 350];
% T = [250];
T = [T 500 625 750 875 1000 1125 1250 1375 1500];
%T = [3300];
molAbbrev = 'CHCl3';
% Have we already computed the DOS for this run?
DOSExists = 1;
runName = '3-30-2017-Raman';
% How many steps do we want to take in the Monte Carlo that computes the spectra?
IRSteps = 50000;
% How many iterations of the WL do we want to perform, and what's the maximum number of steps in each?
DOSIters = 20;
DOSSteps = 2000000;
% Should we compute the Raman spectra instead of the IR spectra?
raman = 1;
% How many independent runs do we want to average for the DOS?
numCores = 5;
WLRun(molAbbrev,T,DOSExists,runName,IRSteps,DOSIters,DOSSteps,raman,numCores);
molAbbrev = '';
