molAbbrev = 'CO2';
%WLRun(molAbbrev,T,0,'10k-sigma',10000,10,2000000);

DOSiters=20;
DOSsteps=2000000;
addpath('../../Macros/Reid/');
%freqRun = '-zero';
freqRun = '';
load([molAbbrev freqRun '-harmFreq']);
load([molAbbrev freqRun '-anharmMatrix']);
load([molAbbrev freqRun '-IRInt']);
zeropoint = floor(getEnergy(harmFrequencies,anharmMatrix,zeros(length(harmFrequencies),1)))-1;
binSize = 40;
upperEnergy = 47200;
energies = zeropoint:binSize:(zeropoint+upperEnergy);
energies_J = wn_to_J(energies);
kb = 1.3806485*10^-23;

% We set upper limits on the maximum values of the occupation vector so
% that the anharmonic energy is always in the increasing regime
maxVec = getMaxOccVec(harmFrequencies,anharmMatrix);

% Average a bunch of runs together to get the variance
% for i=1:DOSiters
%     [DOS] = WLPar(harmFrequencies,anharmMatrix,zeropoint,zeropoint+upperEnergy, binSize, maxVec,i,DOSsteps,30,molAbbrev);
% end

% [DOS] = WLPar(harmFrequencies,anharmMatrix,zeropoint,zeropoint+upperEnergy, binSize, maxVec,20,DOSsteps,1,molAbbrev);