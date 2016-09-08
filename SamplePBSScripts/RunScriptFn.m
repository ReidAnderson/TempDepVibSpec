function [output] = RunScriptFn(molAbbrev,runName,IRSteps,DOSIters,DOSSteps,DOSCores,raman)
% If the temperatures are not presented, then we assume that we want what
% is (realistically) experimentally possible as well as a couple higher
% temperatures for benchmarks.
temps = [250:350];
temps = [temps 500 625 750 875 1000 1125 1250 1375 1500];

IRSteps = str2num(IRSteps);
DOSIters = str2num(DOSIters);
DOSSteps = str2num(DOSSteps);
raman = str2num(raman);
DOSCores = str2num(DOSCores);

% We assume that we want to generate the DOS if we're calling from a PBS
% script
addpath('../Algorithm');
WLRun(molAbbrev,temps,0,runName,IRSteps,DOSIters,DOSSteps,raman,DOSCores);
output = 1;
exit;
end

