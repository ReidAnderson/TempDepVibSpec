function [FinalOmega] = WLPar(harmFreq, anharmMat, Emin, Emax, binSize, maxVec,iters,steps,numPar,molAbbrev)
allDOS = zeros(ceil(Emax-Emin)/binSize,numPar);
parfor i=1:numPar
    allDOS(:,i) = WL(harmFreq,anharmMat,Emin,Emax, binSize, maxVec,iters,steps);
end
%save([molAbbrev '-allDOS'],'allDOS');
FinalOmega = mean(allDOS,2);
end
