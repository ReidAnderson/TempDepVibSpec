function [FinalOmega] = WLPar(harmFreq, anharmMat, Emin, Emax, binSize, maxVec,iters,steps,numPar,molAbbrev)
allDOS = zeros(ceil(Emax-Emin)/binSize,numPar);
parfor i=1:numPar
    allDOS(:,i) = WL(harmFreq,anharmMat,Emin,Emax, binSize, maxVec,iters,steps);
end
% Compute the variance of the sample here
save([molAbbrev '-allDOS'],'allDOS');
FinalOmega = mean(allDOS,2);
if numPar > 10
    runName = [binSize 'bin-numPar' numPar];
    for i = 1:numPar
        DOS = allDOS(:,i);
        zero_indices = DOS == 0;
        DOS(zero_indices) = Emax+Emax;
        DOS = DOS - min(DOS);
        DOS(zero_indices) = 0;
        allDOS(:,i) = DOS;
    end
    uncertainty = zeros(size(allDOS,1),1);
    for i = 1:length(uncertainty)
        uncertainty(i,:) = std(allDOS(i,:));
    end
    uncertainty = uncertainty./sqrt(numPar);
    % For now, don't save uncertainty
    %save(['Results/stdev_' molAbbrev '-' runName 'iter' num2str(iters)],'uncertainty');
end
end
