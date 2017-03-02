function totalEnergy = getEnergy(harmFrequencies, anharmMat, occVec)
% There is inconsistent about the dimensions of occVec that I pass in 
% (Nx1 vs 1xN), so this makes sure that the dimensions are consistent
% and the matrix operations work.
if (size(occVec,1) > size(occVec,2))
    occVecT = occVec';
else
    occVecT = occVec;
    occVec = occVec';
end
negAnharm = 0;
e1 = harmFrequencies.*(occVecT+0.5);
e2 = ((occVec+0.5)*(occVecT+0.5)).*anharmMat;
modeEnergy = e1+sum(e2);
for i = 1:length(modeEnergy)
    if modeEnergy(i) < 0
        negAnharm = 1;
    end
end

if negAnharm ~= 1
    totalEnergy = sum(modeEnergy);
else
    totalEnergy = -1;
end
end

