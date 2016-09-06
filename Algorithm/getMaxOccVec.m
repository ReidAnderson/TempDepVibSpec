function maxVecs = getMaxOccVec(harmFreq,anharmMatrix)
tic
    maxVecs = zeros(length(harmFreq),1);
    for i = 1:length(harmFreq)
        E = 0;
        E_old = 0;
        n = 0;
        cont =1;
        while (cont ==1)
            E = harmFreq(i)*(n+0.5)+anharmMatrix(i,i)*(n+0.5)^2;
            delE = E-E_old;
            if delE < 0 || n > 20000
                maxVecs(i) = n;
                cont = 0;
            end
            n = n +1;
            E_old = E;
        end
    end
    toc
end