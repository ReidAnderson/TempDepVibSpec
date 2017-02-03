function [FinalOmega allDOS] = WL(harmFreq, anharmMat, Emin, Emax, binSize, maxVec,numIters,numSteps)
tic
steps = 1;
iter = 1;
cont = 1;
zero_threshold = 0.3;
e = exp(1);

% n = GetRandomVector(harmFreq,anharmMat,Emin,Emax,binSize,3000,length(harmFreq));
n = zeros(1,length(harmFreq),1);
g = zeros(1,ceil(Emax-Emin)/binSize);
H = zeros(1,ceil(Emax-Emin)/binSize);
f = e;
p = 0.1;
allDOS = zeros(numIters,length(g));

while(cont > 0)
    if (cont == 2)
        allDOS(iter,:) = g;
        iter
        steps
        cont = 1;
        steps = 1;
        f = sqrt(f);
        H = zeros(1,ceil(Emax-Emin)/binSize);
        n = zeros(1,length(harmFreq));
    end

    n_old = n;
    rnums = rand(length(n));
    for i = 1:length(n)
       delta = 0;
       if rnums(i) < p
           delta = -1;
       elseif rnums(i) > 1-p && n(i)<maxVec(i)
           delta = 1;
       end
       n(i) = n(i) + delta;
       if n(i) < 0
           n(i) = 0;
       end
    end
    E_old = floor(getEnergy(harmFreq,anharmMat,n_old))-Emin;
    E_new = floor(getEnergy(harmFreq,anharmMat,n))-Emin;

    E_old_bin = ceil((floor(getEnergy(harmFreq,anharmMat,n_old))-Emin)/binSize);
    E_new_bin = ceil((floor(getEnergy(harmFreq,anharmMat,n))-Emin)/binSize);
    acc_prob = 1;

    if E_new <= 0 || E_new > Emax-Emin-binSize
        acc_prob = 0;
    elseif E_old_bin > 0 && E_new_bin > 0
        if g(E_old_bin) < g(E_new_bin) 
            acc_prob = e^(g(E_old_bin)-g(E_new_bin));
        else
            acc_prob = 1;
        end
    end
    acc_rand = rand();
    
    if acc_rand < acc_prob
        H(E_new_bin) = H(E_new_bin)+1;
        g(E_new_bin) = g(E_new_bin)+f;
    else
        H(E_old_bin) = H(E_old_bin)+1;
        g(E_old_bin) = g(E_old_bin)+f;
        n=n_old;
    end
    
    steps = steps+1;
    
    dE_all = 1;
    avgH = mean(H);
    num0 = 0;
    
    % We can expect some bins to be zero so right now we say that if less
    % than 5 percent of them are we have a sufficiently flat histogram
    for i = 1:length(H)
       if abs((H(i)-avgH)/avgH) > 0.9 && H(i) ~= 0
           dE_all = 0;
       elseif H(i) == 0
           num0 = num0 +1;
       end
    end
    
    if num0 > zero_threshold*length(H)
        dE_all = 0;  
    end

    % Update as the computation runs
    if mod(steps,100000) == 0
        disp(steps);
    end
    
    if steps > numSteps || dE_all == 1 
        iter = iter +1;
        if iter > numIters
            cont = -1;
        else
            cont = 2;
        end
    end
end
FinalOmega = g;
toc
end
