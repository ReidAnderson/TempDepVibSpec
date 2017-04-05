function [I freq spectra] = IntensityWL(harmFreq, anharmMat, IRInt, Emin, Emax, vmin, vmax, DOS,v_gr,binSize,steps)
numSteps = steps;
all_wn = vmin+v_gr:v_gr:vmax;
%n = GetRandomVector(harmFreq,anharmMat,Emin,Emax,binSize,3000,length(harmFreq));
n = zeros(1,length(harmFreq),1);
g = DOS;
e = exp(1);
p = 0.1;
I = zeros(ceil((Emax-Emin)/binSize),(vmax-vmin)/v_gr);
harmModesUsed = [];
E_freq = zeros(ceil((Emax-Emin)/binSize),1);

for i = 1:length(IRInt)
  if IRInt(i) > 0
    harmModesUsed = [harmModesUsed i];
  end
end

spectra = zeros(ceil((Emax-Emin)/binSize),(vmax-vmin)/v_gr,length(harmModesUsed));

tic
steps = 1;
while(steps < numSteps)
    n_old = n;
    rnums = rand(length(n));
    for i = 1:length(n)
       delta = 0;
       if rnums(i) < p
           delta = -1;
       elseif rnums(i) > 1-p
           delta = 1;
       end
       n(i) = n(i) + delta;
       if n(i) < 0
           n(i) = 0;
       end
    end

    E_old_bin = ceil((floor(getEnergy(harmFreq,anharmMat,n_old))-Emin)/binSize);
    E_new_bin = ceil((floor(getEnergy(harmFreq,anharmMat,n))-Emin)/binSize);
    acc_prob = 1;

    % We automatically reject if the bin is above or below our energy
    % range
    if E_new_bin <= 0 || E_new_bin > length(g) || E_old_bin > length(g)
        acc_prob = 0;
    elseif E_old_bin > 0 && E_new_bin > 0
        if g(E_old_bin) < g(E_new_bin) 
            acc_prob = e^(g(E_old_bin)-g(E_new_bin));
        else
            acc_prob = 1;
        end
    end
    
    acc_rand = rand();
    AbsVal = 0;
    if acc_rand < acc_prob
        [AbsVal modesMap] = GetAbsorption(harmFreq,anharmMat,IRInt,n,all_wn,E_new_bin, spectra, harmModesUsed);
        I(E_new_bin,:) = I(E_new_bin,:) + AbsVal;
        E_freq(E_new_bin) = E_freq(E_new_bin)+1;
        allAbs(steps,1) = E_new_bin;
        spectra = spectra+modesMap;
    else
        n=n_old;
        [AbsVal modesMap] = GetAbsorption(harmFreq,anharmMat,IRInt,n,all_wn,E_old_bin, spectra, harmModesUsed);
        I(E_old_bin,:) = I(E_old_bin,:) + AbsVal;
        E_freq(E_old_bin) = E_freq(E_old_bin)+1;
        spectra = spectra+modesMap;
    end
    
    steps = steps+1;

    % Update on how many steps have been performed
    if mod(steps,10000) == 0
       steps
    end
end
freq=E_freq;
toc
end
