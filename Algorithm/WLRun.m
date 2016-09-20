function [I_T I freq] = WLRun(molAbbrev, T, DOSexists, runName, IRsteps,DOSiters,DOSsteps,raman,DOScores)

addpath('../Conversions/');
%freqRun = '-zero';
freqRun = '';
freqDir = ['../Results/' molAbbrev '/'];
resultsDir = ['../Results/' molAbbrev '/' runName '/'];
load([freqDir 'Freqs/' molAbbrev freqRun '-harmFreq']);
load([freqDir 'Freqs/' molAbbrev freqRun '-anharmMatrix']);
load([freqDir 'Freqs/' molAbbrev freqRun '-IRInt']);

if raman == 1
    load([freqDir 'Freqs/' molAbbrev freqRun '-RamanAct.mat']);
    RamanInt = RamanActToInt(harmFrequencies,RamanAct,9398.5,300);
    IRInt = RamanAct;
%     IRInt = RamanInt;
end

% TODO: Just setting the raman activities to equal IRInt is lazy, but it
% makes things really easy. This should change at some point.
zeropoint = floor(getEnergy(harmFrequencies,anharmMatrix,zeros(length(harmFrequencies),1)))-1;
binSize = 40;
upperEnergy = 47200;
energies = zeropoint:binSize:(zeropoint+upperEnergy);
energies_J = wn_to_J(energies);
kb = 1.3806485*10^-23;

% We set upper limits on the maximum values of the occupation vector so
% that the anharmonic energy is always in the increasing regime
maxVec = getMaxOccVec(harmFrequencies,anharmMatrix);

if DOSexists == 0
    % Calculate the density of states
    % second to last argument is number of independent runs to average over
    [DOS] = WLPar(harmFrequencies,anharmMatrix,zeropoint,zeropoint+upperEnergy, binSize, maxVec,DOSiters,DOSsteps,DOScores,molAbbrev);
    
    % Scale the resulting DOS so that it starts at 0
    zero_indices = DOS == 0;
    DOS(zero_indices) = max(DOS)+upperEnergy;
    DOS = DOS - min(DOS);
    DOS(zero_indices) = 0;

    if ~exist([resultsDir '/DOS'],'dir')
        mkdir([resultsDir '/DOS']);
    end
    save([resultsDir '/DOS/' runName '-DOS.mat'],'DOS');
else
    load([resultsDir '/DOS/' runName '-DOS']);
end
V_gr = 0.1;
vmin = 0;
vmax = 3500;
all_wn = vmin:V_gr:vmax;
% Dummy uncertainty values for now
% DOS_err = ones(size(DOS));
[I N h] = IntensityWL(harmFrequencies,anharmMatrix,IRInt,zeropoint,zeropoint+upperEnergy,vmin,vmax,DOS,V_gr,binSize,IRsteps);
normalizedI = zeros(size(I,1),size(I,2));
for i = 1:length(N)
    if N(i) ~= 0
        normalizedI(i,:) = I(i,:)./N(i); 
    end
end

%uncertainties = BootstrapIR(allAbs);

if ~exist([resultsDir '/EnergyDepVibSpec'],'dir')
    mkdir([resultsDir '/EnergyDepVibSpec']);
end
if raman == 0
    save([resultsDir '/EnergyDepVibSpec/' runName '-I_E'],'normalizedI');
else
    save([resultsDir '/EnergyDepVibSpec/' runName '-R_E'],'normalizedI');
end

if raman==0
    % Generate I_T for each of the specified temperatures
    for idx_T = 1:length(T)
        % Now do a Laplace transform to make I(v,E) into I(v,T)
        % Z is partition function
        Z = 0;
        for i = 1:length(DOS)
            Tdep = exp(-energies_J(i)/(kb*T(idx_T)));
            Z = Z+DOS(i)*Tdep;
        end
        
        I_T = zeros(1,length(all_wn)-1);
        for i = 1:length(DOS)
            Tdep = exp(-energies_J(i)/(kb*T(idx_T)));
            next = normalizedI(i,:)*DOS(i)*Tdep;
            I_T = I_T + next;
        end
        I_T = I_T*(1/Z);
        
        if ~exist([resultsDir '/TempDepVibSpec'],'dir')
            mkdir([resultsDir '/TempDepVibSpec']);
        end
        if ~exist([resultsDir '/TempDepVibSpec/IR'],'dir')
            mkdir([resultsDir '/TempDepVibSpec/IR']);
        end
        save([resultsDir '/TempDepVibSpec/IR/' runName '-' num2str(T(idx_T))],'I_T')
    end
    
else
    % Generate I_T for each of the specified temperatures
    for idx_T = 1:length(T)
        h = 6.626*10^-34;
        c_cm = 2.998*10^10;
        B = zeros(length(harmFrequencies),1);
        % Assume B(i)=1 for now, and then add the effect back in at the end
        for i = 1:length(B)
            B(i) = 1-exp(-(h*harmFrequencies(i)*c_cm)/(kb*T(idx_T)));
            %         B(i) = 1;
        end

        % Now do a Laplace transform to make I(v,E) into I(v,T)
        % Z is partition function
        Z = 0;
        for i = 1:length(DOS)
            Tdep = exp(-energies_J(i)/(kb*T(idx_T)));
            Z = Z+DOS(i)*Tdep;
        end
        
        I_T = zeros(1,length(all_wn)-1);
        for i = 1:length(DOS)
            Tdep = exp(-energies_J(i)/(kb*T(idx_T)));
            next = normalizedI(i,:)*DOS(i)*Tdep;
            I_T = I_T + next;
        end
        I_T = I_T*(1/Z);
        
        if ~exist([resultsDir '/TempDepVibSpec'],'dir')
            mkdir([resultsDir '/TempDepVibSpec']);
        end
        if ~exist([resultsDir '/TempDepVibSpec/Raman'],'dir')
            mkdir([resultsDir '/TempDepVibSpec/Raman']);
        end
        save([resultsDir '/TempDepVibSpec/Raman/' runName '-' num2str(T(idx_T))],'I_T')
    end
end
end
