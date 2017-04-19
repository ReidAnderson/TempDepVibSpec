## Author: Reid Anderson <reid@reid-HP-Pavilion-dv6-Notebook-PC>
## Created: 2017-04-17

function [retval] = SaveSpectraByMode (harmFrequencies,anharmMatrix,VibInt,Emin,Emax,vmin,vmax,v_gr,binSize,occVecs,E_freq,resultsDir, runName)
numSteps = length(occVecs);
all_wn = vmin+v_gr:v_gr:vmax;

% Broadening function used in generating the spectrum
lor = @(x,x0,s) (1/pi)*(0.5*s./((x-x0).^2+(0.5*s).^2));

% Speed of light in centimeters
c_cm = 2.9979 *10^10;
output = zeros(1,length(all_wn));
map = 0;

tic
for k = 1:length(harmFrequencies)
  I = zeros(ceil((Emax-Emin)/binSize),(vmax-vmin)/v_gr);
  for steps=1:numSteps
      occVec = double(occVecs(steps,:));
      E = ceil((floor(getEnergy(harmFrequencies,anharmMatrix,occVec))-Emin)/binSize);
      
      del_E = harmFrequencies(k) + 2*anharmMatrix(k,k) + 2*anharmMatrix(k,k)*occVec(k);
      for i = 1:length(harmFrequencies)
          % Assumes the matrix is in reduced form and only has values in the
          % lower left
          if (k < i)
              del_E = del_E + anharmMatrix(i,k)+2*anharmMatrix(i,k)*occVec(i);
          end
      end
      sigma_IR = VibInt(k)*(occVec(k)+1);
      
      Evk_contrib = sigma_IR.*lor(c_cm.*all_wn,wn_to_hz(del_E),3.3E11);
  
      I(E,:) = I(E,:) + Evk_contrib;
  
      % Update on how many steps have been performed
      if mod(steps,10000) == 0
         disp([num2str(k) ':' num2str(steps)])
      end
  end
  normalizedI = zeros(size(I,1),size(I,2));
  for i = 1:length(E_freq)
      if E_freq(i) ~= 0
          normalizedI(i,:) = I(i,:)./E_freq(i);
      end
  end
  save([resultsDir '/EnergyDepVibSpec/' runName '-vibMode-' num2str(k) '-R_E'],'normalizedI');
end
freq=E_freq;
toc
retval = 1;
endfunction
