function [output map] = GetAbsorption(harmFrequencies, anharmMatrix, IRInt, occVec, all_wn, E)
% Broadening function used in generating the spectrum
lor = @(x,x0,s) (1/pi)*(0.5*s./((x-x0).^2+(0.5*s).^2));

% Speed of light in centimeters
c_cm = 2.9979 *10^10;
output = zeros(1,length(all_wn));
map = 0;

vibWeight = squeeze(spectra.*(double(spectraWeights(E,:,:))./65535));

for k=1:length(harmFrequencies)
    del_E = harmFrequencies(k) + 2*anharmMatrix(k,k) + 2*anharmMatrix(k,k)*occVec(k);
    for i = 1:length(harmFrequencies)
        % Assumes the matrix is in reduced form and only has values in the
        % lower left
        if (k < i)
            del_E = del_E + anharmMatrix(i,k)+2*anharmMatrix(i,k)*occVec(i);
        end
    end
    sigma_IR = IRInt(k)*(occVec(k)+1);
    
    Evk_contrib = sigma_IR.*lor(c_cm.*all_wn,wn_to_hz(del_E),3.3E11);
    output = output + Evk_contrib;
end

test = 2+2;
output = sum(vibWeight,2);
EspectraWeights = uint16((vibWeight./output).*65535);
test = 1+1;
output = output';
end