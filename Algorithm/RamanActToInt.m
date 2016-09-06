function ramanInt = RamanActToInt(harmFrequencies, ramanAct, v0, T)
B = zeros(1,length(harmFrequencies));
% v0 = 9398.5; 
ramanInt = zeros(1,length(harmFrequencies));
% 
% This component is temperature dependent which seems like a problem. So we
% need to perform the computations using the scattering values and convert
% at the end.

kb = 1.3806485*10^-23;
% Planck's constant in Js
h = 6.626*10^-34;
c = 2.998*10^8;
C = 10^-12;
% Assume B(i)=1 for now, and then add the effect back in at the end
for i = 1:length(B)
%     B(i) = 1-exp(-(h*wn_to_hz(harmFrequencies(i)*c))/(kb*T));
    B(i) = 1;
end

for i = 1:length(ramanAct)
    ramanInt(i) = C*(v0-harmFrequencies(i))^4*(harmFrequencies(i))^-1*B(i)^-1*ramanAct(i);
end
end