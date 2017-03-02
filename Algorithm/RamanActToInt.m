function ramanInt = RamanActToInt(harmFrequencies, ramanAct, v0, T)
% v0 = 9398.5; 
ramanInt = zeros(1,length(harmFrequencies));

kb = 1.3806485*10^-23;
% Planck's constant in Js
h = 6.626*10^-34;
c = 2.998*10^8;
C = 10^-12;
% Assume B(i)=1 for now, we will add the effect back in for each temperature at the end
B = ones(1,length(harmFrequencies));

for i = 1:length(ramanAct)
    ramanInt(i) = C*(v0-harmFrequencies(i))^4*(harmFrequencies(i))^-1*B(i)^-1*ramanAct(i);
end
end