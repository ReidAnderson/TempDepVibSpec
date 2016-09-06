function [J] = wn_to_J(wn)
%WN_TO_J Convert a value in wavenumbers to Joules
wn_m = wn.*100;
c = 2.998*10^8;
h = 6.626*10^-34;
wavelength = 1./wn_m;
hz = c./wavelength;
E = h.*hz;
J = E;
end
