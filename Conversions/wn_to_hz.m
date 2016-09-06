function [hz] = wn_to_hz(wn)
%WN_TO_HZ Convert a value in wavenumbers to a value in Hertz
wn_m = wn*100;
c = 2.998*10^8;
wavelength = 1./wn_m;
hz = c./wavelength;
end
