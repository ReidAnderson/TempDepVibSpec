function [nm] = eV_to_nm(eV)
%EV_TO_NM Convert a eV value or an array of eV values to nanometers
c = 2.998*10^8;
h = 6.626*10^-34;
J = 1.602*10^-19;
nm = (h*c/(eV.*J)).*10^9;
end

