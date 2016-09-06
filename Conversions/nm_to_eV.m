function [eV] = nm_to_ev(nm)
%NM_TO_EV Convert a nanometer value or array of nanometer values to eV
c = 2.998*10^8;
h = 6.626*10^-34;
J = 1.602*10^-19;
eV = (h*c./(nm.*10^-9))./J;
end

