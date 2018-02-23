function [mag,phase]=impedance_rect2polar(R,X)
%Convert impedance measurements from rectangular to polar coordinates
%--------------------------------------------------------------------------
%Usage:
%[mag,phase]=impedance_rect2polar(R,X)
%--------------------------------------------------------------------------
%Inputs:
%-R resistance:value of real part of impedance [ohms]
%-X reactance: value of imaginary part of impedance [ohms]
%--------------------------------------------------------------------------
%Outputs:
%mag: impedance magnitude [ohms]
%phase: impedance phase [degrees]
%Leopoldo Cendejas Zaragoza-2018
mag=sqrt(R.^2+X.^2);%compute model magnitude
phase=rad2deg(atan(X./R));%compute model phase
