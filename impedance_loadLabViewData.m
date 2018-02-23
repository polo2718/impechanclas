function [f,R,X]=impedance_loadLabViewData(fileref)
%Load impedance data from my Labview program Dataset
%--------------------------------------------------------------------------
%Usage:
%[f,R,X]=impedance_loadLabViewData(fileref)
%--------------------------------------------------------------------------
%Inputs:
%fileref: String pointing to file reference e.g. fileref='data.txt'
%--------------------------------------------------------------------------
%Outputs:
%-f: frequency data points in [Hz], this is the independent non-random
%   variable that will be used in regression. This is as 1D column vector (1XM)
%-R  resistance data points:value of real part of impedance [ohms]. This is as 1D column vector (1XM)
%-X reactance datapoints: value of imaginary part of impedance [ohms]. This is as 1D column vector (1XM)
%--------------------------------------------------------------------------
%Author:Leopoldo Cendejas Zaragoza-2018
data=load(fileref);
data=sortrows(data,1);%sort data from smallest to largest frequency
f=data(:,1)';        %frequency [Hz] (independent, non random variable)
R=data(:,5)';        %observed resistance [ohms]
X=data(:,6)';        %observed reactance [ohms]
