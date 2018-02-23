function impedance_regression_plot(f,R,X, R_model, X_model)
%--------------------------------------------------------------------------
%Inputs:
%-f: 1D column vector with M data points specifying frequency data [Hz].
%   This is the independent non-random variable used for regression.
%-R: column vector containing the observed resitance values [ohms] (This is the real part of the impedance).
%   Note that the number of rows should be equal to the number of rows in the 
%   frequency input vector f (i.e MagZ is a column vector of Mx1).        
%-X: observed reactance data points [ohms] (imaginary part of impedance).
%   Note that X is a column vector with dimensions (Mx1)
%-R_model: estimated resistance using regression model[ohms]
%-X_model: estimated reactance using regression model[degrees]
%
%--------------------------------------------------------------------------
%Leopoldo Cendejas Zaragoza 2018

figure('DefaultAxesFontSize',14);
subplot(1,3,1);
%Graph magnitude
[MagZ,PhaseZ]=impedance_rect2polar(R,X);
[MagZ_model,PhaseZ_model]=impedance_rect2polar(R_model,X_model);
semilogx(f, MagZ,'*',f,MagZ_model, 'LineWidth', 1);
title('Impedance magnitude vs frequency');
xlabel('Frequency [Hz]');
ylabel('Magnitude [Ohms]');
%xlim(flimits); ylim(maglimits);
legend('Data Points', 'Regression Model');

subplot(1,3,2);
semilogx(f, PhaseZ,'*',f,PhaseZ_model,'LineWidth', 1);
title('Impedance phase vs frequency');
xlabel('Frequency [Hz]');
ylabel('Phase [Degrees]');
%xlim(flimits); ylim(phaselimits);
legend('Data Points', 'Regression Model');

subplot(1,3,3);
plot(R,X,'*',R_model, X_model, 'LineWidth', 1);
title('Nyquist plot');
xlabel('Resistance [Ohms]');
ylabel('Reactance [Ohms]');
%xlim(resistancelimits); ylim(reactancelimits);
legend('Data Points', 'Regression Model');