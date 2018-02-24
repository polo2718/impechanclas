%Ejemplo
fileref='cnt_training.txt';
[f,R,X]=impedance_loadLabViewData(fileref);
%% valores iniciales
p0=[400,329,1E-12,1,1E-9,0.5];
%% calibracion
Rcalibration=200;
%% dominio de solucion
lb=[0,0,1E-12,0,1E-12,0];
ub=[1E6,1E6,100E-6,1,100E-6,1];
%% ejecutar la regresion
[p,f_clean,R_clean,X_clean,R_model,X_model]=impedance_regression_Warburg(f,R, X,...
    p0, Rcalibration, lb ,ub);
%% graficar
impedance_regression_plot(f_clean, R_clean, X_clean, R_model, X_model);