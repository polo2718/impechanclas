function [p,f_clean,R_clean,X_clean,R_model,X_model]=impedance_regression_Warburg(f,R, X, p0, Rcalibration, lb,ub)
%Perform a non-linear regression to fit impedance data to an impedance model that includes a Warburg element and a non-ideal capacitor.
%The function returns the calculated model parameters as well as the modeled resistance and reactance.
%--------------------------------------------------------------------------
%Usage
%[p,f_clean,R_clean,X_clean,R_model,X_model]]=impedance_regression_Warburg(f,R, X, p0, Rcalibration, lb ,ub)
%--------------------------------------------------------------------------
%Requirements
%This function requires the global optimization toolbox, and the parallel
% computation toolbox to be able to minimize the mean square error using a multistart algorithm
% using different threads.
%--------------------------------------------------------------------------
%Inputs:
%-f: 1D column vector with M data points specifying frequency data [Hz].
%   This is the independent non-random variable used for regression.
%-R: column vector containing the observed resitance values [ohms] (This is the real part of the impedance).
%   Note that the number of rows should be equal to the number of rows in the 
%   frequency input vector f (i.e MagZ is a column vector of Mx1).        
%-X: observed reactance data points [ohms] (imaginary part of impedance).
%   Note that X is a column vector with dimensions (Mx1)
%-p0 (optional but important): vector containing initial approximation for the model parameters
%   -p0(1)=Rint: initial guess for Rint (resistance associated with charge transfer) [ohms]
%   -p0(2)=Rm: initial guess for Rm (resistance associated with medium) [ohms]
%   -p0(3)=C: initial guess for C [F]
%             This value is associated with the double layer capacitance.
%   -p0(4)=n: initial guess for n
%             This value is associated with the double layer capacitance.
%   -p0(5)=A: initial guess for A [F]
%             This value is associated with a warburg impedance  
%   -p0(6)=m: initial guess with a Warburg impedance
%-Rcalibration(optional but important): value of calibration resistance
%   [ohms] This resistance will be subtracted from the observed data points.
%   (default=0)
%lb (optional): (1x5) vector containing the lower bound for the five parameters (see p in the output section)
%ub (optional): (1x5) vector containing the upper bound for the five parameters (see p in the output section)
%--------------------------------------------------------------------------
%Outputs:
%-R_model: estimated resistance using regression model[ohms]
%-X_model: estimated reactance using regression model[degrees]
%-p_sol: The regression solution will be a parameter vector containing the values for the impedance model
%   -p(1)=Rint: fitted parameter Rint (resistance associated with charge transfer) [ohms]
%   -p(2)=Rm: fitted parameter for Rm (resistance associated with medium) [ohms]
%   -p(3)=C: fitted parameter for C [F]
%             This value is associated with the double layer capacitance.
%   -p(4)=n: fitted parameter for n
%             This value is associated with the double layer capacitance.
%   -p(5)=A: fitted parameter for A
%             This value is associated with a warburg impedance  
%   -p(6)=;: This value is associated with a Warburg impedance
%-f_clean:Vector containing frequency points without outliers
%-R_clean:Vector containing resistance points without outliers
%-X_clean: Vector containing reactance points without outliers
%Author:Leopoldo Cendejas-Zaragoza 2018

%Set default values according to input parameters
if ~exist('Rcalibration', 'var')||isempty(Rcalibration)
    Rcalibration=0; %set default value for Rcalibration
end

%Set default parameter constraintts
if ~exist('lb', 'var')||isempty(lb)
    lb=[0,0,1E-12,1,1E-12,0]; %set default value for lower bound
end
if ~exist('ub', 'var')||isempty(ub)
    %Set parameter constraints
    ub=[1E6,1E6,100E-6,1,100E-6,1];%upper bound; %set default value for upperbound
end
%Set default initial conditions
if ~exist('p0', 'var')||isempty(p0)
    %Set parameter constraints
    p0=[400,329,1E-12,1,1E-9,0.5];%default initial parameter guess
end
%Perform calibration if needed
if Rcalibration<0
    error('Rcalibration should be a positive number');
elseif Rcalibration~=0
    R=R-Rcalibration;
end

% Remove outliers from data
%Detect in which positions we have an outlier (an outlier is defined as a point
%that is more than 3 median absolute deviations away from the median)
outlierR=isoutlier(R,'movmedian',15,...
    'ThresholdFactor', 4);%%Detect outliers that happen in a window of 30 neighbors using the resistance data
outlierX=isoutlier(R,'movmedian',15,...
    'ThresholdFactor', 4);%Detect outliers in reactance data
outliers=outlierR|outlierX; %this variable returns a vector that has a value of 1 in the outlier positions
%Detect very high frequencies
outliers(f>=950E3)=1;
%Detect noisy range
outliers(f>=23 & f<=100)=1;
%Delete outliers from data
f(outliers==1)=[];
R(outliers==1)=[];
X(outliers==1)=[];


% Create and solve optimization problem
%The regression solution will be a parameter vector containing the values
%for the impedance model
%   -p(1)=Rint: model parameter [ohms]
%   -p(2)=Rmedium: model parameter [ohms]
%   -p(3)=C: model parameter [1/s] (this is not an actual capacitor, see my notes)
%   -p(4)=n: Value of a constant
%   -p(5)=A: Value of a constant

%Create optimization problem structure, use least square solver
options= optimoptions('lsqcurvefit', 'FunctionTolerance', 1E-16); %adjust tolerance for matlab least square fit algorithm
problem= createOptimProblem('lsqcurvefit','x0',p0,'objective',@impedance_model_Warburg,...
    'lb', lb, 'ub', ub, 'xdata', f, 'ydata', [R;X], 'options', options);

%Set a multistart problem to find global minimum
par=gcp('nocreate');
if isempty(par)
    parpool('local');%Use parallel computation
end
ms=MultiStart('PlotFcns', [], 'XTolerance', 1E-16, 'UseParallel', true);
%ms=MultiStart('PlotFcns', @gsplotbestf, 'XTolerance', 1E-16, 'UseParallel', true);
[p, errorOpt]=run(ms, problem,200); %p is a vector tha contains the best fit parameters, errorOpt contains the squared norm error at this point
disp('error:')
disp(errorOpt/size(f,1));
%delete(gcp('nocreate'))
%Compute R_model and X_model vectors 
temp=impedance_model_Warburg(p, f);%evaluate the model with the calculated parameters
R_model=temp(1,:); %compute resistance points according to model parameters
X_model=temp(2,:); %compute reactance points according to model parameters
f_clean=f;
R_clean=R;
X_clean=X;
end
