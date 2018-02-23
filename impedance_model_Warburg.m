%Determines the theoretical impedance model taking
%into consideration a double layer capacitor and a Warburg impedance element
%--------------------------------------------------------------------------
%Usage:
%function z=modelWarburg(p,f)
%All parameters should be defined in SI units
%--------------------------------------------------------------------------
%Inputs:
%p: parameter vector containing the values to fit into the model
%   -p(1)=Rint: Resistance associated with the double layer [ohms]
%   -p(2)=Rmedium: Resistance associated with the medium [ohms]
%   -p(3)=C:  Parameter that defines the double layer capacitance effect [1/s] (this is not an actual capacitor, see my notes)
%   -p(4)=n:  Parameter that defines the double layer capacitance effect 
%   $\frac{1}{(j\omega)^nC}$
%   -p(5)=A: Warburg coefficient
%-f: frequency values (in a one dimensional column vector) [Hz]
%Outputs:
%Calculated impedance values
%z(1,:)=Calculated resistance R [ohms] (real part of the impedance)
%z(2,:)=Calculated Reactance X [ohms] (imaginary part of impedance)
%
%Author:Leopoldo Cendejas Zaragoza 2018

function z=modelWarburg(p,f)
  %extract parameters from vector
  Rint=p(1);
  Rmedium=p(2);
  C=p(3);
  n=p(4);
  A=p(5);
  
  z=zeros(2, length(f)); %initialize model output 
  %Compute model 
  w=2*pi*f; %angular frequency
  
  %Warburg impedance
  z_w=A./sqrt(w)+A./(1i*sqrt(w));
  %z_w=A./sqrt(1i*w);
  %Double layer capacitance
  z_c=1./((1i*w).^n*C);
  
  %Compute model impedance
  z_eq=((Rint+z_w).*z_c)./(Rint+z_w+z_c)+Rmedium;
  
  %Extract real and imaginary components to use them in 
  z(1,:)= real(z_eq); %value of resistance
  z(2,:)= imag(z_eq); %value of reactance
end