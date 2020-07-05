function [Y] = m16lines_port_pfe(arg1,arg2,arg3)
% This function takes the coefficients from m16lines_port_pfe and
% either calculates the Y matrix at a given frequency or
% calculates a vector of the Y_ij over a frequency range.
% Frequency is assumed to be in GHz.
% 
% Example 1: Y = m16lines_port_pfe(f) returns
%            at the SINGLE frequency f.
% Example 2: Y_12 = m16lines_port_pfe(1,2,f) returns
%            the element Y_{1,2} at the VECTOR frequency f.
% 
% This M-File brought to you by ElectronicsDesktop

if nargin == 1
  Y = CalculateYMatrix(arg1);
elseif nargin == 3
  Y = CalculateYVector(arg1,arg2,arg3);
else
  Y = [];
  disp(' Incorrect Number Of Arguments to m16lines_port_pfe ') 
end

function [Y] = CalculateYMatrix(freq)
% This function takes the coefficients from coeff.m and
% and calculates the Y matrix at a given frequency.
m16lines_port_pfe_sol
s = i * freq/70;
N = 32;
Y = zeros(N);

% set the constant term
Y = A1;

% scroll through the poles
for ii = 1:(N-1)
   Y = Y + A2./(s - P(1));
   Y = Y + A3./(s - P(2));
end

% now the pole at infinity
Y = Y + A4.*s;


function [Y] = CalculateYVector(ii,jj,freq)
% This function takes the coefficients from coeff.m and
% and calculates the Y matrix at a given frequency.
m16lines_port_pfe_sol
s = i * freq/70;
N = 32;
Y = zeros(size(freq));

% set the constant term
Y = A1(ii,jj)*ones(size(freq));

% scroll through the poles
   Y = Y + A2(ii,jj)./(s - P(1));
   Y = Y + A3(ii,jj)./(s - P(2));
% now the pole at infinity
Y = Y + A4(ii,jj).*s;


