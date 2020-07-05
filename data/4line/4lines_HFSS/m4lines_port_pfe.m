function [Y] = m4lines_port_pfe(arg1,arg2,arg3)
% This function takes the coefficients from m4lines_port_pfe and
% either calculates the Y matrix at a given frequency or
% calculates a vector of the Y_ij over a frequency range.
% Frequency is assumed to be in GHz.
% 
% Example 1: Y = m4lines_port_pfe(f) returns
%            at the SINGLE frequency f.
% Example 2: Y_12 = m4lines_port_pfe(1,2,f) returns
%            the element Y_{1,2} at the VECTOR frequency f.
% 
% This M-File brought to you by ElectronicsDesktop

if nargin == 1
  Y = CalculateYMatrix(arg1);
elseif nargin == 3
  Y = CalculateYVector(arg1,arg2,arg3);
else
  Y = [];
  disp(' Incorrect Number Of Arguments to m4lines_port_pfe ') 
end

function [Y] = CalculateYMatrix(freq)
% This function takes the coefficients from coeff.m and
% and calculates the Y matrix at a given frequency.
m4lines_port_pfe_sol
s = i * freq/70;
N = 8;
Y = zeros(N);

% set the constant term
Y = A1;

% scroll through the poles
for ii = 1:(N-1)
end

% now the pole at infinity
Y = Y + A2.*s;


function [Y] = CalculateYVector(ii,jj,freq)
% This function takes the coefficients from coeff.m and
% and calculates the Y matrix at a given frequency.
m4lines_port_pfe_sol
s = i * freq/70;
N = 8;
Y = zeros(size(freq));

% set the constant term
Y = A1(ii,jj)*ones(size(freq));

% scroll through the poles
% now the pole at infinity
Y = Y + A2(ii,jj).*s;


