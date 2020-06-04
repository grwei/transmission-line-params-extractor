function [Y] = m4lines_HFSS_pfe(arg1,arg2,arg3)
% This function takes the coefficients from m4lines_HFSS_pfe and
% either calculates the Y matrix at a given frequency or
% calculates a vector of the Y_ij over a frequency range.
% Frequency is assumed to be in GHz.
% 
% Example 1: Y = m4lines_HFSS_pfe(f) returns
%            at the SINGLE frequency f.
% Example 2: Y_12 = m4lines_HFSS_pfe(1,2,f) returns
%            the element Y_{1,2} at the VECTOR frequency f.
% 
% This M-File brought to you by ElectronicsDesktop

if nargin == 1
  Y = CalculateYMatrix(arg1);
elseif nargin == 3
  Y = CalculateYVector(arg1,arg2,arg3);
else
  Y = [];
  disp(' Incorrect Number Of Arguments to m4lines_HFSS_pfe ') 
end

function [Y] = CalculateYMatrix(freq)
% This function takes the coefficients from coeff.m and
% and calculates the Y matrix at a given frequency.
m4lines_HFSS_pfe_sol
s = i * freq/70;
N = 8;
Y = zeros(N);

% set the constant term
Y = A1;

% scroll through the poles
for ii = 1:(N-1)
   Y = Y + A2./(s - P(1));
   Y = Y + A3./(s - P(2));
   Y = Y + A4./(s - P(3));
   Y = Y + A5./(s - P(4));
   Y = Y + A6./(s - P(5));
   Y = Y + A7./(s - P(6));
   Y = Y + A8./(s - P(7));
   Y = Y + A9./(s - P(8));
   Y = Y + A10./(s - P(9));
   Y = Y + A11./(s - P(10));
   Y = Y + A12./(s - P(11));
   Y = Y + A13./(s - P(12));
   Y = Y + A14./(s - P(13));
   Y = Y + A15./(s - P(14));
   Y = Y + A16./(s - P(15));
   Y = Y + A17./(s - P(16));
   Y = Y + A18./(s - P(17));
   Y = Y + A19./(s - P(18));
   Y = Y + A20./(s - P(19));
   Y = Y + A21./(s - P(20));
   Y = Y + A22./(s - P(21));
   Y = Y + A23./(s - P(22));
   Y = Y + A24./(s - P(23));
   Y = Y + A25./(s - P(24));
   Y = Y + A26./(s - P(25));
   Y = Y + A27./(s - P(26));
   Y = Y + A28./(s - P(27));
   Y = Y + A29./(s - P(28));
end

% now the pole at infinity
Y = Y + A30.*s;


function [Y] = CalculateYVector(ii,jj,freq)
% This function takes the coefficients from coeff.m and
% and calculates the Y matrix at a given frequency.
m4lines_HFSS_pfe_sol
s = i * freq/70;
N = 8;
Y = zeros(size(freq));

% set the constant term
Y = A1(ii,jj)*ones(size(freq));

% scroll through the poles
   Y = Y + A2(ii,jj)./(s - P(1));
   Y = Y + A3(ii,jj)./(s - P(2));
   Y = Y + A4(ii,jj)./(s - P(3));
   Y = Y + A5(ii,jj)./(s - P(4));
   Y = Y + A6(ii,jj)./(s - P(5));
   Y = Y + A7(ii,jj)./(s - P(6));
   Y = Y + A8(ii,jj)./(s - P(7));
   Y = Y + A9(ii,jj)./(s - P(8));
   Y = Y + A10(ii,jj)./(s - P(9));
   Y = Y + A11(ii,jj)./(s - P(10));
   Y = Y + A12(ii,jj)./(s - P(11));
   Y = Y + A13(ii,jj)./(s - P(12));
   Y = Y + A14(ii,jj)./(s - P(13));
   Y = Y + A15(ii,jj)./(s - P(14));
   Y = Y + A16(ii,jj)./(s - P(15));
   Y = Y + A17(ii,jj)./(s - P(16));
   Y = Y + A18(ii,jj)./(s - P(17));
   Y = Y + A19(ii,jj)./(s - P(18));
   Y = Y + A20(ii,jj)./(s - P(19));
   Y = Y + A21(ii,jj)./(s - P(20));
   Y = Y + A22(ii,jj)./(s - P(21));
   Y = Y + A23(ii,jj)./(s - P(22));
   Y = Y + A24(ii,jj)./(s - P(23));
   Y = Y + A25(ii,jj)./(s - P(24));
   Y = Y + A26(ii,jj)./(s - P(25));
   Y = Y + A27(ii,jj)./(s - P(26));
   Y = Y + A28(ii,jj)./(s - P(27));
   Y = Y + A29(ii,jj)./(s - P(28));
% now the pole at infinity
Y = Y + A30(ii,jj).*s;


