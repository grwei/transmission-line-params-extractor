function [Y] = m16lines_HFSS_pfe(arg1,arg2,arg3)
% This function takes the coefficients from m16lines_HFSS_pfe and
% either calculates the Y matrix at a given frequency or
% calculates a vector of the Y_ij over a frequency range.
% Frequency is assumed to be in GHz.
% 
% Example 1: Y = m16lines_HFSS_pfe(f) returns
%            at the SINGLE frequency f.
% Example 2: Y_12 = m16lines_HFSS_pfe(1,2,f) returns
%            the element Y_{1,2} at the VECTOR frequency f.
% 
% This M-File brought to you by ElectronicsDesktop

if nargin == 1
  Y = CalculateYMatrix(arg1);
elseif nargin == 3
  Y = CalculateYVector(arg1,arg2,arg3);
else
  Y = [];
  disp(' Incorrect Number Of Arguments to m16lines_HFSS_pfe ') 
end

function [Y] = CalculateYMatrix(freq)
% This function takes the coefficients from coeff.m and
% and calculates the Y matrix at a given frequency.
m16lines_HFSS_pfe_sol
s = i * freq/70;
N = 32;
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
   Y = Y + A30./(s - P(29));
   Y = Y + A31./(s - P(30));
   Y = Y + A32./(s - P(31));
   Y = Y + A33./(s - P(32));
   Y = Y + A34./(s - P(33));
   Y = Y + A35./(s - P(34));
   Y = Y + A36./(s - P(35));
   Y = Y + A37./(s - P(36));
   Y = Y + A38./(s - P(37));
   Y = Y + A39./(s - P(38));
   Y = Y + A40./(s - P(39));
   Y = Y + A41./(s - P(40));
   Y = Y + A42./(s - P(41));
   Y = Y + A43./(s - P(42));
   Y = Y + A44./(s - P(43));
   Y = Y + A45./(s - P(44));
   Y = Y + A46./(s - P(45));
   Y = Y + A47./(s - P(46));
   Y = Y + A48./(s - P(47));
   Y = Y + A49./(s - P(48));
   Y = Y + A50./(s - P(49));
   Y = Y + A51./(s - P(50));
   Y = Y + A52./(s - P(51));
   Y = Y + A53./(s - P(52));
   Y = Y + A54./(s - P(53));
   Y = Y + A55./(s - P(54));
   Y = Y + A56./(s - P(55));
   Y = Y + A57./(s - P(56));
   Y = Y + A58./(s - P(57));
   Y = Y + A59./(s - P(58));
   Y = Y + A60./(s - P(59));
   Y = Y + A61./(s - P(60));
   Y = Y + A62./(s - P(61));
   Y = Y + A63./(s - P(62));
   Y = Y + A64./(s - P(63));
   Y = Y + A65./(s - P(64));
   Y = Y + A66./(s - P(65));
   Y = Y + A67./(s - P(66));
   Y = Y + A68./(s - P(67));
   Y = Y + A69./(s - P(68));
   Y = Y + A70./(s - P(69));
   Y = Y + A71./(s - P(70));
   Y = Y + A72./(s - P(71));
   Y = Y + A73./(s - P(72));
   Y = Y + A74./(s - P(73));
   Y = Y + A75./(s - P(74));
   Y = Y + A76./(s - P(75));
   Y = Y + A77./(s - P(76));
   Y = Y + A78./(s - P(77));
   Y = Y + A79./(s - P(78));
   Y = Y + A80./(s - P(79));
   Y = Y + A81./(s - P(80));
end

% now the pole at infinity
Y = Y + A82.*s;


function [Y] = CalculateYVector(ii,jj,freq)
% This function takes the coefficients from coeff.m and
% and calculates the Y matrix at a given frequency.
m16lines_HFSS_pfe_sol
s = i * freq/70;
N = 32;
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
   Y = Y + A30(ii,jj)./(s - P(29));
   Y = Y + A31(ii,jj)./(s - P(30));
   Y = Y + A32(ii,jj)./(s - P(31));
   Y = Y + A33(ii,jj)./(s - P(32));
   Y = Y + A34(ii,jj)./(s - P(33));
   Y = Y + A35(ii,jj)./(s - P(34));
   Y = Y + A36(ii,jj)./(s - P(35));
   Y = Y + A37(ii,jj)./(s - P(36));
   Y = Y + A38(ii,jj)./(s - P(37));
   Y = Y + A39(ii,jj)./(s - P(38));
   Y = Y + A40(ii,jj)./(s - P(39));
   Y = Y + A41(ii,jj)./(s - P(40));
   Y = Y + A42(ii,jj)./(s - P(41));
   Y = Y + A43(ii,jj)./(s - P(42));
   Y = Y + A44(ii,jj)./(s - P(43));
   Y = Y + A45(ii,jj)./(s - P(44));
   Y = Y + A46(ii,jj)./(s - P(45));
   Y = Y + A47(ii,jj)./(s - P(46));
   Y = Y + A48(ii,jj)./(s - P(47));
   Y = Y + A49(ii,jj)./(s - P(48));
   Y = Y + A50(ii,jj)./(s - P(49));
   Y = Y + A51(ii,jj)./(s - P(50));
   Y = Y + A52(ii,jj)./(s - P(51));
   Y = Y + A53(ii,jj)./(s - P(52));
   Y = Y + A54(ii,jj)./(s - P(53));
   Y = Y + A55(ii,jj)./(s - P(54));
   Y = Y + A56(ii,jj)./(s - P(55));
   Y = Y + A57(ii,jj)./(s - P(56));
   Y = Y + A58(ii,jj)./(s - P(57));
   Y = Y + A59(ii,jj)./(s - P(58));
   Y = Y + A60(ii,jj)./(s - P(59));
   Y = Y + A61(ii,jj)./(s - P(60));
   Y = Y + A62(ii,jj)./(s - P(61));
   Y = Y + A63(ii,jj)./(s - P(62));
   Y = Y + A64(ii,jj)./(s - P(63));
   Y = Y + A65(ii,jj)./(s - P(64));
   Y = Y + A66(ii,jj)./(s - P(65));
   Y = Y + A67(ii,jj)./(s - P(66));
   Y = Y + A68(ii,jj)./(s - P(67));
   Y = Y + A69(ii,jj)./(s - P(68));
   Y = Y + A70(ii,jj)./(s - P(69));
   Y = Y + A71(ii,jj)./(s - P(70));
   Y = Y + A72(ii,jj)./(s - P(71));
   Y = Y + A73(ii,jj)./(s - P(72));
   Y = Y + A74(ii,jj)./(s - P(73));
   Y = Y + A75(ii,jj)./(s - P(74));
   Y = Y + A76(ii,jj)./(s - P(75));
   Y = Y + A77(ii,jj)./(s - P(76));
   Y = Y + A78(ii,jj)./(s - P(77));
   Y = Y + A79(ii,jj)./(s - P(78));
   Y = Y + A80(ii,jj)./(s - P(79));
   Y = Y + A81(ii,jj)./(s - P(80));
% now the pole at infinity
Y = Y + A82(ii,jj).*s;


