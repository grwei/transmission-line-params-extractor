function [s_params,rlgc_struct] = rlgc2s_t(resistance, inductance, conductance, capacitance, linelength, freq, z0)
% RLGC2S Converts RLGC-parameters of transmission lines to S-parameters
%   [S_PARAMS,GAMMA] = RLGC2S_T(RESISTANCE,INDUCTANCE,CONDUCTANCE,CAPACITANCE,LINELENGTH,FREQ,Z0)
%   converts the RLGC-matrices of a transmission line into its scattering
%   parameters
%   
%   RESISTANCE is a real N-by-N-by-M Resistance matrix (ohm/m) 
%   INDUCTANCE is a real N-by-N-by-M Inductance matrix (H/m) 
%   CONDUCTANCE is a real N-by-N-by-M Conductance matrix (S/m) 
%   CAPACITANCE is a real N-by-N-by-M Capacitance matrix (F/m) 
%   LINELENGTH is the length of the transmission line 
%   FREQ is the Mx1 frequency vector 
%   Z0 is the reference impedance, the default is 50 ohms.
%   
%   The output S_PARAMS is a complex 2N-by-2N-by-M array, where M is the
%   number of frequency points at which the RLGC-matrices are specified and
%   N is the number of transmission lines.
% 
% Properties of RLGC matrices
%   R,L,G,C matrices are symmetric. 
%   Diagonal terms of L and C are positive, non-zero.
% 	Diagonal terms of R and G are non-negative (can be zero).
% 	Off-diagonal terms of the L matrix are non-negative. 
%   Off-diagonal terms of C and G matrices are non-positive. 
%   Off-diagonal terms of all matrices can be zero.
% 
%   See also RLGC2S, ABCD2S, S2Y, S2Z, S2H, Y2ABCD, Z2ABCD, H2ABCD, S2RLGC

narginchk(6,7);

if nargin < 7
    z0 = 50;
else
    if(~isscalar(z0))
        error(message('rf:rlgc2s:InvalidInputZ0NonScalar'));
    end
    if(isnan(z0) || isinf(z0))
        error(message('rf:rlgc2s:InvalidInputZ0NanInf'));
    end
end

if ~all(imag(z0)==0)
    error(message('rflib:shared:ComplexZ0'))
end

freqpts  = size(freq(:),1);         % Number of frequency points
num_lines = size(resistance,1);     % Number of transmission lines
% Allocate memory
s_params            = zeros(2*num_lines, 2*num_lines, freqpts); % 
z0_matrix           = z0 * eye(2*num_lines, 2*num_lines);       % reference characteristic impedance matrix
transform_matrix    = zeros(num_lines, num_lines, freqpts);     % 
square_of_gammaEig  = zeros(num_lines, freqpts);                % 
gammaEig            = square_of_gammaEig; % eigen-mode propagation constant
gamma               = transform_matrix;   % complex propagation constant matrix
Zc                  = transform_matrix;   % characteristic impedance matrix
A                   = transform_matrix;   % A term of ABCD(chain) matrix
B                   = transform_matrix;   % B term of ABCD(chain) matrix
C                   = transform_matrix;   % C term of ABCD(chain) matrix
D                   = transform_matrix;   % D term of ABCD(chain) matrix
Z11                 = transform_matrix;   % 11 term of Z-params matrix
Z12                 = transform_matrix;   % 12 term of Z-params matrix
Z21                 = transform_matrix;   % 21 term of Z-params matrix
Z22                 = transform_matrix;   % 22 term of Z-params matrix
Z_params            = s_params;           % Z-params matrix
cosh_gamma_length   = transform_matrix;   % cosh(gamma * linelength)
sinh_gamma_length   = transform_matrix;   % sinh(gamma * linelength)

for freqidx=1:freqpts    
    Z = complex(resistance(:,:,freqidx), 2*pi*freq(freqidx).*inductance(:,:,freqidx));   % p.u.l. impedance
    Y = complex(conductance(:,:,freqidx), 2*pi*freq(freqidx).*capacitance(:,:,freqidx)); % p.u.l. admittance

    %%% 1. Calculate Complex Propagation Constants
    % [V,D] = eig(A,'vector') returns colume vector D of eigenvalues 
    % and matrix V whose columns are the corresponding right eigenvectors,
    % so that A*V = V*diag(D).
    % The eigenvectors in V are normalized so that the 2-norm of each is 1.
    [transform_matrix(:,:,freqidx),square_of_gammaEig(:,freqidx)] =     ...
                                                       eig(Z * Y,'vector');                  
    % In general, a complex number has two square roots. 
    % We choose the one with a non-negative real part, 
    % which is also the one returned by the MATLAB function sqrt().
    gammaEig(:,freqidx) = sqrt(square_of_gammaEig(:,freqidx));
    gamma(:,:,freqidx) = transform_matrix(:,:,freqidx)              *   ...
                         diag(gammaEig(:,freqidx))                  /   ...
                         transform_matrix(:,:,freqidx);
    %%% 2. Calculate Characteristic Impedance matrix
    Zc(:,:,freqidx) = gamma(:,:,freqidx) \ Z;
    %%% 3. Calculate the ABCD matices [Paul2007]
    cosh_gamma_length(:,:,freqidx)                                  =   ... 
                      transform_matrix(:,:,freqidx)                 *   ...
                      diag(cosh(gammaEig(:,freqidx) * linelength))  /   ...
                      transform_matrix(:,:,freqidx);
    sinh_gamma_length(:,:,freqidx)                                  =   ... 
                      transform_matrix(:,:,freqidx)                 *   ...
                      diag(sinh(gammaEig(:,freqidx) * linelength))  /   ...
                      transform_matrix(:,:,freqidx);
    A(:,:,freqidx) = cosh_gamma_length(:,:,freqidx);
    B(:,:,freqidx) = sinh_gamma_length(:,:,freqidx) * Zc(:,:,freqidx);
    C(:,:,freqidx) = Zc(:,:,freqidx) \ sinh_gamma_length(:,:,freqidx);
    D(:,:,freqidx) = Zc(:,:,freqidx)                                \   ...
                     cosh_gamma_length(:,:,freqidx)                 *   ...
                     Zc(:,:,freqidx);
    %%% test [Sampath2008] ABCD<-RLGC
    % It seems that this formula is bad.
    % should be removed when released
%     D(:,:,freqidx) = cosh_gamma_length(:,:,freqidx);
%     B(:,:,freqidx) = Zc(:,:,freqidx) * sinh_gamma_length(:,:,freqidx);
%     C(:,:,freqidx) = sinh_gamma_length(:,:,freqidx) / Zc(:,:,freqidx);
%     A(:,:,freqidx) = Zc(:,:,freqidx)                                *   ...
%                      cosh_gamma_length(:,:,freqidx)                 /   ...
%                      Zc(:,:,freqidx);
    %%% 4. Convert ABCD matrix to Z-parameter(impedance) matrix[Reveyrand2018]
    Z11(:,:,freqidx) = A(:,:,freqidx) / C(:,:,freqidx);
    Z12(:,:,freqidx) = Z11(:,:,freqidx) * D(:,:,freqidx) - B(:,:,freqidx);
    Z21(:,:,freqidx) = eye(num_lines,num_lines) / C(:,:,freqidx);
    Z22(:,:,freqidx) = Z21(:,:,freqidx) * D(:,:,freqidx);
    Z_params(:,:,freqidx) = [Z11(:,:,freqidx),Z12(:,:,freqidx);         ...
                             Z21(:,:,freqidx),Z22(:,:,freqidx)];
    %%% 5. Convert Z-parameter(impedance) matrix to S matrix [Reveyrand2018]
%     s_params(:,:,freqidx) = (Z_params(:,:,freqidx) + z0_matrix)       \ ...
%                             (Z_params(:,:,freqidx) - z0_matrix);
    %%% test [Reveyrand2018]: S->Z
    % It seems that both formula works well.
    % should be removed when released
    s_params(:,:,freqidx) = (Z_params(:,:,freqidx) - z0_matrix)       / ...
                            (Z_params(:,:,freqidx) + z0_matrix);
end

% The following output struct can be removed when released.
rlgc_struct = struct('R',resistance,'L',inductance,'G',conductance,'C',capacitance,'alpha',real(gamma),'beta',imag(gamma),'Zc',Zc);

end
