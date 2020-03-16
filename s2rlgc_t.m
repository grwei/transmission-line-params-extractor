function output = s2rlgc_t(s_params,linelength,freq,z0,port_reorder)
%S2RLGC Converts S-parameters of a transmission line to RLGC-parameters
%   OUTPUT = S2RLGC(S_PARAMS, LINELENGTH, FREQ, Z0, PORT_REORDER) converts
%   the scattering parameters S_PARAMS of a transmission line into
%   RLGC-matrices. 
%
%   S_PARAMS is a complex 2N-by-2N-by-M array, where M is the number of
%   frequency points at which the S-parameters are specified and N is
%   the number of transmission lines.
%   LINELENGTH is the length of the transmission line 
%   FREQ is a real Mx1 frequency vector
%   Z0 is the reference impedance, the default is 50 ohms.
%   PORT_REORDER is a 2Nx1 vector indicating the input and output ports
%   [IP... OP ...]. Ports on one side of the lines are numbered first 
%   followed by the ports on the other side; thus, if i refers to a port 
%   of a line at z=0, then the port on the other side is numbered  N+i.
%
%   The outputs are per unit length transmission line parameters
%   OUTPUT.R is a real N-by-N-by-M Resistance matrix (ohm/m)
%   OUTPUT.L is a real N-by-N-by-M Inductance matrix (H/m)
%   OUTPUT.C is a real N-by-N-by-M Capacitance matrix (F/m)
%   OUTPUT.G is a real N-by-N-by-M Conductance matrix (S/m)
%   OUTPUT.Zc is a complex N-by-N-by-M Characteristic line impedance(ohm)
%   OUTPUT.alpha is a real N-by-N-by-M attenuation constant (Nepers/m)
%   OUTPUT.beta is a real N-by-N-by-M phase constant (radians/m)

%% Input validity check and initialization

narginchk(3,5)

% freq should be a column vector
if isvector(freq) && isrow(freq)
    freq = transpose(freq);
end

freqpts  = size(freq(:),1);        % Number of frequency points
numLines = size(s_params,1)/2;     % Number of transmission lines

% the default reference impedance is 50 ohms
if nargin < 4
    z0 = 50;
end

% Reorder the port
if nargin < 5
    port_reorder = [];
end

if ~isempty(port_reorder)
    s_params = snp2smp(s_params,z0,port_reorder);
end

%% Convert S-parameters to ABCD-Parameters[Sampath2008]

Z_params = zeros(2*numLines,2*numLines,freqpts);    % impedance matrix
TA = zeros(numLines,numLines,freqpts);              % the A term of transmission(ABCD) matrix:[A B;C D]
% The B, C, and D terms are not actually used.
TB = TA;
TC = TA;
TD = TA;
I = eye(2*numLines,2*numLines);                     % 2Nx2N identity matrix

for idx=1:freqpts
    Z_params(:,:,idx) = z0*(I+s_params(:,:,idx)) / (I-s_params(:,:,idx));
    TA(:,:,idx) =   Z_params(1:numLines,1:numLines,idx) /               ...
                    Z_params(numLines+1:end,1:numLines,idx);
    % The B, C, and D terms are not actually used.
    TB(:,:,idx) =   TA * Z_params(numLines+1:end,numLines+1:end) -      ...
                    Z_params(1:numLines,numLines+1:end,idx);
    TC(:,:,idx) =   I / (I-s_params(:,:,idx));
    TD(:,:,idx) =   TC(:,:,idx) * Z_params(numLines+1:end,numLines+1:end);
end

%% Extract Complex Propagation Constants while preserving relative Eigenvalue position

% Eigenvalue decomposition(similarity transformation) of TA
eigVal = zeros(numLines,freqpts);           % Eigenvalues of TA
eigVec = zeros(numLines,numLines,freqpts);  % Columns are the corresponding Right Eigenvectors of TA
for idx = 1:freqpts
    % [V,D] = eig(TA,'vector') --> TA*diag(V) = diag(V)*D 
    % Note: The eigenvectors should be normalized[Braunisch1998] 
    % so that the 2-norm of each is 1, which is done by default 
    % by the MATLAB function eig().
    [eigVec(:,:,idx),eigVal(:,idx)] = eig(TA(:,:,idx),'vector');
end

% Adjust the order of eigenvalues and eigenvectors [Braunisch1998, Chu2015]
prodTable = nan(numLines,numLines,freqpts);     % Hermitian Inner Product recoder.
newIndex  = nan(numLines,freqpts);              % Correct order of eigVal and corresponding eigVec
for freqidx = 2:freqpts                         % Index of frequency point
    prodTable(:,:,freqidx) =  ctranspose(eigVec(:,:,freqidx)) * eigVec(:,:,freqidx-1);
    % ²¹³ä×¢ÊÍ! --grwei, 20200316
    newIndex(:,freqidx) = max(prodTable(:,:,freqidx),[],2);
end

%% 

output = struct('R',R,'L',L,'G',G,'C',C,'alpha',alpha,'beta',beta,'Zc',Zc);

end

