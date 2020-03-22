function output = s2rlgc_t(s_params,linelength,freq,z0,port_reorder)
%S2RLGC Converts S-parameters of a transmission line to RLGC-parameters
%   OUTPUT = S2RLGC_T(S_PARAMS, LINELENGTH, FREQ, Z0, PORT_REORDER) converts
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

%% Configuration

debugFlag = true;       % Debug mode
EigSortMethod = 0;      % 0(default)- abs(real(...))  1- abs(...)

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
TB = TA;
% The C and D terms are not actually used.
TC = TA;
TD = TA;
I = eye(2*numLines,2*numLines);                     % 2Nx2N identity matrix

for idx=1:freqpts
    Z_params(:,:,idx) = z0*(I+s_params(:,:,idx))                /       ...
        (I-s_params(:,:,idx));
    TA(:,:,idx) =   Z_params(1:numLines,1:numLines,idx)         /       ...
        Z_params(numLines+1:end,1:numLines,idx);
    TB(:,:,idx) =   TA(:,:,idx)                                 *       ...
        Z_params(numLines+1:end,numLines+1:end,idx)             -       ...
        Z_params(1:numLines,numLines+1:end,idx);
    % The C and D terms are not actually used.
    TC(:,:,idx) =   eye(numLines,numLines)                      /       ...
        Z_params(numLines+1:end,1:numLines,idx);
    TD(:,:,idx) =   TC(:,:,idx)                                 *       ...
        Z_params(numLines+1:end,numLines+1:end,idx);
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

% Adjust the order of Eigenvalues and Eigenvectors [Braunisch1998, Chu2015]
prodTable   = nan(numLines,numLines,freqpts);   % Hermitian Inner Product recoder.
CorrectPos  = nan(numLines,freqpts);            % Correct position of eigVal and corresponding eigVec
newIndex    = nan(numLines,freqpts);            % New index of eigVal and corresponding eigVec
for freqidx = 2:freqpts                         % Index of frequency point
    % For each Eigenvector at the current frequency, calculate the
    % Hermitian inner product of each Eigenvector at the previous
    % frequency. Since the Eigenvector have been normalized in advance,
    % there should be only one of these inner products closest to 1,
    % which indicates that the two Eigenvectors involved correspond to
    % the Eigenvalue of the same position.[Braunisch1998]
    prodTable(:,:,freqidx) =  ctranspose(eigVec(:,:,freqidx)) *     ...
        eigVec(:,:,freqidx-1);
    
    % Determine the correct position by Hermitian Inner Product.
    % WARNNING! Duplicate serial numbers may be generated in extreme
    % cases.
    
    % !!!Two choices: 0. max(abs(real(...)))   1. max(abs(...))
    switch(EigSortMethod)
        case 0
            [~,CorrectPos(:,freqidx)] = max(abs(real(prodTable(:,:,freqidx))),[],2);
        case 1
            [~,CorrectPos(:,freqidx)] = max(abs(prodTable(:,:,freqidx)),[],2);
    end
    [~,newIndex(:,freqidx)] = sort(CorrectPos(:,freqidx));
    eigVec(:,:,freqidx)     = eigVec(:,newIndex(:,freqidx),freqidx);
    eigVal(:,freqidx)       = eigVal(newIndex(:,freqidx),freqidx);
end

%% Extract Attenuation Constants and Unwrapped Phase Constants

gammaLenEigWrap = acosh(eigVal);  % Principle Value of gammaEig*linelength:
% Real part is non-negative, imag part in [-pi,pi]
betaLenEigWrapDiff(:,2:freqpts) = diff(imag(gammaLenEigWrap),1,2);
discontCount = cumsum(abs(betaLenEigWrapDiff) > pi,2);

%% !!!Phase-unwrapping algorithm is unreliable near singular frequency! --grwei,20200318
%%Maybe position-tracking algorithm instead???
% discontCount(:,491:end) = 1;

%%
betaLenEigUnwrap = imag(gammaLenEigWrap) + 2*pi*discontCount;
gammaEigUnwrap = complex(real(gammaLenEigWrap),betaLenEigUnwrap) / linelength;

gamma = nan(numLines,numLines,freqpts);
for idx = 1:freqpts
    gamma(:,:,idx) = eigVec(:,:,idx)        *                           ...
        diag(gammaEigUnwrap(:,idx))         /                           ...
        eigVec(:,:,idx);
end

alpha = real(gamma);
beta  = imag(gamma);

%% Extract Characteristic Impedance Matrix

Zc = zeros(numLines,numLines,freqpts);  % Characteristic impedance
% The phase contant need not be unwrapped to compute Zc[]
for idx = 1:freqpts
    Zc(:,:,idx) = eigVec(:,:,idx)               *               ...
        diag(1./sinh(gammaLenEigWrap(:,idx)))   /               ...
        eigVec(:,:,idx)                         *               ...
        TB(:,:,idx);
end

%% Extract RLGC

R = zeros(numLines,numLines,freqpts);       % Resistance matrix
L = R;                                      % Inductance matrix
C = R;                                      % Capacitance matrix
G = R;                                      % Conductance matrix
Z_pul = R;                                  % p.u.l impedance matrix
Y_pul = R;                                  % p.u.l admittance matrix

for idx = 1:freqidx
    Z_pul(:,:,idx) =  gamma(:,:,idx) * Zc(:,:,idx);
    Y_pul(:,:,idx) =  Zc(:,:,idx) \ gamma(:,:,idx);
    R(:,:,idx) = real(Z_pul(:,:,idx));
    L(:,:,idx) = imag(Z_pul(:,:,idx)) / (2*pi*freq(idx));
    G(:,:,idx) = real(Y_pul(:,:,idx));
    C(:,:,idx) = imag(Y_pul(:,:,idx)) / (2*pi*freq(idx));
end

%%

output = struct('R',R,'L',L,'G',G,'C',C,'alpha',alpha,'beta',beta,'Zc',Zc);

% end

%%
%The following code snippets are for debugging purposes only.

if debugFlag == false
    return
end
%%

figure('Name','prodTable(Real-part)')
sgtitle('Hermitian Inner Product(Real-part)')
sz = 10;
for idx_cur = 1:numLines
    subplot(2,ceil(numLines/2),idx_cur)
    scatter(freq,abs(real(squeeze(prodTable(idx_cur,1,:)))),sz,'o');
    hold on
    for idx_pre = 2:numLines
        scatter(freq,abs(real(squeeze(prodTable(idx_cur,idx_pre,:)))),sz,'o');
    end
    hold off
    txt = cell(1,numLines);
    for idx = 1:numLines
        txt{1,idx} = ['prePos-',sprintf('%u',idx)];
    end
    legend(txt,'Location','best')
    legend('boxoff')
    title(['curPos-',sprintf('%u',idx_cur)])
    xlabel('Frequency(Hz)')
    ylabel('Abs(Real-part)')
    %     ylim([0 1.1])
end

%%


figure('Name','prodTable(Magnitude)')
sgtitle('Hermitian Inner Product(Magnitude)')
sz = 10;
for idx_cur = 1:numLines
    subplot(2,ceil(numLines/2),idx_cur)
    scatter(freq,abs(squeeze(prodTable(idx_cur,1,:))),sz,'o');
    hold on
    for idx_pre = 2:numLines
        scatter(freq,abs(squeeze(prodTable(idx_cur,idx_pre,:))),sz,'o');
    end
    hold off
    txt = cell(1,numLines);
    for idx = 1:numLines
        txt{1,idx} = ['prePos-',sprintf('%u',idx)];
    end
    legend(txt,'Location','best')
    legend('boxoff')
    title(['curPos-',sprintf('%u',idx_cur)])
    xlabel('Frequency(Hz)')
    ylabel('Value(Magnitude)')
    %     ylim([0 1.1])
end


%%

figure('Name','gammaLenEigWrap')
sgtitle('Gamma(wrap) of each eigen mode')
subplot(121)
plot(freq,real(gammaLenEigWrap(1,:))/linelength)
hold on
for idx = 2:numLines
    plot(freq,real(gammaLenEigWrap(idx,:))/linelength)
end
hold off
grid on
xlim([4.5e9 5.5e9])
xlabel('Frequency(Hz)')
ylabel('\alpha(Np/m)')
txt = cell(1,numLines);
for idx = 1:numLines
    txt{1,idx} = ['\alpha_',sprintf('%u',idx)];
end
legend(txt,'Location','best','NumColumns',2)
legend('boxoff')
title('\alpha')

subplot(122)
plot(freq,imag(gammaLenEigWrap(1,:)))
hold on
for idx = 2:numLines
    plot(freq,imag(gammaLenEigWrap(idx,:)))
end
hold off
grid on
xlim([4.5e9 5.5e9])
xlabel('Frequency(Hz)')
ylabel('\betaL(rad)')
txt = cell(1,numLines);
for idx = 1:numLines
    txt{1,idx} = ['\beta_',sprintf('%u',idx)];
end
legend(txt,'Location','best','NumColumns',2)
legend('boxoff')
title('\betaL')

%%

figure('Name','Zc(Real part)')
sgtitle('Charateristic Impedance Matrix(Real part)')
for idx = 1:numLines
    subplot(2,ceil(numLines/2),idx)
    plot(freq,real(squeeze(Zc(idx,1,:))))
    hold on
    for idx_col = 2:numLines
        plot(freq,real(squeeze(Zc(idx,idx_col,:))))
    end
    hold off
    txt = cell(1,numLines);
    for idx_col = 1:numLines
        txt{1,idx_col} = ['Zc_{',sprintf('%u',idx),sprintf('%u',idx_col),'}'];
    end
    legend(txt,'Location','best','NumColumns',2)
    legend('boxoff')
    xlabel('Frequency(Hz)')
    ylabel('Real-part(Ohms)')
    title(['Zc_{',sprintf('%u',idx),'X}'])
end

%
figure('Name','Zc(Imag part)')
sgtitle('Charateristic Impedance Matrix(Imag part)')
for idx = 1:numLines
    subplot(2,ceil(numLines/2),idx)
    plot(freq,imag(squeeze(Zc(idx,1,:))))
    hold on
    for idx_col = 2:numLines
        plot(freq,imag(squeeze(Zc(idx,idx_col,:))))
    end
    hold off
    txt = cell(1,numLines);
    for idx_col = 1:numLines
        txt{1,idx_col} = ['Zc_{',sprintf('%u',idx),sprintf('%u',idx_col),'}'];
    end
    legend(txt,'Location','best','NumColumns',2)
    legend('boxoff')
    xlabel('Frequency(Hz)')
    ylabel('Imag-part(Ohms)')
    title(['Zc_{',sprintf('%u',idx),'X}'])
end

%%

figure('Name','gammaEigUnwrap')
sgtitle('Gamma(unwrap) of each eigen mode')
subplot(121)
plot(freq,real(gammaEigUnwrap(1,:)))
hold on
for idx = 2:numLines
    plot(freq,real(gammaEigUnwrap(idx,:)))
end
hold off
grid on
xlim([4.5e9 5.5e9])
xlabel('Frequency(Hz)')
ylabel('\alpha(Np/m)')
txt = cell(1,numLines);
for idx = 1:numLines
    txt{1,idx} = ['\alpha_',sprintf('%u',idx)];
end
legend(txt,'Location','best','NumColumns',2)
legend('boxoff')
title('\alpha')

subplot(122)
plot(freq,imag(gammaEigUnwrap(1,:)))
hold on
for idx = 2:numLines
    plot(freq,imag(gammaEigUnwrap(idx,:)))
end
hold off
grid on
xlim([4.5e9 5.5e9])
xlabel('Frequency(Hz)')
ylabel('\betaL(rad)')
txt = cell(1,numLines);
for idx = 1:numLines
    txt{1,idx} = ['\beta_',sprintf('%u',idx)];
end
legend(txt,'Location','best','NumColumns',2)
legend('boxoff')
title('\beta')

%%

figure('Name','RLGC matrix')
sgtitle({'Extracted RLGC matrix'})
% R
subplot(221)
plot(freq,squeeze(R(1,1,:)))
hold on
for idx = 2:numLines
    plot(freq,squeeze(R(1,idx,:)))
end
hold off
grid on
xlabel('Frequency(Hz)')
ylabel('Value(Ohms/m)')
txt = cell(1,numLines);
for idx = 1:numLines
    txt{1,idx} = sprintf('R1%u',idx);
end
legend(txt,'Location','best','NumColumns',2)
legend('boxoff')
title('R matrix')
% L
subplot(222)
plot(freq,squeeze(L(1,1,:)))
hold on
for idx = 2:numLines
    plot(freq,squeeze(L(1,idx,:)))
end
hold off
grid on
xlabel('Frequency(Hz)')
ylabel('Value(H/m)')
txt = cell(1,numLines);
for idx = 1:numLines
    txt{1,idx} = sprintf('L1%u',idx);
end
legend(txt,'Location','best','NumColumns',2)
legend('boxoff')
title('L matrix')
% G
subplot(223)
plot(freq,squeeze(G(1,1,:)))
hold on
for idx = 2:numLines
    plot(freq,squeeze(G(1,idx,:)))
end
hold off
grid on
xlabel('Frequency(Hz)')
ylabel('Value(S/m)')
txt = cell(1,numLines);
for idx = 1:numLines
    txt{1,idx} = sprintf('G1%u',idx);
end
legend(txt,'Location','best','NumColumns',2)
legend('boxoff')
title('G matrix')
% C
subplot(224)
plot(freq,squeeze(C(1,1,:)))
hold on
for idx = 2:numLines
    plot(freq,squeeze(C(1,idx,:)))
end
hold off
grid on
xlabel('Frequency(Hz)')
ylabel('Value(F/m)')
txt = cell(1,numLines);
for idx = 1:numLines
    txt{1,idx} = sprintf('C1%u',idx);
end
legend(txt,'Location','best','NumColumns',2)
legend('boxoff')
title('C matrix')

end % end for function
