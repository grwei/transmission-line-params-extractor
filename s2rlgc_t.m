function output = s2rlgc_t(s_params,linelength,freq,z0,port_reorder,debug_mode)
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
%
%   See also S2RLGC, ABCD2S, S2Y, S2Z, S2H, Y2ABCD, Z2ABCD, H2ABCD, RLGC2S

%% Input validity check and initialization

narginchk(3,6)

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

%% Configuration

debugFlag = true;      % Debug flag: 0-released; 1-debug
if nargin > 5 && debug_mode == false
    debugFlag = false;
end
EigSortMethod = 0;      % 0(default)- abs(real(...))  1- abs(...)

%% Convert S-parameters to ABCD-Parameters [Sampath2008]

Z_params    = zeros(2*numLines,2*numLines,freqpts); % impedance matrix
TA          = zeros(numLines,numLines,freqpts);     % the A term of transmission(ABCD) matrix:[A B;C D]
TB          = TA;
% The C and D terms are not actually used.
TC          = TA;
TD          = TA;
I           = eye(2*numLines,2*numLines);           % 2Nx2N identity matrix
z0_matrix   = z0 * I;                               % reference characteristic impedance matrix

for idx=1:freqpts
    Z_params(:,:,idx) = z0_matrix*(I+s_params(:,:,idx))         /       ...
                        (I-s_params(:,:,idx));
    %%% test [Reveyrand2018]: S->Z
    % It seems that both formula works well.
    % should be removed when released
%     Z_params(:,:,idx) =    (I-s_params(:,:,idx))                \       ...
%                            (I+s_params(:,:,idx))*z0_matrix;
    TA(:,:,idx) =   Z_params(1:numLines,1:numLines,idx)         /       ...
                    Z_params(numLines+1:end,1:numLines,idx);
    TB(:,:,idx) =   TA(:,:,idx)                                 *       ...
                    Z_params(numLines+1:end,numLines+1:end,idx) -       ...
                    Z_params(1:numLines,numLines+1:end,idx);
    % The C and D terms are not actually used.
    TC(:,:,idx) =   eye(numLines,numLines)                      /       ...
                    Z_params(numLines+1:end,1:numLines,idx);
    TD(:,:,idx) =   TC(:,:,idx)                                 *       ...
                    Z_params(numLines+1:end,numLines+1:end,idx);
end

%% Extract Complex Propagation Constants while preserving relative Eigenvalue position

%%% Eigenvalue decomposition(similarity transformation) of TA
eigVal = zeros(numLines,freqpts);           % Eigenvalues of TA
eigVec = zeros(numLines,numLines,freqpts);  % Columns are the corresponding Right Eigenvectors of TA
for idx = 1:freqpts
    % [V,D] = eig(TA,'vector') --> TA*V = V*diag(D)
    % Note: The eigenvectors should be normalized[Braunisch1998]
    % so that the 2-norm of each is 1, which is done by default
    % by the MATLAB function eig().
    [eigVec(:,:,idx),eigVal(:,idx)] = eig(TA(:,:,idx),'vector');
end

%%% Adjust the order of Eigenvalues and Eigenvectors [Braunisch1998, Chu2015]
prodTable   = nan(numLines,numLines,freqpts);   % Hermitian Inner Product recorder.
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
    
    %%% Determine the correct position by Hermitian Inner Product.
    % WARNNING! In extreme cases, duplicate serial numbers may appear, 
    % which can lead to fatal errors.
    
    % !!!Two choices: 
    % 0. using the absolute value of the real part, i.e. max(abs(real(...)));   
    % 1. using the modulus of complex number, i.e. max(abs(...))
    switch(EigSortMethod)
        case 0
            [~,CorrectPos(:,freqidx)] = max(abs(real(prodTable(:,:,freqidx))),[],2);
        case 1
            [~,CorrectPos(:,freqidx)] = max(abs(prodTable(:,:,freqidx)),[],2);
    end
    %%% Reorder the Eigenvectors and the corresponding Eigenvalues
    [~,newIndex(:,freqidx)] = sort(CorrectPos(:,freqidx));
    %%% test: 不排序
    % 结论：不影响S重建，但会导致RLGC非物理
%     newIndex(:,freqidx) = 1:numLines;
    
    eigVec(:,:,freqidx)     = eigVec(:,newIndex(:,freqidx),freqidx);
    eigVal(:,freqidx)       = eigVal(newIndex(:,freqidx),freqidx);
end

%% Extract Attenuation Constants and Unwrapped Phase Constants

gammaLenEigWrap = acosh(eigVal);  % Principle Value of gammaEig*linelength:
%%% test: arccosh(z) = log(z +/- sqrt(z^2 - 1));
% Very bad
% gammaLenEigWrap_1 = log(eigVal + sqrt(eigVal.^2 - 1));
% gammaLenEigWrap_2 = log(eigVal - sqrt(eigVal.^2 - 1));
% alphaLenEigWrap = abs(real(gammaLenEigWrap_1 - gammaLenEigWrap_2))/2;
% betaLenEigWrap = abs(imag(gammaLenEigWrap_1 - gammaLenEigWrap_2))/2;
% betaLenEigWrap = cumsum([betaLenEigWrap(:,1),abs(diff(betaLenEigWrap,1,2))],2); % continious phase[0,+inf)
% gammaLenEigWrap = complex(alphaLenEigWrap,betaLenEigWrap);

% Real part should be non-negative, imag part in (-pi,pi]
betaLenEigWrapDiff(:,2:freqpts) = diff(imag(gammaLenEigWrap),1,2);
discontCount = cumsum(abs(betaLenEigWrapDiff) > pi,2);

%%% !!!Phase-unwrapping algorithm is unreliable near singular frequency! --grwei,20200318
% Maybe position-tracking algorithm instead???
% discontCount(:,491:end) = 1;

%%% discontinuity-detection-based phase unwrapping

%%% Test: 不作解折叠
% 结论：求出的RLGC非物理，但不影响重建S参数
% discontCount(:) = 0;

betaLenEigUnwrap = imag(gammaLenEigWrap) + 2*pi*discontCount;

%%% Test: 首次解折叠后，再去除多余的上跳点 (Good -grwei, 20200604)
% 见[US8892414B1]第10-14行
betaLenEigUnwrap = abs(betaLenEigUnwrap);   % beta >= 0 
discontCount(1) = 1;                        % 初始化
while sum(discontCount(:))
    betaLenEigUnWrapDiff(:,2:freqpts) = diff(betaLenEigUnwrap,1,2);
    discontCount = cumsum(betaLenEigUnWrapDiff > pi,2);
    betaLenEigUnwrap = betaLenEigUnwrap - 2*pi*discontCount;
end

%%% Test: 再去除多余的下跳点 (Good -grwei, 20200604)
% 见[US8892414B1]第10-14行
discontCount(1) = 1;            % 初始化
while sum(discontCount(:))
    betaLenEigUnWrapDiff(:,2:freqpts) = diff(betaLenEigUnwrap,1,2);
    discontCount = cumsum(betaLenEigUnWrapDiff < -pi,2);
    betaLenEigUnwrap = betaLenEigUnwrap + 2*pi*discontCount;
end

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

%%% Test: 改用解折叠后的gamma计算Zc
% 见[US8892414B1]第10-43行
for idx = 1:freqpts
    Zc(:,:,idx) = eigVec(:,:,idx)                           *           ...
        diag(1./sinh(gammaEigUnwrap(:,idx) * linelength))   /           ...
        eigVec(:,:,idx)                                     *           ...
        TB(:,:,idx);
    
    %%% Test: 强制令Zc的实部非负
%     for i = 1:numLines
%         for j = 1:numLines
%             if real(Zc(i,j,idx)) < 0
%                 Zc(i,j,idx) = -Zc(i,j,idx);
%             end
%         end
%     end
    
end

%% Extract RLGC

%%% Test: RLGC -> ABCD,S不是单射！
% 结论：gammaEig的j2pi周期会导致提取的RLGC非物理，但不影响S重建
% for idx = 1:freqpts
%     gammaEigUnwrap(1:2,idx)  =  gammaEigUnwrap(1:2,idx) + 200*pi*1i/linelength;
%     gammaEigUnwrap(3:4,idx)  =  gammaEigUnwrap(3:4,idx) + 100*pi*1i/linelength;
%     gamma(:,:,idx) = eigVec(:,:,idx)        *                           ...
%         diag(gammaEigUnwrap(:,idx))         /                           ...
%         eigVec(:,:,idx);
% end

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

%% debug
%The following code snippets are for debugging purposes only.

if debugFlag == false
    return
end

%% Hermitian inner product table
% the basis of Eigenvalues(and their corresponding Eigenvectors)
% reordering, and discontinuity-detection-based phase-unwrapping algorithm

figure('Name','prodTable (Real-part)')
sgtitle('Hermitian Inner Product (Real-part)')
sz = 10; % Marker area
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
    ylabel('Abs (Real-part)')
    %     ylim([0 1.1])
end

%% Hermitian inner product table
% the basis of Eigenvalues(and their corresponding Eigenvectors)
% reordering, and discontinuity-detection-based phase-unwrapping algorithm

figure('Name','prodTable (Magnitude)')
sgtitle('Hermitian Inner Product (Magnitude)')
sz = 10; % Marker area
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

%% Propagation constant(before unwrapping) of each eigen-mode

figure('Name','Propagation constant (before unwrapping) of each eigen-mode')
sgtitle('Gamma (wrap) of each eigen mode')
sz = 10; % Marker area
subplot(121)
% plot(freq,real(gammaLenEigWrap(1,:))/linelength)
scatter(freq,real(gammaLenEigWrap(1,:))/linelength,sz,'o')
hold on
for idx = 2:numLines
%     plot(freq,real(gammaLenEigWrap(idx,:))/linelength)
    scatter(freq,real(gammaLenEigWrap(idx,:))/linelength,sz,'o');
end
hold off
grid on
% xlim([4.5e9 5.5e9])
xlabel('Frequency(Hz)')
ylabel('\alpha(Np/m)')
txt = cell(1,numLines);
for idx = 1:numLines
    txt{1,idx} = ['\alpha_{',sprintf('%u}',idx)];
end
legend(txt,'Location','best','NumColumns',2)
legend('boxoff')
title('\alpha')

subplot(122)
% plot(freq,imag(gammaLenEigWrap(1,:)))
scatter(freq,imag(gammaLenEigWrap(1,:)),sz,'o')
hold on
for idx = 2:numLines
%     plot(freq,imag(gammaLenEigWrap(idx,:)))
    scatter(freq,imag(gammaLenEigWrap(idx,:)),sz,'o')
end
hold off
grid on
% xlim([4.5e9 5.5e9])
xlabel('Frequency(Hz)')
ylabel('\betaL(rad)')
txt = cell(1,numLines);
for idx = 1:numLines
    txt{1,idx} = ['\beta_{',sprintf('%u}',idx)];
end
legend(txt,'Location','best','NumColumns',2)
legend('boxoff')
title('\betaL')

%% Extracted characteristic impedance matrix

figure('Name','Zc (Real part)')
sgtitle('Charateristic Impedance Matrix (Real part)')
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
        txt{1,idx_col} = ['Zc_{',sprintf('%u,',idx),sprintf('%u',idx_col),'}'];
    end
    legend(txt,'Location','best','NumColumns',2)
    legend('boxoff')
    xlabel('Frequency (Hz)')
    ylabel('Real-part (Ohms)')
    title(['Zc_{',sprintf('%u',idx),'X}'])
end

%
figure('Name','Zc (Imag part)')
sgtitle('Charateristic Impedance Matrix (Imag part)')
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
        txt{1,idx_col} = ['Zc_{',sprintf('%u,',idx),sprintf('%u',idx_col),'}'];
    end
    legend(txt,'Location','best','NumColumns',2)
    legend('boxoff')
    xlabel('Frequency(Hz)')
    ylabel('Imag-part(Ohms)')
    title(['Zc_{',sprintf('%u',idx),'X}'])
end

%% Propagation constant of each eigen-mode

figure('Name','gammaEigUnwrap')
sgtitle('Gamma (unwrap) of each eigen mode')
subplot(121)
% plot(freq,real(gammaEigUnwrap(1,:)))
scatter(freq,real(gammaEigUnwrap(1,:)),sz,'o')
hold on
for idx = 2:numLines
%     plot(freq,real(gammaEigUnwrap(idx,:)))
    scatter(freq,real(gammaEigUnwrap(idx,:)),sz,'o')
end
hold off
grid on
% xlim([4.5e9 5.5e9])
xlabel('Frequency(Hz)')
ylabel('\alpha(Np/m)')
txt = cell(1,numLines);
for idx = 1:numLines
    txt{1,idx} = ['\alpha_{',sprintf('%u}',idx)];
end
legend(txt,'Location','best','NumColumns',2)
legend('boxoff')
title('\alpha')

subplot(122)
% plot(freq,imag(gammaEigUnwrap(1,:)))
scatter(freq,imag(gammaEigUnwrap(1,:)),sz,'o')
hold on
for idx = 2:numLines
%     plot(freq,imag(gammaEigUnwrap(idx,:)))
    scatter(freq,imag(gammaEigUnwrap(idx,:)),sz,'o')
end
hold off
grid on
% xlim([4.5e9 5.5e9])
xlabel('Frequency(Hz)')
ylabel('\beta(rad/m)')
txt = cell(1,numLines);
for idx = 1:numLines
    txt{1,idx} = ['\beta_{',sprintf('%u}',idx)];
end
legend(txt,'Location','best','NumColumns',2)
legend('boxoff')
title('\beta')

%% Extracted RLGC (Overview)

figure('Name','RLGC matrix (Overview)')
sgtitle({'Extracted RLGC matrix (Overview)'})
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
    txt{1,idx} = sprintf('R(1,%u)',idx);
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
    txt{1,idx} = sprintf('L(1,%u)',idx);
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
    txt{1,idx} = sprintf('G(1,%u)',idx);
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
    txt{1,idx} = sprintf('C(1,%u)',idx);
end
legend(txt,'Location','best','NumColumns',2)
legend('boxoff')
title('C matrix')

%% Rebuilt S-parameters (dB) using extracted RLGC
% Expected to be consistent with the original S-parameters
% 
% <<doc\pic\port-ordering.png>>
% 

% Rebuild S-parameters using extracted RLGC
[s_params_rebuilt,~] = rlgc2s_t(R,L,G,C,linelength,freq,z0);

% external<-external
figure('Name','Rebuilt S (Extracted-RLGC): See') 
sgtitle({'Comparison Between Rebuilt S-parameters and';'Original S-parameters: See'})
num_of_columes = ceil(numLines/2);
for idx = 1:numLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,db(squeeze(s_params_rebuilt(1,idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(s_params(1,idx,:)),'voltage'),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(1,%u)(dB)',idx));
    title(sprintf('S(1,%u)',idx));
    legend({'Extracted-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% external<-internal
figure('Name','Rebuilt S (Extracted-RLGC): Sei')
sgtitle({'Comparison Between Rebuilt S-parameters and';'Original S-parameters: Sei'})
num_of_columes = ceil(numLines/2);
for idx = 1:numLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,db(squeeze(s_params_rebuilt(1,idx+numLines,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(s_params(1,idx+numLines,:)), 'voltage'),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(1,%u)(dB)',idx+numLines));
    title(sprintf('S(1,%u)',idx+numLines));
    legend({'Extracted-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% internal<-external
figure('Name','Rebuilt S (Extracted-RLGC): Sie') 
sgtitle({'Comparison Between Rebuilt S-parameters and';'Original S-parameters: Sie'})
num_of_columes = ceil(numLines/2);
for idx = 1:numLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,db(squeeze(s_params_rebuilt(numLines+1,idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(s_params(numLines+1,idx,:)),'voltage'),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(%u,%u)(dB)',numLines+1,idx));
    title(sprintf('S(%u,%u)',numLines+1,idx));
    legend({'Extracted-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% internal<-internal
figure('Name','Rebuilt S (Extracted-RLGC): Sii') 
sgtitle({'Comparison Between Rebuilt S-parameters and';'Original S-parameters: Sii'})
num_of_columes = ceil(numLines/2);
for idx = 1:numLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,db(squeeze(s_params_rebuilt(numLines+1,numLines+idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(s_params(numLines+1,numLines+idx,:)),'voltage'),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(%u,%u)(dB)',numLines+1,numLines+idx));
    title(sprintf('S(%u,%u)',numLines+1,numLines+idx));
    legend({'Extracted-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
    legend('boxoff')
end

%% Rebuilt S-parameters (phase) using extracted RLGC
% Expected to be consistent with the original S-parameters
% 
% <<doc\pic\port-ordering.png>>
% 

% external<-external
figure('Name','Rebuilt S (phase) (Extracted-RLGC): See') 
sgtitle({'Comparison Between Rebuilt S-parameters and';'Original S-parameters: See'})
num_of_columes = ceil(numLines/2);
for idx = 1:numLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,angle(squeeze(s_params_rebuilt(1,idx,:))),'k-')
    hold on
    plot(freq/1e9,angle(squeeze(s_params(1,idx,:))),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(1,%u) (rad)',idx));
    title(sprintf('S(1,%u)',idx));
    legend({'Extracted-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% external<-internal
figure('Name','Rebuilt S (phase) (Extracted-RLGC): Sei')
sgtitle({'Comparison Between Rebuilt S-parameters and';'Original S-parameters: Sei'})
num_of_columes = ceil(numLines/2);
for idx = 1:numLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,angle(squeeze(s_params_rebuilt(1,idx+numLines,:))),'k-')
    hold on
    plot(freq/1e9,angle(squeeze(s_params(1,idx+numLines,:))),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(1,%u) (rad)',idx+numLines));
    title(sprintf('S(1,%u)',idx+numLines));
    legend({'Extracted-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% internal<-external
figure('Name','Rebuilt S (phase) (Extracted-RLGC): Sie') 
sgtitle({'Comparison Between Rebuilt S-parameters and';'Original S-parameters: Sie'})
num_of_columes = ceil(numLines/2);
for idx = 1:numLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,angle(squeeze(s_params_rebuilt(numLines+1,idx,:))),'k-')
    hold on
    plot(freq/1e9,angle(squeeze(s_params(numLines+1,idx,:))),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(%u,%u) (rad)',numLines+1,idx));
    title(sprintf('S(%u,%u)',numLines+1,idx));
    legend({'Extracted-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% internal<-internal
figure('Name','Rebuilt S (phase) (Extracted-RLGC): Sii') 
sgtitle({'Comparison Between Rebuilt S-parameters and';'Original S-parameters: Sii'})
num_of_columes = ceil(numLines/2);
for idx = 1:numLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,angle(squeeze(s_params_rebuilt(numLines+1,numLines+idx,:))),'k-')
    hold on
    plot(freq/1e9,angle(squeeze(s_params(numLines+1,numLines+idx,:))),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(%u,%u) (rad)',numLines+1,numLines+idx));
    title(sprintf('S(%u,%u)',numLines+1,numLines+idx));
    legend({'Extracted-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
    legend('boxoff')
end

%% Rebuilt S-parameters (phase) using extracted RLGC
% Expected to be consistent with the original S-parameters
% 
% <<doc\pic\port-ordering.png>>
% 

% external<-external
figure('Name','Rebuilt S (phase) (Extracted-RLGC): See') 
sgtitle({'Comparison Between Rebuilt S-parameters and';'Original S-parameters: See'})
num_of_columes = ceil(numLines/2);
for idx = 1:numLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,angle(squeeze(s_params_rebuilt(1,idx,:))),'k-')
    hold on
    plot(freq/1e9,angle(squeeze(s_params(1,idx,:))),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(1,%u) (rad)',idx));
    title(sprintf('S(1,%u)',idx));
    legend({'Extracted-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% external<-internal
figure('Name','Rebuilt S (phase) (Extracted-RLGC): Sei')
sgtitle({'Comparison Between Rebuilt S-parameters and';'Original S-parameters: Sei'})
num_of_columes = ceil(numLines/2);
for idx = 1:numLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,angle(squeeze(s_params_rebuilt(1,idx+numLines,:))),'k-')
    hold on
    plot(freq/1e9,angle(squeeze(s_params(1,idx+numLines,:))),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(1,%u) (rad)',idx+numLines));
    title(sprintf('S(1,%u)',idx+numLines));
    legend({'Extracted-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% internal<-external
figure('Name','Rebuilt S (phase) (Extracted-RLGC): Sie') 
sgtitle({'Comparison Between Rebuilt S-parameters and';'Original S-parameters: Sie'})
num_of_columes = ceil(numLines/2);
for idx = 1:numLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,angle(squeeze(s_params_rebuilt(numLines+1,idx,:))),'k-')
    hold on
    plot(freq/1e9,angle(squeeze(s_params(numLines+1,idx,:))),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(%u,%u) (rad)',numLines+1,idx));
    title(sprintf('S(%u,%u)',numLines+1,idx));
    legend({'Extracted-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% internal<-internal
figure('Name','Rebuilt S (phase) (Extracted-RLGC): Sii') 
sgtitle({'Comparison Between Rebuilt S-parameters and';'Original S-parameters: Sii'})
num_of_columes = ceil(numLines/2);
for idx = 1:numLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,angle(squeeze(s_params_rebuilt(numLines+1,numLines+idx,:))),'k-')
    hold on
    plot(freq/1e9,angle(squeeze(s_params(numLines+1,numLines+idx,:))),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(%u,%u) (rad)',numLines+1,numLines+idx));
    title(sprintf('S(%u,%u)',numLines+1,numLines+idx));
    legend({'Extracted-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
    legend('boxoff')
end

end % end for function
