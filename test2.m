%% Basic Information
% Overview
% Transmission-line parameters extractor
% MATLAB implementation of Patent US8892414B1
% Author Name: Guorui Wei
% Created in: 2020-03-15 12:45

clc; clear; close all;

%% Import data

% Import simulated data
lineLength = 0.02; % Line Length(meters)
filename_4line = 'data/4line_4linetoRLCG_10201603.s8p';
SingleEnded8PortData = read(rfdata.data,filename_4line);
freq = SingleEnded8PortData.Freq;
freqPts = length(freq);
z0 = SingleEnded8PortData.Z0; % Reference Impedance
SingleEnded8PortData.S_Parameters = snp2smp(SingleEnded8PortData.S_Parameters,...
    z0,[4 2 1 3 8 6 5 7]); % Classic style
numOfLines = size(SingleEnded8PortData.S_Parameters,1)/2;

% Import Cadence-PowerSI-extracted params
% allocate memory
rlgc_PowerSI.R = zeros(numOfLines,numOfLines);
rlgc_PowerSI.L = rlgc_PowerSI.R;
rlgc_PowerSI.C = rlgc_PowerSI.R;
rlgc_PowerSI.G = rlgc_PowerSI.R;
% load data
filename_PowerSI = 'data/4line_4linetoRLCG_10201603.csv';
opts = detectImportOptions(filename_PowerSI);
rlgc_PowerSI_mat = readtable(filename_PowerSI);
for freqIdx = 1:freqPts
    rlgc_PowerSI.R(3,3,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,3}/lineLength;
    rlgc_PowerSI.R(2,3,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,4}/lineLength;
    rlgc_PowerSI.R(3,4,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,5}/lineLength;
    rlgc_PowerSI.R(1,3,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,6}/lineLength;
    rlgc_PowerSI.R(2,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,7}/lineLength;
    rlgc_PowerSI.R(2,4,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,8}/lineLength;
    rlgc_PowerSI.R(1,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,9}/lineLength;
    rlgc_PowerSI.R(4,4,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,10}/lineLength;
    rlgc_PowerSI.R(1,4,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,11}/lineLength;
    rlgc_PowerSI.R(1,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,12}/lineLength;
    rlgc_PowerSI.R(2,1,freqIdx) = rlgc_PowerSI.R(1,2,freqIdx);
    rlgc_PowerSI.R(3,1,freqIdx) = rlgc_PowerSI.R(1,3,freqIdx);
    rlgc_PowerSI.R(3,2,freqIdx) = rlgc_PowerSI.R(2,3,freqIdx);
    rlgc_PowerSI.R(4,1,freqIdx) = rlgc_PowerSI.R(1,4,freqIdx);
    rlgc_PowerSI.R(4,2,freqIdx) = rlgc_PowerSI.R(2,4,freqIdx);
    rlgc_PowerSI.R(4,3,freqIdx) = rlgc_PowerSI.R(3,4,freqIdx);
    
    rlgc_PowerSI.L(3,3,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,3}/lineLength;
    rlgc_PowerSI.L(2,3,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,4}/lineLength;
    rlgc_PowerSI.L(3,4,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,5}/lineLength;
    rlgc_PowerSI.L(1,3,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,6}/lineLength;
    rlgc_PowerSI.L(2,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,7}/lineLength;
    rlgc_PowerSI.L(2,4,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,8}/lineLength;
    rlgc_PowerSI.L(1,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,9}/lineLength;
    rlgc_PowerSI.L(4,4,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,10}/lineLength;
    rlgc_PowerSI.L(1,4,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,11}/lineLength;
    rlgc_PowerSI.L(1,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,12}/lineLength;
    rlgc_PowerSI.L(2,1,freqIdx) = rlgc_PowerSI.L(1,2,freqIdx);
    rlgc_PowerSI.L(3,1,freqIdx) = rlgc_PowerSI.L(1,3,freqIdx);
    rlgc_PowerSI.L(3,2,freqIdx) = rlgc_PowerSI.L(2,3,freqIdx);
    rlgc_PowerSI.L(4,1,freqIdx) = rlgc_PowerSI.L(1,4,freqIdx);
    rlgc_PowerSI.L(4,2,freqIdx) = rlgc_PowerSI.L(2,4,freqIdx);
    rlgc_PowerSI.L(4,3,freqIdx) = rlgc_PowerSI.L(3,4,freqIdx);
    
    rlgc_PowerSI.G(3,3,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,3}/lineLength;
    rlgc_PowerSI.G(2,3,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,4}/lineLength;
    rlgc_PowerSI.G(3,4,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,5}/lineLength;
    rlgc_PowerSI.G(1,3,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,6}/lineLength;
    rlgc_PowerSI.G(2,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,7}/lineLength;
    rlgc_PowerSI.G(2,4,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,8}/lineLength;
    rlgc_PowerSI.G(1,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,9}/lineLength;
    rlgc_PowerSI.G(4,4,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,10}/lineLength;
    rlgc_PowerSI.G(1,4,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,11}/lineLength;
    rlgc_PowerSI.G(1,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,12}/lineLength;
    rlgc_PowerSI.G(2,1,freqIdx) = rlgc_PowerSI.G(1,2,freqIdx);
    rlgc_PowerSI.G(3,1,freqIdx) = rlgc_PowerSI.G(1,3,freqIdx);
    rlgc_PowerSI.G(3,2,freqIdx) = rlgc_PowerSI.G(2,3,freqIdx);
    rlgc_PowerSI.G(4,1,freqIdx) = rlgc_PowerSI.G(1,4,freqIdx);
    rlgc_PowerSI.G(4,2,freqIdx) = rlgc_PowerSI.G(2,4,freqIdx);
    rlgc_PowerSI.G(4,3,freqIdx) = rlgc_PowerSI.G(3,4,freqIdx);
    
    rlgc_PowerSI.C(3,3,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,3}/lineLength;
    rlgc_PowerSI.C(2,3,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,4}/lineLength;
    rlgc_PowerSI.C(3,4,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,5}/lineLength;
    rlgc_PowerSI.C(1,3,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,6}/lineLength;
    rlgc_PowerSI.C(2,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,7}/lineLength;
    rlgc_PowerSI.C(2,4,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,8}/lineLength;
    rlgc_PowerSI.C(1,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,9}/lineLength;
    rlgc_PowerSI.C(4,4,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,10}/lineLength;
    rlgc_PowerSI.C(1,4,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,11}/lineLength;
    rlgc_PowerSI.C(1,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,12}/lineLength;
    rlgc_PowerSI.C(2,1,freqIdx) = rlgc_PowerSI.C(1,2,freqIdx);
    rlgc_PowerSI.C(3,1,freqIdx) = rlgc_PowerSI.C(1,3,freqIdx);
    rlgc_PowerSI.C(3,2,freqIdx) = rlgc_PowerSI.C(2,3,freqIdx);
    rlgc_PowerSI.C(4,1,freqIdx) = rlgc_PowerSI.C(1,4,freqIdx);
    rlgc_PowerSI.C(4,2,freqIdx) = rlgc_PowerSI.C(2,4,freqIdx);
    rlgc_PowerSI.C(4,3,freqIdx) = rlgc_PowerSI.C(3,4,freqIdx);    
end

%% Extract RLGC params using proposed method

rlgc_t = s2rlgc_t(SingleEnded8PortData.S_Parameters,lineLength,freq,z0);

%% 3. Extract RLGC params using proposed method
figure('Name','Extracted RLGC using proposed method')
sgtitle('Extracted RLGC using proposed algorithm')
subplot(2,2,1)
plot(freq,squeeze(rlgc_t.R(1,1,:)))
grid on
hold on
plot(freq,squeeze(rlgc_t.R(2,1,:)))
plot(freq,squeeze(rlgc_t.R(2,2,:)))
plot(freq,squeeze(rlgc_t.R(3,1,:)))
plot(freq,squeeze(rlgc_t.R(3,2,:)))
plot(freq,squeeze(rlgc_t.R(4,1,:)))
hold off
xlabel('Frequency')
xlim([0 3e9])
ylabel('R(Ohms/m)')
legend({'R11','R21','R22','R31','R32','R41'},'Location','best','NumColumns',2)
legend('boxoff')
title('R matrix')

subplot(2,2,2)
plot(freq,squeeze(rlgc_t.L(1,1,:)))
grid on
hold on
plot(freq,squeeze(rlgc_t.L(2,1,:)))
plot(freq,squeeze(rlgc_t.L(2,2,:)))
plot(freq,squeeze(rlgc_t.L(3,1,:)))
plot(freq,squeeze(rlgc_t.L(3,2,:)))
plot(freq,squeeze(rlgc_t.L(4,1,:)))
hold off
xlabel('Frequency')
xlim([0 3e9])
ylabel('L(H/m)')
legend({'L11','L21','L22','L31','L32','L41'},'Location','best','NumColumns',2)
legend('boxoff')
title('L matrix')

subplot(2,2,3)
plot(freq,squeeze(rlgc_t.C(1,1,:)))
grid on
hold on
plot(freq,squeeze(rlgc_t.C(2,1,:)))
plot(freq,squeeze(rlgc_t.C(2,2,:)))
plot(freq,squeeze(rlgc_t.C(3,1,:)))
plot(freq,squeeze(rlgc_t.C(3,2,:)))
plot(freq,squeeze(rlgc_t.C(4,1,:)))
hold off
xlabel('Frequency')
xlim([0 3e9])
ylabel('C(F/m)')
legend({'C11','C21','C22','C31','C32','C41'},'Location','best','NumColumns',2)
legend('boxoff')
title('C matrix')

subplot(2,2,4)
plot(freq,squeeze(rlgc_t.G(1,1,:)))
grid on
hold on
plot(freq,squeeze(rlgc_t.G(2,1,:)))
plot(freq,squeeze(rlgc_t.G(2,2,:)))
plot(freq,squeeze(rlgc_t.G(3,1,:)))
plot(freq,squeeze(rlgc_t.G(3,2,:)))
plot(freq,squeeze(rlgc_t.G(4,1,:)))
hold off
xlabel('Frequency')
xlim([0 3e9])
ylabel('G(S/m)')
legend({'G11','G21','G22','G31','G32','G41'},'Location','best','NumColumns',2)
legend('boxoff')
title('G matrix')

%% Candence Sigrity PowerSI extracted RLGC
figure('Name','Candence Sigrity PowerSI extracted RLGC')
sgtitle('Candence Sigrity PowerSI extracted RLGC')
subplot(2,2,1)
plot(freq,squeeze(rlgc_PowerSI.R(1,1,:)))
grid on
hold on
plot(freq,squeeze(rlgc_PowerSI.R(2,1,:)))
plot(freq,squeeze(rlgc_PowerSI.R(2,2,:)))
plot(freq,squeeze(rlgc_PowerSI.R(3,1,:)))
plot(freq,squeeze(rlgc_PowerSI.R(3,2,:)))
plot(freq,squeeze(rlgc_PowerSI.R(4,1,:)))
hold off
xlabel('Frequency')
xlim([0 3e9])
ylabel('R(Ohms/m)')
legend({'R11','R21','R22','R31','R32','R41'},'Location','best','NumColumns',2)
legend('boxoff')
title('R matrix')

subplot(2,2,2)
plot(freq,squeeze(rlgc_PowerSI.L(1,1,:)))
grid on
hold on
plot(freq,squeeze(rlgc_PowerSI.L(2,1,:)))
plot(freq,squeeze(rlgc_PowerSI.L(2,2,:)))
plot(freq,squeeze(rlgc_PowerSI.L(3,1,:)))
plot(freq,squeeze(rlgc_PowerSI.L(3,2,:)))
plot(freq,squeeze(rlgc_PowerSI.L(4,1,:)))
hold off
xlabel('Frequency')
xlim([0 3e9])
ylabel('L(H/m)')
legend({'L11','L21','L22','L31','L32','L41'},'Location','best','NumColumns',2)
legend('boxoff')
title('L matrix')

subplot(2,2,3)
plot(freq,squeeze(rlgc_PowerSI.C(1,1,:)))
grid on
hold on
plot(freq,squeeze(rlgc_PowerSI.C(2,1,:)))
plot(freq,squeeze(rlgc_PowerSI.C(2,2,:)))
plot(freq,squeeze(rlgc_PowerSI.C(3,1,:)))
plot(freq,squeeze(rlgc_PowerSI.C(3,2,:)))
plot(freq,squeeze(rlgc_PowerSI.C(4,1,:)))
hold off
xlabel('Frequency')
xlim([0 3e9])
ylabel('C(F/m)')
legend({'C11','C21','C22','C31','C32','C41'},'Location','best','NumColumns',2)
legend('boxoff')
title('C matrix')

subplot(2,2,4)
plot(freq,squeeze(rlgc_PowerSI.G(1,1,:)))
grid on
hold on
plot(freq,squeeze(rlgc_PowerSI.G(2,1,:)))
plot(freq,squeeze(rlgc_PowerSI.G(2,2,:)))
plot(freq,squeeze(rlgc_PowerSI.G(3,1,:)))
plot(freq,squeeze(rlgc_PowerSI.G(3,2,:)))
plot(freq,squeeze(rlgc_PowerSI.G(4,1,:)))
hold off
xlabel('Frequency')
xlim([0 3e9])
ylabel('G(S/m)')
legend({'G11','G21','G22','G31','G32','G41'},'Location','best','NumColumns',2)
legend('boxoff')
title('G matrix')

%% Extracted RLGC compared with Cadence PowerSI

% R, L
figure('Name','R, L matrix')
sgtitle({'Comparison Between Proposed Algorithm and';' Cadence Sigrity PowerSI: R, L Matrix'})
subplot(221)
plot(freq/1e9,squeeze(rlgc_t.R(1,1,:)),'k-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.R(1,1,:)),'g--')
hold off
grid on
xlabel('Freq(GHz)');
ylabel('R11(Ohms/m)');
title('R11 Comparison');
legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

subplot(222)
plot(freq/1e9,squeeze(rlgc_t.R(1,2,:)),'k-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.R(1,2,:)),'g--')
hold off
grid on
xlabel('Freq(GHz)');
ylabel('R12(Ohms/m)');
title('R12 Comparison');
legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

subplot(223)
plot(freq/1e9,squeeze(rlgc_t.L(1,1,:)),'k-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.L(1,1,:)),'g--')
hold off
grid on
xlabel('Freq(GHz)');
ylabel('L11(H/m)');
title('L11 Comparison');
legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

subplot(224)
plot(freq/1e9,squeeze(rlgc_t.L(1,2,:)),'k-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.L(1,2,:)),'g--')
hold off
grid on
xlabel('Freq(GHz)');
ylabel('L12(H/m)');
title('L12 Comparison');
legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

% C, G
figure('Name','C, G matrix')
sgtitle({'Comparison Between Proposed Algorithm and';' Cadence Sigrity PowerSI: C, G Matrix'})
subplot(221)
plot(freq/1e9,squeeze(rlgc_t.C(1,1,:)),'m-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.C(1,1,:)),'k--')
hold off
grid on
% ave_y = mean(squeeze(rlgc_t.C(1,1,:)));
% del_r = 4e-3;
% ylim([ave_y*(1-del_r), ave_y*(1+del_r)])
xlabel('Freq(GHz)');
ylabel('C11(F/m)');
title('C11 Comparison');
legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

subplot(222)
plot(freq/1e9,squeeze(rlgc_t.C(1,2,:)),'m-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.C(1,2,:)),'k--')
hold off
grid on
% ave_y = mean(squeeze(rlgc_t.C(1,2,:)));
% del_r = 5e-4;
% ylim([ave_y*(1+del_r), ave_y*(1-del_r)])
xlabel('Freq(GHz)');
ylabel('C12(F/m)');
title('C12 Comparison');
legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

subplot(223)
plot(freq/1e9,squeeze(rlgc_t.G(1,1,:)),'m-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.G(1,1,:)),'k--')
hold off
grid on
xlabel('Freq(GHz)');
ylabel('G11(S/m)');
title('G11 Comparison');
legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

subplot(224)
plot(freq/1e9,squeeze(rlgc_t.G(1,2,:)),'m-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.G(1,2,:)),'k--')
hold off
grid on
xlabel('Freq(GHz)');
ylabel('G12(S/m)');
title('G12 Comparison');
legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')
