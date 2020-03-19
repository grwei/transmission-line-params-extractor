%% Basic Information

% Overview
% Transmission-line parameters extractor
% MATLAB implementation of Patent US8892414B1
% Author Name: Guorui Wei
% Created in: 2020-03-15 12:45

clc; clear; close all;

%% Import data

lineLength = 0.0254; % Line Length(meters)

% Import simulated data
filename_2line = 'data/CPL_1in_20200202.s4p';
SingleEnded4PortData = read(rfdata.data,filename_2line);
freq = SingleEnded4PortData.Freq;
freqPts = length(freq);
z0 = SingleEnded4PortData.Z0; % Reference Impedance
SingleEnded4PortData.S_Parameters = snp2smp(SingleEnded4PortData.S_Parameters,...
    z0,[1 2 3 4]); % Classic style
numOfLines = size(SingleEnded4PortData.S_Parameters,1)/2;

% Import Cadence-PowerSI-extracted params
filename_PowerSI = 'data/CPL_1in_20200202_PowerSI.csv';
opts = detectImportOptions(filename_PowerSI);
rlgc_PowerSI_mat = readtable(filename_PowerSI);
% Allocate memory
rlgc_PowerSI.R = zeros(numOfLines,numOfLines,freqPts);
rlgc_PowerSI.L = rlgc_PowerSI.R;
rlgc_PowerSI.C = rlgc_PowerSI.R;
rlgc_PowerSI.G = rlgc_PowerSI.R;
% Load data
for freqIdx = 1:freqPts
    rlgc_PowerSI.R(1,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,3}/lineLength;
    rlgc_PowerSI.R(1,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,4}/lineLength;
    rlgc_PowerSI.R(2,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,4}/lineLength;
    rlgc_PowerSI.R(2,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,5}/lineLength;
    rlgc_PowerSI.L(1,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,3}/lineLength;
    rlgc_PowerSI.L(1,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,4}/lineLength;
    rlgc_PowerSI.L(2,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,4}/lineLength;
    rlgc_PowerSI.L(2,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,5}/lineLength;    
    rlgc_PowerSI.G(1,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,3}/lineLength;
    rlgc_PowerSI.G(1,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,4}/lineLength;
    rlgc_PowerSI.G(2,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,4}/lineLength;
    rlgc_PowerSI.G(2,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,5}/lineLength;
    rlgc_PowerSI.C(1,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,3}/lineLength;
    rlgc_PowerSI.C(1,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,4}/lineLength;
    rlgc_PowerSI.C(2,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,4}/lineLength;
    rlgc_PowerSI.C(2,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,5}/lineLength;
end

%% Extract RLGC params using proposed method

rlgc_t = s2rlgc_t(SingleEnded4PortData.S_Parameters,lineLength,freq,z0);

%% Extracted RLGC compared with Cadence Sigrity PowerSI

% R, L
figure('Name','R, L matrix')
sgtitle({'Comparison Between Proposed Algorithm and';'Cadence Sigrity PowerSI: R, L Matrix'})
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
ave_y = mean(squeeze(rlgc_t.C(1,1,:)));
del_r = 4e-5;
ylim([ave_y*(1-del_r), ave_y*(1+del_r)])
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
ave_y = mean(squeeze(rlgc_t.C(1,2,:)));
del_r = 5e-4;
ylim([ave_y*(1+del_r), ave_y*(1-del_r)])
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
