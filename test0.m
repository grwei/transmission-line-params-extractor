%% Basic Information
%%% Overview
% Transmission-line parameters extractor
% MATLAB implementation of Patent US8892414B1
% Author Name: Guorui Wei
% Created in: 2020-03-15 12:45

clc; clear; close all;

%% Import data

lineLength = 0.0254; % Line Length(meters)

% Import simulated data
filename_1line = 'data/MSL_1in_20200319.s2p';
SingleEnded2PortData = read(rfdata.data,filename_1line);
freq = SingleEnded2PortData.Freq;
freqPts = length(freq);
z0 = SingleEnded2PortData.Z0; % Reference Impedance
SingleEnded2PortData.S_Parameters = snp2smp(SingleEnded2PortData.S_Parameters,...
    z0,[1 2]); % Classic style
numOfLines = size(SingleEnded2PortData.S_Parameters,1)/2;

% Import Cadence-PowerSI-extracted params
filename_PowerSI = 'data/MSL_1in_20200319_PowerSI.csv';
opts = detectImportOptions(filename_PowerSI);
rlgc_PowerSI_mat = readtable(filename_PowerSI);
% Allocate memory
rlgc_PowerSI.R = zeros(numOfLines,numOfLines,freqPts);
rlgc_PowerSI.L = rlgc_PowerSI.R;
rlgc_PowerSI.C = rlgc_PowerSI.R;
rlgc_PowerSI.G = rlgc_PowerSI.R;
% Import data
for idx = 1:freqPts
    rlgc_PowerSI.R(:,:,idx) = rlgc_PowerSI_mat{4*idx-3,3}/lineLength;
    rlgc_PowerSI.L(:,:,idx) = rlgc_PowerSI_mat{4*idx-2,3}/lineLength;
    rlgc_PowerSI.G(:,:,idx) = rlgc_PowerSI_mat{4*idx-1,3}/lineLength;
    rlgc_PowerSI.C(:,:,idx) = rlgc_PowerSI_mat{4*idx-0,3}/lineLength;
end

%% Extract RLGC params using proposed method

rlgc_t = s2rlgc_t(SingleEnded2PortData.S_Parameters,lineLength,freq,z0);
check_consistence(rlgc_t.R, rlgc_t.L, rlgc_t.G, rlgc_t.C, lineLength, freq, z0);

%% Extracted RLGC compared with Cadence Sigrity PowerSI

figure('Name','RLGC Comparison')
sgtitle({'Comparison Between Proposed Algorithm and';'Cadence Sigrity PowerSI: Extracted RLGC'})
subplot(221)
plot(freq/1e9,squeeze(rlgc_t.R(1,1,:)),'k-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.R(1,1,:)),'g--')
hold off
grid on
xlabel('Freq(GHz)');
ylabel('R(Ohms/m)');
title('R Comparison');
legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

subplot(222)
plot(freq/1e9,squeeze(rlgc_t.L(1,1,:)),'k-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.L(1,1,:)),'g--')
hold off
grid on
xlabel('Freq(GHz)');
ylabel('L(H/m)');
title('L Comparison');
legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

subplot(223)
plot(freq/1e9,squeeze(rlgc_t.G(1,1,:)),'k-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.G(1,1,:)),'g--')
hold off
grid on
xlabel('Freq(GHz)');
ylabel('G(S/m)');
title('G Comparison');
legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

subplot(224)
plot(freq/1e9,squeeze(rlgc_t.C(1,1,:)),'k-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.C(1,1,:)),'g--')
hold off
grid on
xlabel('Freq(GHz)');
ylabel('C(F/m)');
title('C Comparison');
legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')
