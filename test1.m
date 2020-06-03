%% Basic Information
%%% Overview
% Transmission-line parameters extractor
% MATLAB implementation of Patent US8892414B1
% Author Name: Guorui Wei
% Created in: 2020-03-15 12:45
% Example: Edge-Coupled Surface Microstrip 

clc; clear; close all;

%% Import data

% Import simulated data
lineLength = 0.0127; % Line Length(meters)
filename_2line = 'data/2line/2lines_Polar_500mil.s4p';
SingleEnded4PortData = read(rfdata.data,filename_2line);
freq = SingleEnded4PortData.Freq;
freqPts = length(freq);
z0 = SingleEnded4PortData.Z0; % Reference Impedance
numOfLines = size(SingleEnded4PortData.S_Parameters,1)/2;
SingleEnded4PortData.S_Parameters = snp2smp(SingleEnded4PortData.S_Parameters,...
    z0,1:1:2*numOfLines); % Classic style

% Import Cadence-PowerSI-extracted params
% Allocate memory
rlgc_PowerSI.R = zeros(numOfLines,numOfLines,freqPts);
rlgc_PowerSI.L = rlgc_PowerSI.R;
rlgc_PowerSI.C = rlgc_PowerSI.R;
rlgc_PowerSI.G = rlgc_PowerSI.R;
% Load data
filename_PowerSI = 'data/2line/Transmission_RLGC_res.csv';
opts = detectImportOptions(filename_PowerSI);
rlgc_PowerSI_mat = readtable(filename_PowerSI);
for freqIdx = 1:freqPts
    for i = 1:numOfLines
        for j = i:numOfLines
            rlgc_PowerSI.R(i,j,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,(2*numOfLines+2-i)*(i-1)/2+j-i+3}/lineLength;
            rlgc_PowerSI.L(i,j,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,(2*numOfLines+2-i)*(i-1)/2+j-i+3}/lineLength;
            rlgc_PowerSI.G(i,j,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,(2*numOfLines+2-i)*(i-1)/2+j-i+3}/lineLength;
            rlgc_PowerSI.C(i,j,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,(2*numOfLines+2-i)*(i-1)/2+j-i+3}/lineLength;
        end
    end
    % RLGC是对称阵
    for i = 1:numOfLines
        for j = i+1:numOfLines
            rlgc_PowerSI.R(j,i,freqIdx) = rlgc_PowerSI.R(i,j,freqIdx);
            rlgc_PowerSI.L(j,i,freqIdx) = rlgc_PowerSI.L(i,j,freqIdx);
            rlgc_PowerSI.G(j,i,freqIdx) = rlgc_PowerSI.G(i,j,freqIdx);
            rlgc_PowerSI.C(j,i,freqIdx) = rlgc_PowerSI.C(i,j,freqIdx);
        end
    end
end

%% Extract RLGC params using proposed method

rlgc_t = s2rlgc_t(SingleEnded4PortData.S_Parameters,lineLength,freq,z0);
check_consistence(rlgc_t.R, rlgc_t.L, rlgc_t.G, rlgc_t.C, lineLength, freq, z0);

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
% ave_y = mean(squeeze(rlgc_t.C(1,1,:)));
% del_r = 4e-5;
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
