%% Basic Information
%%% Overview
% Transmission-line parameters extractor
% MATLAB implementation of Patent US8892414B1
% Author Name: Guorui Wei
% Created in: 2020-03-15 12:45

clc; clear; close all;

%% Import data

% Import simulated data
lineLength = 0.02; % Line Length(meters)
filename_4line = 'data/Four-line_20mm_20191020.s8p';
SingleEnded8PortData = read(rfdata.data,filename_4line);
freq = SingleEnded8PortData.Freq;
freqPts = length(freq);
z0 = SingleEnded8PortData.Z0; % Reference Impedance
SingleEnded8PortData.S_Parameters = snp2smp(SingleEnded8PortData.S_Parameters,...
    z0,[4 2 1 3 8 6 5 7]); % Classic style
numOfLines = size(SingleEnded8PortData.S_Parameters,1)/2;

% Import Cadence-PowerSI-extracted params
% Allocate memory
rlgc_PowerSI.R = zeros(numOfLines,numOfLines,freqPts);
rlgc_PowerSI.L = rlgc_PowerSI.R;
rlgc_PowerSI.C = rlgc_PowerSI.R;
rlgc_PowerSI.G = rlgc_PowerSI.R;
% Load data
filename_PowerSI = 'data/Four-line_20mm_20191020_PowerSI.csv';
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
% [s_params_rebuilt, rlgc_rebuilt] = rlgc2s_t(rlgc_t.R,rlgc_t.L,rlgc_t.G,rlgc_t.C,lineLength,freq,z0);

%% Extracted RLGC using proposed method(before 1st resonance frequency)

figure('Name','Extracted RLGC using proposed method(before 1st resonance frequency)')
sgtitle({'Extracted RLGC using proposed algorithm','(before 1st resonance frequency)'})
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

%% Candence Sigrity PowerSI extracted RLGC(before 1st resonance frequency)

figure('Name','Candence Sigrity PowerSI extracted RLGC(before 1st resonance frequency)')
sgtitle({'Candence Sigrity PowerSI extracted RLGC','(before 1st resonance frequency)'})
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

% R
figure('Name','R matrix')
sgtitle({'Comparison Between Proposed Algorithm and';' Cadence Sigrity PowerSI: R Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_t.R(1,idx,:)),'k-')
    hold on
    plot(freq/1e9,squeeze(rlgc_PowerSI.R(1,idx,:)),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('R1%u(Ohms/m)',idx));
    title(sprintf('R1%u',idx));
    legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% L
figure('Name','L matrix')
sgtitle({'Comparison Between Proposed Algorithm and';' Cadence Sigrity PowerSI: L Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_t.L(1,idx,:)),'k-')
    hold on
    plot(freq/1e9,squeeze(rlgc_PowerSI.L(1,idx,:)),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('L1%u(H/m)',idx));
    title(sprintf('L1%u',idx));
    legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% G
figure('Name','G matrix')
sgtitle({'Comparison Between Proposed Algorithm and';' Cadence Sigrity PowerSI: G Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_t.G(1,idx,:)),'k-')
    hold on
    plot(freq/1e9,squeeze(rlgc_PowerSI.G(1,idx,:)),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('G1%u(S/m)',idx));
    title(sprintf('G1%u',idx));
    legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% C
figure('Name','C matrix')
sgtitle({'Comparison Between Proposed Algorithm and';' Cadence Sigrity PowerSI: C Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_t.C(1,idx,:)),'k-')
    hold on
    plot(freq/1e9,squeeze(rlgc_PowerSI.C(1,idx,:)),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('C1%u(F/m)',idx));
    title(sprintf('C1%u',idx));
    legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
    legend('boxoff')
end
