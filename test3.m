%% Basic Information
%%% Overview
% Transmission-line parameters extractor
% MATLAB implementation of Patent US8892414B1
% Author Name: Guorui Wei
% Created in: 2020-05-27 12:45

clc; clear; close all;

%% Import data

% Import simulated data
lineLength = 0.00508; % Line Length(meters)
filename_16line = 'data/16line/16-lines_8-diff-pairs.s32p';
SingleEnded32PortData = read(rfdata.data,filename_16line);
freq = SingleEnded32PortData.Freq;
freqPts = length(freq);
z0 = SingleEnded32PortData.Z0; % Reference Impedance
port_order = nan(1,32);
port_order(2:2:16) = 1:1:8;
port_order(1:2:15) = 9:1:16;
port_order(17:1:32) = port_order(1:1:16) + 16;
SingleEnded32PortData.S_Parameters = snp2smp(SingleEnded32PortData.S_Parameters,...
    z0,port_order); % Classic style
numOfLines = size(SingleEnded32PortData.S_Parameters,1)/2;

% Import Cadence-PowerSI-extracted params
% Allocate memory
rlgc_PowerSI.R = zeros(numOfLines,numOfLines,freqPts);
rlgc_PowerSI.L = rlgc_PowerSI.R;
rlgc_PowerSI.C = rlgc_PowerSI.R;
rlgc_PowerSI.G = rlgc_PowerSI.R;
% Load data
filename_PowerSI = 'data/16line/RLGC_PowerSI.csv';
opts = detectImportOptions(filename_PowerSI);
rlgc_PowerSI_mat = readtable(filename_PowerSI);
for freqIdx = 1:freqPts
    for i = 1:numOfLines
        for j = i:numOfLines
            rlgc_PowerSI.R(i,j,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,(34-i)*(i-1)/2+j-i+3}/lineLength;
            rlgc_PowerSI.L(i,j,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,(34-i)*(i-1)/2+j-i+3}/lineLength;
            rlgc_PowerSI.G(i,j,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,(34-i)*(i-1)/2+j-i+3}/lineLength;
            rlgc_PowerSI.C(i,j,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,(34-i)*(i-1)/2+j-i+3}/lineLength;
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
    % 重排端口顺序
    rlgc_PowerSI.R([2:2:16,1:2:15],:,freqIdx) = rlgc_PowerSI.R(:,:,freqIdx);
    rlgc_PowerSI.R(:,[2:2:16,1:2:15],freqIdx) = rlgc_PowerSI.R(:,:,freqIdx);
    rlgc_PowerSI.L([2:2:16,1:2:15],:,freqIdx) = rlgc_PowerSI.L(:,:,freqIdx);
    rlgc_PowerSI.L(:,[2:2:16,1:2:15],freqIdx) = rlgc_PowerSI.L(:,:,freqIdx);
    rlgc_PowerSI.G([2:2:16,1:2:15],:,freqIdx) = rlgc_PowerSI.G(:,:,freqIdx);
    rlgc_PowerSI.G(:,[2:2:16,1:2:15],freqIdx) = rlgc_PowerSI.G(:,:,freqIdx);
    rlgc_PowerSI.C([2:2:16,1:2:15],:,freqIdx) = rlgc_PowerSI.C(:,:,freqIdx);
    rlgc_PowerSI.C(:,[2:2:16,1:2:15],freqIdx) = rlgc_PowerSI.C(:,:,freqIdx);
end

%% Extract RLGC params using proposed method

rlgc_t = s2rlgc_t(SingleEnded32PortData.S_Parameters,lineLength,freq,z0,[],true);
check_consistence(rlgc_t.R, rlgc_t.L, rlgc_t.G, rlgc_t.C, lineLength, freq, z0);

%% Extracted RLGC compared with Cadence PowerSI

% R
figure('Name','R matrix (compared with PowerSI)')
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
    ylabel(sprintf('R(1,%u)(Ohms/m)',idx));
    title(sprintf('R(1,%u)',idx));
    legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% L
figure('Name','L matrix (compared with PowerSI)')
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
    ylabel(sprintf('L(1,%u)(H/m)',idx));
    title(sprintf('L(1,%u)',idx));
    legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% G
figure('Name','G matrix (compared with PowerSI)')
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
    ylabel(sprintf('G(1,%u)(S/m)',idx));
    title(sprintf('G(1,%u)',idx));
    legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% C
figure('Name','C matrix (compared with PowerSI)')
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
    ylabel(sprintf('C(1,%u)(F/m)',idx));
    title(sprintf('C(1,%u)',idx));
    legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
    legend('boxoff')
end

%% Calculated S-parameters using powerSI-extracted RLGC
% 结论：与原始S参数高度一致。用重建S参数再提取RLGC，此RLGC可以准确恢复重建S参数。
% 说明PowerSI是对原始S参数作了前处理

[s_params_using_cadence_extracted_rlgc,~] = rlgc2s_t(rlgc_PowerSI.R,rlgc_PowerSI.L,rlgc_PowerSI.G,rlgc_PowerSI.C,lineLength,freq,z0);
% external<-external
figure('Name','Rebuilt S (using PowerSI-RLGC): See') 
sgtitle({'Comparison Between Calculated S-parameters using';'PowerSI-RLGC and Original S-parameters: See'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,db(squeeze(s_params_using_cadence_extracted_rlgc(1,idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(SingleEnded32PortData.S_Parameters(1,idx,:)),'voltage'),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(1,%u)(Ohms/m)',idx));
    title(sprintf('S(1,%u)',idx));
    legend({'PowerSI-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% external<-internal
figure('Name','Rebuilt S (using PowerSI-RLGC): Sei')
sgtitle({'Comparison Between Calculated S-parameters using';'PowerSI-RLGC and Original S-parameters: Sei'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,db(squeeze(s_params_using_cadence_extracted_rlgc(1,idx+numOfLines,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(SingleEnded32PortData.S_Parameters(1,idx+numOfLines,:)), 'voltage'),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S1%u(dB)',idx+numOfLines));
    title(sprintf('S1%u',idx+numOfLines));
    legend({'PowerSI-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% internal<-external
figure('Name','Rebuilt S (using PowerSI-RLGC): Sie') 
sgtitle({'Comparison Between Calculated S-parameters using';'PowerSI-RLGC and Original S-parameters: Sie'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,db(squeeze(s_params_using_cadence_extracted_rlgc(numOfLines+1,idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(SingleEnded32PortData.S_Parameters(numOfLines+1,idx,:)),'voltage'),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S%u%u(Ohms/m)',numOfLines+1,idx));
    title(sprintf('S%u%u',numOfLines+1,idx));
    legend({'PowerSI-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% internal<-internal
figure('Name','Rebuilt S (using PowerSI-RLGC): Sii') 
sgtitle({'Comparison Between Calculated S-parameters using';'PowerSI-RLGC and Original S-parameters: Sii'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,db(squeeze(s_params_using_cadence_extracted_rlgc(numOfLines+1,numOfLines+idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(SingleEnded32PortData.S_Parameters(numOfLines+1,numOfLines+idx,:)),'voltage'),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S%u%u(Ohms/m)',numOfLines+1,numOfLines+idx));
    title(sprintf('S%u%u',numOfLines+1,numOfLines+idx));
    legend({'PowerSI-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
    legend('boxoff')
end

%% 初始S参数经PowerSI提取RLGC，用此RLGC重建S以实现对S参数的“合理化”，再用本文方法提取RLGC
% 测试结论：完全一致！说明PowerSI是对S参数作了前处理！

rlgc_from_refined_S = s2rlgc_t(s_params_using_cadence_extracted_rlgc,lineLength,freq,z0,[],false);
%%% Extracted RLGC compared with Cadence PowerSI

% R
figure('Name','R matrix (using refined S)')
sgtitle({'Comparison Between Proposed Algorithm and';' Cadence Sigrity PowerSI: R Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_from_refined_S.R(1,idx,:)),'k-')
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
figure('Name','L matrix (using refined S)')
sgtitle({'Comparison Between Proposed Algorithm and';' Cadence Sigrity PowerSI: L Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_from_refined_S.L(1,idx,:)),'k-')
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
figure('Name','G matrix (using refined S)')
sgtitle({'Comparison Between Proposed Algorithm and';' Cadence Sigrity PowerSI: G Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_from_refined_S.G(1,idx,:)),'k-')
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
figure('Name','C matrix (using refined S)')
sgtitle({'Comparison Between Proposed Algorithm and';' Cadence Sigrity PowerSI: C Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_from_refined_S.C(1,idx,:)),'k-')
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

%% Rebuilt S-parameters (phase) using PowerSI-RLGC
% Expected to be consistent with the original S-parameters
% 
% external<-external

figure('Name','Rebuilt S (phase) (PowerSI-RLGC): See') 
sgtitle({'Comparison Between Rebuilt S-parameters and';'Original S-parameters: See'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,angle(squeeze(s_params_using_cadence_extracted_rlgc(1,idx,:))),'k-')
    hold on
    plot(freq/1e9,angle(squeeze(SingleEnded32PortData.S_Parameters(1,idx,:))),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(1,%u)(dB)',idx));
    title(sprintf('S(1,%u)',idx));
    legend({'PowerSI-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% external<-internal
figure('Name','Rebuilt S (phase) (PowerSI-RLGC): Sei')
sgtitle({'Comparison Between Rebuilt S-parameters and';'Original S-parameters: Sei'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,angle(squeeze(s_params_using_cadence_extracted_rlgc(1,idx+numOfLines,:))),'k-')
    hold on
    plot(freq/1e9,angle(squeeze(SingleEnded32PortData.S_Parameters(1,idx+numOfLines,:))),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(1,%u)(dB)',idx+numOfLines));
    title(sprintf('S(1,%u)',idx+numOfLines));
    legend({'PowerSI-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% internal<-external
figure('Name','Rebuilt S (phase) (PowerSI-RLGC): Sie') 
sgtitle({'Comparison Between Rebuilt S-parameters and';'Original S-parameters: Sie'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,angle(squeeze(s_params_using_cadence_extracted_rlgc(numOfLines+1,idx,:))),'k-')
    hold on
    plot(freq/1e9,angle(squeeze(SingleEnded32PortData.S_Parameters(numOfLines+1,idx,:))),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(%u,%u)(dB)',numOfLines+1,idx));
    title(sprintf('S(%u,%u)',numOfLines+1,idx));
    legend({'PowerSI-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% internal<-internal
figure('Name','Rebuilt S (phase) (PowerSI-RLGC): Sii') 
sgtitle({'Comparison Between Rebuilt S-parameters and';'Original S-parameters: Sii'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,angle(squeeze(s_params_using_cadence_extracted_rlgc(numOfLines+1,numOfLines+idx,:))),'k-')
    hold on
    plot(freq/1e9,angle(squeeze(SingleEnded32PortData.S_Parameters(numOfLines+1,numOfLines+idx,:))),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(%u,%u)(dB)',numOfLines+1,numOfLines+idx));
    title(sprintf('S(%u,%u)',numOfLines+1,numOfLines+idx));
    legend({'PowerSI-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
    legend('boxoff')
end