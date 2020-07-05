%% Basic Information
%%% Overview
% Transmission-line parameters extractor
% MATLAB implementation of Patent US8892414B1
% Author Name: Guorui Wei
% Created in: 2020-05-27 12:45

clc; clear; close all;

%% initialize
% 此节只需在每次数据改变后执行一次。
% 读取仿真数据，然后以.mat文件存储到工程根目录，以缩短程序多次运行时读数据时间。

% %%% Import simulated S
% filename_16line = 'data/16line/16lines_HFSS/16lines_HFSS_200mil.s32p';
% SingleEnded32PortData = read(rfdata.data,filename_16line);
% numOfLines = size(SingleEnded32PortData.S_Parameters,1)/2;
% freq = SingleEnded32PortData.Freq;
% freqPts = length(freq);
% lineLength = 0.00508; % Line Length(meters)
%  
% %%% Import Cadence-PowerSI-extracted params
% % Allocate memory
% rlgc_PowerSI.R = zeros(numOfLines,numOfLines,freqPts);
% rlgc_PowerSI.L = rlgc_PowerSI.R;
% rlgc_PowerSI.C = rlgc_PowerSI.R;
% rlgc_PowerSI.G = rlgc_PowerSI.R;
% % Load data
% filename_PowerSI = 'data/16line/16lines_HFSS/Transmission_RLGC_res.csv';
% opts = detectImportOptions(filename_PowerSI);
% rlgc_PowerSI_mat = readtable(filename_PowerSI);
% for freqIdx = 1:freqPts
%     for i = 1:numOfLines
%         for j = i:numOfLines
%             rlgc_PowerSI.R(i,j,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,(2*numOfLines+2-i)*(i-1)/2+j-i+3}/lineLength;
%             rlgc_PowerSI.L(i,j,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,(2*numOfLines+2-i)*(i-1)/2+j-i+3}/lineLength;
%             rlgc_PowerSI.G(i,j,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,(2*numOfLines+2-i)*(i-1)/2+j-i+3}/lineLength;
%             rlgc_PowerSI.C(i,j,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,(2*numOfLines+2-i)*(i-1)/2+j-i+3}/lineLength;
%         end
%     end
%     % RLGC是对称阵
%     for i = 1:numOfLines
%         for j = i+1:numOfLines
%             rlgc_PowerSI.R(j,i,freqIdx) = rlgc_PowerSI.R(i,j,freqIdx);
%             rlgc_PowerSI.L(j,i,freqIdx) = rlgc_PowerSI.L(i,j,freqIdx);
%             rlgc_PowerSI.G(j,i,freqIdx) = rlgc_PowerSI.G(i,j,freqIdx);
%             rlgc_PowerSI.C(j,i,freqIdx) = rlgc_PowerSI.C(i,j,freqIdx);
%         end
%     end
% end
% 
% %%% Read W-element (HFSS) file
% % Allocate memory
% rlgc_HFSSW.R = zeros(numOfLines,numOfLines,freqPts);
% rlgc_HFSSW.L = rlgc_HFSSW.R;
% rlgc_HFSSW.C = rlgc_HFSSW.R;
% rlgc_HFSSW.G = rlgc_HFSSW.R;
% % load data
% filename_HFSSW = 'data/16line/16lines_HFSS/m16lines_HFSS_W.csv';
% HFSSW_mat = readtable(filename_HFSSW);
% for freqIdx = 1:freqPts
%     for i = 1:numOfLines
%         for j = 1:i
%             rlgc_HFSSW.R(i,j,freqIdx) = HFSSW_mat{(freqIdx-1)*(1+(1+numOfLines)*numOfLines/2)+i*(i-1)/2+j+1,1};
%             rlgc_HFSSW.L(i,j,freqIdx) = HFSSW_mat{(freqIdx-1)*(1+(1+numOfLines)*numOfLines/2)+i*(i-1)/2+j+1,2};
%             rlgc_HFSSW.G(i,j,freqIdx) = HFSSW_mat{(freqIdx-1)*(1+(1+numOfLines)*numOfLines/2)+i*(i-1)/2+j+1,3};
%             rlgc_HFSSW.C(i,j,freqIdx) = HFSSW_mat{(freqIdx-1)*(1+(1+numOfLines)*numOfLines/2)+i*(i-1)/2+j+1,4};
%         end
%     end
%     
%     % 对称
%     for i = 1:numOfLines
%         for j = 1:i-1
%             rlgc_HFSSW.R(j,i,freqIdx) = rlgc_HFSSW.R(i,j,freqIdx);
%             rlgc_HFSSW.L(j,i,freqIdx) = rlgc_HFSSW.L(i,j,freqIdx);
%             rlgc_HFSSW.G(j,i,freqIdx) = rlgc_HFSSW.G(i,j,freqIdx);
%             rlgc_HFSSW.C(j,i,freqIdx) = rlgc_HFSSW.C(i,j,freqIdx);
%         end
%     end
% end
% 
% %%% Save some data
% save('test3','SingleEnded32PortData','rlgc_PowerSI','rlgc_HFSSW');

%% process data

load('test3','SingleEnded32PortData','rlgc_PowerSI','rlgc_HFSSW');
% process simulated data
lineLength = 0.00508; % Line Length(meters)
freq = SingleEnded32PortData.Freq;
freqPts = length(freq);
z0 = SingleEnded32PortData.Z0; % Reference Impedance
SingleEnded32PortData.S_Parameters = snp2smp(SingleEnded32PortData.S_Parameters,...
    z0,1:1:32); % Classic style
numOfLines = size(SingleEnded32PortData.S_Parameters,1)/2;

%% Extract RLGC params using proposed method

rlgc_t = s2rlgc_t(SingleEnded32PortData.S_Parameters,lineLength,freq,z0,[],false);
% check_consistence(rlgc_t.R, rlgc_t.L, rlgc_t.G, rlgc_t.C, lineLength, freq, z0);

%% Extracted RLGC compared with Cadence PowerSI, HFSS

% R
figure('Name','R (compared with PowerSI, HFSS)')
sgtitle({'Comparison Between Proposed Algorithm and';' PowerSI and HFSS: R Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_t.R(1,idx,:)),'k-')
    hold on
    plot(freq/1e9,squeeze(rlgc_PowerSI.R(1,idx,:)),'g--')
    plot(freq/1e9,squeeze(rlgc_HFSSW.R(1,idx,:)),'m--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('R(1,%u) (Ohms/m)',idx));
    title(sprintf('R(1,%u)',idx));
    legend({'Proposed','PowerSI','HFSS'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% L
figure('Name','L (compared with PowerSI, HFSS)')
sgtitle({'Comparison Between Proposed Algorithm and';' Cadence Sigrity PowerSI: L Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_t.L(1,idx,:)),'k-')
    hold on
    plot(freq/1e9,squeeze(rlgc_PowerSI.L(1,idx,:)),'g--')
    plot(freq/1e9,squeeze(rlgc_HFSSW.L(1,idx,:)),'m--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('L(1,%u) (H/m)',idx));
    title(sprintf('L(1,%u)',idx));
    legend({'Proposed','PowerSI','HFSS'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% G
figure('Name','G (compared with PowerSI, HFSS)')
sgtitle({'Comparison Between Proposed Algorithm and';' PowerSI and HFSS: G Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_t.G(1,idx,:)),'k-')
    hold on
    plot(freq/1e9,squeeze(rlgc_PowerSI.G(1,idx,:)),'g--')
    plot(freq/1e9,squeeze(rlgc_HFSSW.G(1,idx,:)),'m--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('G(1,%u) (S/m)',idx));
    title(sprintf('G(1,%u)',idx));
    legend({'Proposed','PowerSI','HFSS'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% C
figure('Name','C (compared with PowerSI, HFSS)')
sgtitle({'Comparison Between Proposed Algorithm and';' PowerSI and HFSS: C Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_t.C(1,idx,:)),'k-')
    hold on
    plot(freq/1e9,squeeze(rlgc_PowerSI.C(1,idx,:)),'g--')
    plot(freq/1e9,squeeze(rlgc_HFSSW.C(1,idx,:)),'m--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('C(1,%u) (F/m)',idx));
    title(sprintf('C(1,%u)',idx));
    legend({'Proposed','PowerSI','HFSS'},'Location','best','NumColumns',1)
    legend('boxoff')
end

%% Calculated S-parameters (dB) using powerSI-extracted RLGC
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
