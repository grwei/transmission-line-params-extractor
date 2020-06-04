%% Basic Information
%%% Overview
% Transmission-line parameters extractor
% MATLAB implementation of Patent US8892414B1
% Author Name: Guorui Wei
% Created in: 2020-03-15 12:45

clc; clear; close all;

%% Import data

% Import simulated data
lineLength = 0.00508; % Line Length(meters)
filename_4line = 'data/4line/4lines_HFSS/4lines_HFSS_200mil.s8p';
SingleEnded8PortData = read(rfdata.data,filename_4line);
freq = SingleEnded8PortData.Freq;
freqPts = length(freq);
z0 = SingleEnded8PortData.Z0; % Reference Impedance
SingleEnded8PortData.S_Parameters = snp2smp(SingleEnded8PortData.S_Parameters,...
    z0,1:1:8); % Classic style
numOfLines = size(SingleEnded8PortData.S_Parameters,1)/2;

% Import Cadence-PowerSI-extracted params
% Allocate memory
rlgc_PowerSI.R = zeros(numOfLines,numOfLines,freqPts);
rlgc_PowerSI.L = rlgc_PowerSI.R;
rlgc_PowerSI.C = rlgc_PowerSI.R;
rlgc_PowerSI.G = rlgc_PowerSI.R;
% Load data
filename_PowerSI = 'data/4line/4lines_HFSS/Transmission_RLGC_res.csv';
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

rlgc_t = s2rlgc_t(SingleEnded8PortData.S_Parameters,lineLength,freq,z0,[],true);
check_consistence(rlgc_t.R, rlgc_t.L, rlgc_t.G, rlgc_t.C, lineLength, freq, z0);
%%% Test: rational fit
% 结论：无区别？
% data = SingleEnded8PortData.S_Parameters; 
% fit_data = rationalfit(freq, data);
% [resp,freq]=freqresp(fit_data,freq);
% rlgc_t = s2rlgc_t(resp,lineLength,freq,z0,[],true);

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
    ylabel(sprintf('R1%u(Ohms/m)',idx));
    title(sprintf('R1%u',idx));
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
    ylabel(sprintf('L1%u(H/m)',idx));
    title(sprintf('L1%u',idx));
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
    ylabel(sprintf('G1%u(S/m)',idx));
    title(sprintf('G1%u',idx));
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
    ylabel(sprintf('C1%u(F/m)',idx));
    title(sprintf('C1%u',idx));
    legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
    legend('boxoff')
end

%% Calculated S-parameters using powerSI-extracted RLGC
% 结论：与原始S参数高度一致。用重建S参数再提取RLGC，此RLGC可以准确恢复重建S参数。
% 说明PowerSI是对原始S参数作了前处理

[s_params_using_cadence_extracted_rlgc,~] = rlgc2s_t(rlgc_PowerSI.R,rlgc_PowerSI.L,rlgc_PowerSI.G,rlgc_PowerSI.C,lineLength,freq,z0);
% % external<-external
% figure('Name','Rebuilt S (using PowerSI-RLGC): See') 
% sgtitle({'Comparison Between Calculated S-parameters using';'PowerSI-RLGC and Original S-parameters: See'})
% num_of_columes = ceil(numOfLines/2);
% for idx = 1:numOfLines
%     subplot(2,num_of_columes,idx)
%     plot(freq/1e9,db(squeeze(s_params_using_cadence_extracted_rlgc(1,idx,:)),'voltage'),'k-')
%     hold on
%     plot(freq/1e9,db(squeeze(SingleEnded8PortData.S_Parameters(1,idx,:)),'voltage'),'g--')
%     hold off
%     grid on
%     xlabel('Freq(GHz)');
%     ylabel(sprintf('S1%u(Ohms/m)',idx));
%     title(sprintf('S1%u',idx));
%     legend({'PowerSI-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
%     legend('boxoff')
% end
% 
% % external<-internal
% figure('Name','Rebuilt S (using PowerSI-RLGC): Sei')
% sgtitle({'Comparison Between Calculated S-parameters using';'PowerSI-RLGC and Original S-parameters: Sei'})
% num_of_columes = ceil(numOfLines/2);
% for idx = 1:numOfLines
%     subplot(2,num_of_columes,idx)
%     plot(freq/1e9,db(squeeze(s_params_using_cadence_extracted_rlgc(1,idx+numOfLines,:)),'voltage'),'k-')
%     hold on
%     plot(freq/1e9,db(squeeze(SingleEnded8PortData.S_Parameters(1,idx+numOfLines,:)), 'voltage'),'g--')
%     hold off
%     grid on
%     xlabel('Freq(GHz)');
%     ylabel(sprintf('S1%u(dB)',idx+numOfLines));
%     title(sprintf('S1%u',idx+numOfLines));
%     legend({'PowerSI-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
%     legend('boxoff')
% end
% 
% % internal<-external
% figure('Name','Rebuilt S (using PowerSI-RLGC): Sie') 
% sgtitle({'Comparison Between Calculated S-parameters using';'PowerSI-RLGC and Original S-parameters: Sie'})
% num_of_columes = ceil(numOfLines/2);
% for idx = 1:numOfLines
%     subplot(2,num_of_columes,idx)
%     plot(freq/1e9,db(squeeze(s_params_using_cadence_extracted_rlgc(numOfLines+1,idx,:)),'voltage'),'k-')
%     hold on
%     plot(freq/1e9,db(squeeze(SingleEnded8PortData.S_Parameters(numOfLines+1,idx,:)),'voltage'),'g--')
%     hold off
%     grid on
%     xlabel('Freq(GHz)');
%     ylabel(sprintf('S%u%u(Ohms/m)',numOfLines+1,idx));
%     title(sprintf('S%u%u',numOfLines+1,idx));
%     legend({'PowerSI-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
%     legend('boxoff')
% end
% 
% % internal<-internal
% figure('Name','Rebuilt S (using PowerSI-RLGC): Sii') 
% sgtitle({'Comparison Between Calculated S-parameters using';'PowerSI-RLGC and Original S-parameters: Sii'})
% num_of_columes = ceil(numOfLines/2);
% for idx = 1:numOfLines
%     subplot(2,num_of_columes,idx)
%     plot(freq/1e9,db(squeeze(s_params_using_cadence_extracted_rlgc(numOfLines+1,numOfLines+idx,:)),'voltage'),'k-')
%     hold on
%     plot(freq/1e9,db(squeeze(SingleEnded8PortData.S_Parameters(numOfLines+1,numOfLines+idx,:)),'voltage'),'g--')
%     hold off
%     grid on
%     xlabel('Freq(GHz)');
%     ylabel(sprintf('S%u%u(Ohms/m)',numOfLines+1,numOfLines+idx));
%     title(sprintf('S%u%u',numOfLines+1,numOfLines+idx));
%     legend({'PowerSI-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
%     legend('boxoff')
% end

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
