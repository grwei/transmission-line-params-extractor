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
filename_4line = 'data/4line/4lines_ADS/data/4lines_HFSSW_200mil.s8p';
SingleEnded8PortData = read(rfdata.data,filename_4line);
freq = SingleEnded8PortData.Freq;
freqPts = length(freq);
z0 = SingleEnded8PortData.Z0; % Reference Impedance
SingleEnded8PortData.S_Parameters = snp2smp(SingleEnded8PortData.S_Parameters,...
    z0,1:1:8); % Classic style
numOfLines = size(SingleEnded8PortData.S_Parameters,1)/2;

%% initialize
% 此节在每个调试周期只需执行一次

% %%% Import Cadence-PowerSI-extracted params
% % Allocate memory
% rlgc_PowerSI.R = zeros(numOfLines,numOfLines,freqPts);
% rlgc_PowerSI.L = rlgc_PowerSI.R;
% rlgc_PowerSI.C = rlgc_PowerSI.R;
% rlgc_PowerSI.G = rlgc_PowerSI.R;
% % Load data
% filename_PowerSI = 'data/4line/4lines_ADS/data/Transmission_RLGC_HFSSW.csv';
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
% filename_HFSSW = 'data/4line/4lines_HFSS/m4lines_HFSS_W.csv';
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
% save('test2_3','rlgc_PowerSI','rlgc_HFSSW');

%%
load('test2_3','rlgc_PowerSI','rlgc_HFSSW');

%% Extract RLGC params using proposed method

rlgc_t = s2rlgc_t(SingleEnded8PortData.S_Parameters,lineLength,freq,z0,[],false);

%% Extracted RLGC compared with HFSS

% R
figure('Name','R matrix (compared with HFSS)')
sgtitle({'Comparison Between Proposed Algorithm and';' HFSS: R Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_t.R(1,idx,:)),'k-')
    hold on
    plot(freq/1e9,squeeze(rlgc_HFSSW.R(1,idx,:)),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('R(1,%u) (Ohms/m)',idx));
    title(sprintf('R(1,%u)',idx));
    legend({'Proposed','HFSS'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% L
figure('Name','L matrix (compared with HFSS)')
sgtitle({'Comparison Between Proposed Algorithm and';'HFSS: L Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_t.L(1,idx,:)),'k-')
    hold on
    plot(freq/1e9,squeeze(rlgc_HFSSW.L(1,idx,:)),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('L(1,%u) (H/m)',idx));
    title(sprintf('L(1,%u)',idx));
    legend({'Proposed','HFSS'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% G
figure('Name','G matrix (compared with HFSS)')
sgtitle({'Comparison Between Proposed Algorithm and';'HFSS: G Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_t.G(1,idx,:)),'k-')
    hold on
    plot(freq/1e9,squeeze(rlgc_HFSSW.G(1,idx,:)),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('G(1,%u) (S/m)',idx));
    title(sprintf('G(1,%u)',idx));
    legend({'Proposed','HFSS'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% C
figure('Name','C matrix (compared with HFSS)')
sgtitle({'Comparison Between Proposed Algorithm and';'HFSS: C Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_t.C(1,idx,:)),'k-')
    hold on
    plot(freq/1e9,squeeze(rlgc_HFSSW.C(1,idx,:)),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('C(1,%u) (F/m)',idx));
    title(sprintf('C(1,%u)',idx));
    legend({'Proposed','HFSS'},'Location','best','NumColumns',1)
    legend('boxoff')
end
