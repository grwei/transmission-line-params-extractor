%% Basic Information
%%% Overview
% Transmission-line parameters extractor
% MATLAB implementation of Patent US8892414B1
% Author Name: Guorui Wei
% Created in: 2020-06-04 00:35
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
s_original = SingleEnded4PortData.S_Parameters;

%% Add noise to original S-params

% initialize the random number generator to make the results  repeatable.
rng(0,'twister');

noise_mean = 1;   % mean
noise_std = 0.05; % standard deviation
noise_cof_real = noise_std .* randn(size(s_original)) + noise_mean;
noise_cof_imag = noise_std .* randn(size(s_original)) + noise_mean;
s_noisy = real(s_original) .* noise_cof_real                            ...
        + 1i * imag(s_original) .* noise_cof_imag;
rfwrite(sparameters(s_noisy,freq,z0),'data/2line/2lines_noisy_500mil.s4p');

%%% Load PowerSI extracted noisy-RLGC
% Allocate memory
rlgc_noisy_PowerSI.R = zeros(numOfLines,numOfLines,freqPts);
rlgc_noisy_PowerSI.L = rlgc_noisy_PowerSI.R;
rlgc_noisy_PowerSI.C = rlgc_noisy_PowerSI.R;
rlgc_noisy_PowerSI.G = rlgc_noisy_PowerSI.R;
% Load data
filename_PowerSI = 'data/2line/Transmission_RLGC_noisy.csv';
opts = detectImportOptions(filename_PowerSI);
rlgc_PowerSI_mat = readtable(filename_PowerSI);
for freqIdx = 1:freqPts
    for i = 1:numOfLines
        for j = i:numOfLines
            rlgc_noisy_PowerSI.R(i,j,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,(2*numOfLines+2-i)*(i-1)/2+j-i+3}/lineLength;
            rlgc_noisy_PowerSI.L(i,j,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,(2*numOfLines+2-i)*(i-1)/2+j-i+3}/lineLength;
            rlgc_noisy_PowerSI.G(i,j,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,(2*numOfLines+2-i)*(i-1)/2+j-i+3}/lineLength;
            rlgc_noisy_PowerSI.C(i,j,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,(2*numOfLines+2-i)*(i-1)/2+j-i+3}/lineLength;
        end
    end
    % RLGC是对称阵
    for i = 1:numOfLines
        for j = i+1:numOfLines
            rlgc_noisy_PowerSI.R(j,i,freqIdx) = rlgc_noisy_PowerSI.R(i,j,freqIdx);
            rlgc_noisy_PowerSI.L(j,i,freqIdx) = rlgc_noisy_PowerSI.L(i,j,freqIdx);
            rlgc_noisy_PowerSI.G(j,i,freqIdx) = rlgc_noisy_PowerSI.G(i,j,freqIdx);
            rlgc_noisy_PowerSI.C(j,i,freqIdx) = rlgc_noisy_PowerSI.C(i,j,freqIdx);
        end
    end
end

%% Extract RLGC params using proposed method

rlgc_original = s2rlgc_t(s_original,lineLength,freq,z0,[],false);
rlgc_noisy = s2rlgc_t(s_noisy,lineLength,freq,z0,[],false);

%%% Check

% 1. rlgc_noisy -> s_extracted_rlgc -> rlgc'
check_consistence(rlgc_noisy.R, rlgc_noisy.L, rlgc_noisy.G, rlgc_noisy.C, lineLength, freq, z0);
% 2. s_extracted_rlgc -> rlgc' -> s'
[s_extracted_rlgc,~] = rlgc2s_t(rlgc_noisy.R,rlgc_noisy.L,rlgc_noisy.G,rlgc_noisy.C,lineLength,freq,z0);
rlgc_reextract = s2rlgc_t(s_extracted_rlgc,lineLength,freq,z0,[],true);

%% Figure

%% Comparison between noisy-S (dB) and original-S

% external<-external
figure('Name','noisy-S (dB): See') 
sgtitle({'Comparison Between noisy-S and original-S: See'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,db(squeeze(s_noisy(1,idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(s_original(1,idx,:)),'voltage'),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(1,%u)(dB)',idx));
    title(sprintf('S(1,%u)',idx));
    legend({'noisy','original'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% external<-internal
figure('Name','noisy-S (dB): Sei')
sgtitle({'Comparison Between noisy-S and original-S: Sei'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,db(squeeze(s_noisy(1,idx+numOfLines,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(s_original(1,idx+numOfLines,:)), 'voltage'),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(1,%u)(dB)',idx+numOfLines));
    title(sprintf('S(1,%u)',idx+numOfLines));
    legend({'noisy','original'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% internal<-external
figure('Name','noisy-S (dB): Sie') 
sgtitle({'Comparison Between noisy-S and original-S: Sie'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,db(squeeze(s_noisy(numOfLines+1,idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(s_original(numOfLines+1,idx,:)),'voltage'),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(%u,%u)(dB)',numOfLines+1,idx));
    title(sprintf('S(%u,%u)',numOfLines+1,idx));
    legend({'noisy','original'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% internal<-internal
figure('Name','noisy-S (dB): Sii') 
sgtitle({'Comparison Between noisy-S and original-S: Sii'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,db(squeeze(s_noisy(numOfLines+1,numOfLines+idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(s_original(numOfLines+1,numOfLines+idx,:)),'voltage'),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(%u,%u)(dB)',numOfLines+1,numOfLines+idx));
    title(sprintf('S(%u,%u)',numOfLines+1,numOfLines+idx));
    legend({'noisy','original'},'Location','best','NumColumns',1)
    legend('boxoff')
end

%% Comparison between noisy-S (phase) and original-S

% external<-external
figure('Name','noisy-S (phase): See') 
sgtitle({'Comparison Between noisy-S and original-S: See'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,angle(squeeze(s_noisy(1,idx,:))),'k-')
    hold on
    plot(freq/1e9,angle(squeeze(s_original(1,idx,:))),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(1,%u) (rad)',idx));
    title(sprintf('S(1,%u)',idx));
    legend({'noisy','original'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% external<-internal
figure('Name','noisy-S (phase): Sei')
sgtitle({'Comparison Between noisy-S and original-S: Sei'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,angle(squeeze(s_noisy(1,idx+numOfLines,:))),'k-')
    hold on
    plot(freq/1e9,angle(squeeze(s_original(1,idx+numOfLines,:))),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(1,%u) (rad)',idx+numOfLines));
    title(sprintf('S(1,%u)',idx+numOfLines));
    legend({'noisy','original'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% internal<-external
figure('Name','noisy-S (phase): Sie') 
sgtitle({'Comparison Between noisy-S and original-S: Sie'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,angle(squeeze(s_noisy(numOfLines+1,idx,:))),'k-')
    hold on
    plot(freq/1e9,angle(squeeze(s_original(numOfLines+1,idx,:))),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(%u,%u) (rad)',numOfLines+1,idx));
    title(sprintf('S(%u,%u)',numOfLines+1,idx));
    legend({'noisy','original'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% internal<-internal
figure('Name','noisy-S (phase): Sii') 
sgtitle({'Comparison Between noisy-S and original-S: Sii'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,angle(squeeze(s_noisy(numOfLines+1,numOfLines+idx,:))),'k-')
    hold on
    plot(freq/1e9,angle(squeeze(s_original(numOfLines+1,numOfLines+idx,:))),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S(%u,%u) (rad)',numOfLines+1,numOfLines+idx));
    title(sprintf('S(%u,%u)',numOfLines+1,numOfLines+idx));
    legend({'noisy','original'},'Location','best','NumColumns',1)
    legend('boxoff')
end

%% Reconstruct noisy S using powerSI-noisy RLGC
% 结论：与原始S参数高度一致。用重建S参数再提取RLGC，此RLGC可以准确恢复重建S参数。
% 说明PowerSI是对原始S参数作了前处理

[s_rlgc_noisy_PowerSI,~] = rlgc2s_t(rlgc_noisy_PowerSI.R,rlgc_noisy_PowerSI.L,rlgc_noisy_PowerSI.G,rlgc_noisy_PowerSI.C,lineLength,freq,z0);
% external<-external
figure('Name','Rebuilt S (using PowerSI-RLGC): See') 
sgtitle({'Comparison Between Calculated S-parameters using';'PowerSI-RLGC and Original S-parameters: See'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,db(squeeze(s_rlgc_noisy_PowerSI(1,idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(s_noisy(1,idx,:)),'voltage'),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S1%u(Ohms/m)',idx));
    title(sprintf('S1%u',idx));
    legend({'PowerSI-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% external<-internal
figure('Name','Rebuilt S (using PowerSI-RLGC): Sei')
sgtitle({'Comparison Between Calculated S-parameters using';'PowerSI-RLGC and Original S-parameters: Sei'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    plot(freq/1e9,db(squeeze(s_rlgc_noisy_PowerSI(1,idx+numOfLines,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(s_noisy(1,idx+numOfLines,:)), 'voltage'),'g--')
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
    plot(freq/1e9,db(squeeze(s_rlgc_noisy_PowerSI(numOfLines+1,idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(s_noisy(numOfLines+1,idx,:)),'voltage'),'g--')
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
    plot(freq/1e9,db(squeeze(s_rlgc_noisy_PowerSI(numOfLines+1,numOfLines+idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(s_noisy(numOfLines+1,numOfLines+idx,:)),'voltage'),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('S%u%u(Ohms/m)',numOfLines+1,numOfLines+idx));
    title(sprintf('S%u%u',numOfLines+1,numOfLines+idx));
    legend({'PowerSI-RLGC','Original S-parameters'},'Location','best','NumColumns',1)
    legend('boxoff')
end
