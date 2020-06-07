%% Basic Information
%%% Overview
% Transmission-line parameters extractor
% MATLAB implementation of Patent US8892414B1
% Author Name: Guorui Wei
% Created in: 2020-03-15 12:45

clc; clear; close all;

%% initialize
% 此节只需在每次数据改变后执行一次。
% 读取仿真数据，然后以.mat文件存储到工程根目录，以缩短程序多次运行时读数据时间。

% %%% Import simulated data
% lineLength = 0.00508; % Line Length(meters)
% filename_4line = 'data/4line/4lines_ADS/data/4lines_HFSSW_200mil.s8p';
% SingleEnded8PortData = read(rfdata.data,filename_4line);
% freq = SingleEnded8PortData.Freq;
% freqPts = length(freq);
% numOfLines = size(SingleEnded8PortData.S_Parameters,1)/2;
% 
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
% save('test2_3','SingleEnded8PortData','rlgc_PowerSI','rlgc_HFSSW');

%% process data

load('test2_3','SingleEnded8PortData','rlgc_PowerSI','rlgc_HFSSW');
full_wave = load('test2','SingleEnded8PortData');

% process simulated data
lineLength = 0.00508; % Line Length(meters)
freq = SingleEnded8PortData.Freq;
freqPts = length(freq);
z0 = SingleEnded8PortData.Z0; % Reference Impedance
SingleEnded8PortData.S_Parameters = snp2smp(SingleEnded8PortData.S_Parameters,...
    z0,1:1:8); % Classic style
numOfLines = size(SingleEnded8PortData.S_Parameters,1)/2;

%% Extract RLGC params using proposed method

% extract RLGC
rlgc_t = s2rlgc_t(SingleEnded8PortData.S_Parameters,lineLength,freq,z0,[],false);
rlgc_full_wave = s2rlgc_t(full_wave.SingleEnded8PortData.S_Parameters,lineLength,freq,z0,[],false);

%% Calculate S from HFSS_W-element using proposed method
% compare with the result using ADS 2020.2.2

[s_HFSSW_prop,~] = rlgc2s_t(rlgc_HFSSW.R,rlgc_HFSSW.L,rlgc_HFSSW.G,rlgc_HFSSW.C,lineLength,freq,z0);
% external<-external
figure('Name','Rebuilt s-HFSSW (cmp. with ADS): See') 
sgtitle({'Cmp. between s-HFSSW using prop. method and ADS: See'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    % left axis: dB(S(i,j))
    yyaxis left
    plot(freq/1e9,db(squeeze(s_HFSSW_prop(1,idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(SingleEnded8PortData.S_Parameters(1,idx,:)),'voltage'),'g--')
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('dB(S(1,%u))',idx));
    title(sprintf('S(1,%u)',idx));
    
    % right axis: arg(S(i,j))
    yyaxis right
    plot(freq/1e9,angle(squeeze(s_HFSSW_prop(1,idx,:))),'m-')
    hold on
    plot(freq/1e9,angle(squeeze(SingleEnded8PortData.S_Parameters(1,idx,:))),'c--')
    hold off
    ylabel(sprintf('arg(S(1,%u))',idx));
    legend({'prop.(dB)','ADS(dB)','prop.(arg)','ADS(arg)'},'Location','best','NumColumns',2)
    legend('boxoff')
end

% external<-internal
figure('Name','Rebuilt s-HFSSW (cmp. with ADS): Sei')
sgtitle({'Cmp. between s-HFSSW using prop. method and ADS: Sei'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    % left axis: dB(S(i,j))
    yyaxis left
    plot(freq/1e9,db(squeeze(s_HFSSW_prop(1,idx+numOfLines,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(SingleEnded8PortData.S_Parameters(1,idx+numOfLines,:)), 'voltage'),'g--')
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('dB(S(1,%u))',idx+numOfLines));
    title(sprintf('S(1,%u)',idx+numOfLines));
    
    % right axis: arg(S(i,j))
    yyaxis right
    plot(freq/1e9,angle(squeeze(s_HFSSW_prop(1,idx+numOfLines,:))),'m-')
    hold on
    plot(freq/1e9,angle(squeeze(SingleEnded8PortData.S_Parameters(1,idx+numOfLines,:))),'c--')
    hold off
    ylabel(sprintf('arg(S(1,%u))',idx+numOfLines));
    legend({'prop.(dB)','ADS(dB)','prop.(arg)','ADS(arg)'},'Location','best','NumColumns',2)
    legend('boxoff')
end

% internal<-external
figure('Name','Rebuilt s-HFSSW (cmp. with ADS): Sie') 
sgtitle({'Cmp. between s-HFSSW using prop. method and ADS: Sie'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    % left axis: dB(S(i,j))
    yyaxis left
    plot(freq/1e9,db(squeeze(s_HFSSW_prop(numOfLines+1,idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(SingleEnded8PortData.S_Parameters(numOfLines+1,idx,:)),'voltage'),'g--')
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('dB(S(%u,%u))',numOfLines+1,idx));
    title(sprintf('S(%u,%u)',numOfLines+1,idx));
    
    % right axis: arg(S(i,j))
    yyaxis right
    plot(freq/1e9,angle(squeeze(s_HFSSW_prop(numOfLines+1,idx,:))),'m-')
    hold on
    plot(freq/1e9,angle(squeeze(SingleEnded8PortData.S_Parameters(numOfLines+1,idx,:))),'c--')
    hold off
    ylabel(sprintf('arg(S(%u,%u))',numOfLines+1,idx));
    legend({'prop.(dB)','ADS(dB)','prop.(arg)','ADS(arg)'},'Location','best','NumColumns',2)
    legend('boxoff')
end

% internal<-internal
figure('Name','Rebuilt s-HFSSW (cmp. with ADS): Sii') 
sgtitle({'Cmp. between s-HFSSW using prop. method and ADS: Sii'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    % left axis: dB(S(i,j))
    yyaxis left
    plot(freq/1e9,db(squeeze(s_HFSSW_prop(numOfLines+1,numOfLines+idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(SingleEnded8PortData.S_Parameters(numOfLines+1,numOfLines+idx,:)),'voltage'),'g--')
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('dB(S(%u,%u))',numOfLines+1,numOfLines+idx));
    title(sprintf('S(%u,%u)',numOfLines+1,numOfLines+idx));
    
    % right axis: arg(S(i,j))
    yyaxis right
    plot(freq/1e9,angle(squeeze(s_HFSSW_prop(numOfLines+1,numOfLines+idx,:))),'m-')
    hold on
    plot(freq/1e9,angle(squeeze(SingleEnded8PortData.S_Parameters(numOfLines+1,numOfLines+idx,:))),'c--')
    hold off
    ylabel(sprintf('arg(S(%u,%u))',numOfLines+1,numOfLines+idx));
    legend({'prop.(dB)','ADS(dB)','prop.(arg)','ADS(arg)'},'Location','best','NumColumns',2)
    legend('boxoff')
end

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

%% cmp.of s_full-wave and s_HFSSW

% external<-external
figure('Name','cmp.of s_full-wave and s_HFSSW: See') 
% sgtitle({'cmp.of s_full-wave and s_HFSSW: See'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    % left axis: dB(S(i,j))
    yyaxis left
    plot(freq/1e9,db(squeeze(full_wave.SingleEnded8PortData.S_Parameters(1,idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(SingleEnded8PortData.S_Parameters(1,idx,:)),'voltage'),'g--')
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('dB(S(1,%u))',idx));
    title(sprintf('S(1,%u)',idx));
    
    % right axis: arg(S(i,j))
    yyaxis right
    plot(freq/1e9,angle(squeeze(full_wave.SingleEnded8PortData.S_Parameters(1,idx,:))),'m-')
    hold on
    plot(freq/1e9,angle(squeeze(SingleEnded8PortData.S_Parameters(1,idx,:))),'c--')
    hold off
    ylabel(sprintf('arg(S(1,%u))',idx));
    legend({'full-wave(dB)','W-element(dB)','full-wave(arg)','W-element(arg)'},'Location','best','NumColumns',2)
    legend('boxoff')
end

% external<-internal
figure('Name','cmp.of s_full-wave and s_HFSSW: Sei')
% sgtitle({'cmp.of s_full-wave and s_HFSSW: Sei'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    % left axis: dB(S(i,j))
    yyaxis left
    plot(freq/1e9,db(squeeze(full_wave.SingleEnded8PortData.S_Parameters(1,idx+numOfLines,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(SingleEnded8PortData.S_Parameters(1,idx+numOfLines,:)), 'voltage'),'g--')
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('dB(S(1,%u))',idx+numOfLines));
    title(sprintf('S(1,%u)',idx+numOfLines));
    
    % right axis: arg(S(i,j))
    yyaxis right
    plot(freq/1e9,angle(squeeze(full_wave.SingleEnded8PortData.S_Parameters(1,idx+numOfLines,:))),'m-')
    hold on
    plot(freq/1e9,angle(squeeze(SingleEnded8PortData.S_Parameters(1,idx+numOfLines,:))),'c--')
    hold off
    ylabel(sprintf('arg(S(1,%u))',idx+numOfLines));
    legend({'full-wave(dB)','W-element(dB)','full-wave(arg)','W-element(arg)'},'Location','best','NumColumns',2)
    legend('boxoff')
end

% internal<-external
figure('Name','cmp.of s_full-wave and s_HFSSW: Sie') 
% sgtitle({'cmp.of s_full-wave and s_HFSSW: Sie'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    % left axis: dB(S(i,j))
    yyaxis left
    plot(freq/1e9,db(squeeze(full_wave.SingleEnded8PortData.S_Parameters(numOfLines+1,idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(SingleEnded8PortData.S_Parameters(numOfLines+1,idx,:)),'voltage'),'g--')
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('dB(S(%u,%u))',numOfLines+1,idx));
    title(sprintf('S(%u,%u)',numOfLines+1,idx));
    
    % right axis: arg(S(i,j))
    yyaxis right
    plot(freq/1e9,angle(squeeze(full_wave.SingleEnded8PortData.S_Parameters(numOfLines+1,idx,:))),'m-')
    hold on
    plot(freq/1e9,angle(squeeze(SingleEnded8PortData.S_Parameters(numOfLines+1,idx,:))),'c--')
    hold off
    ylabel(sprintf('arg(S(%u,%u))',numOfLines+1,idx));
    legend({'full-wave(dB)','W-element(dB)','full-wave(arg)','W-element(arg)'},'Location','best','NumColumns',2)
    legend('boxoff')
end

% internal<-internal
figure('Name','cmp.of s_full-wave and s_HFSSW: Sii') 
% sgtitle({'cmp.of s_full-wave and s_HFSSW: Sii'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    % left axis: dB(S(i,j))
    yyaxis left
    plot(freq/1e9,db(squeeze(full_wave.SingleEnded8PortData.S_Parameters(numOfLines+1,numOfLines+idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(SingleEnded8PortData.S_Parameters(numOfLines+1,numOfLines+idx,:)),'voltage'),'g--')
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('dB(S(%u,%u))',numOfLines+1,numOfLines+idx));
    title(sprintf('S(%u,%u)',numOfLines+1,numOfLines+idx));
    
    % right axis: arg(S(i,j))
    yyaxis right
    plot(freq/1e9,angle(squeeze(full_wave.SingleEnded8PortData.S_Parameters(numOfLines+1,numOfLines+idx,:))),'m-')
    hold on
    plot(freq/1e9,angle(squeeze(SingleEnded8PortData.S_Parameters(numOfLines+1,numOfLines+idx,:))),'c--')
    hold off
    ylabel(sprintf('arg(S(%u,%u))',numOfLines+1,numOfLines+idx));
    legend({'full-wave(dB)','W-element(dB)','full-wave(arg)','W-element(arg)'},'Location','best','NumColumns',2)
    legend('boxoff')
end

%%  cmp. of ext. RLGC from full-wave S and W-elem S

% R
figure('Name','R matrix (cmp. full-wave-s and W-elem-s)')
% sgtitle({'Comparison Between Proposed Algorithm and';' HFSS: R Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_t.R(1,idx,:)),'k-')
    hold on
    plot(freq/1e9,squeeze(rlgc_full_wave.R(1,idx,:)),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('R(1,%u) (Ohms/m)',idx));
    title(sprintf('R(1,%u)',idx));
    legend({'W-element S','full-wave S'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% L
figure('Name','L matrix (cmp. full-wave-s and W-elem-s)')
% sgtitle({'Comparison Between Proposed Algorithm and';'HFSS: L Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_t.L(1,idx,:)),'k-')
    hold on
    plot(freq/1e9,squeeze(rlgc_full_wave.L(1,idx,:)),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('L(1,%u) (H/m)',idx));
    title(sprintf('L(1,%u)',idx));
    legend({'W-element S','full-wave S'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% G
figure('Name','G matrix (cmp. full-wave-s and W-elem-s)')
% sgtitle({'Comparison Between Proposed Algorithm and';'HFSS: G Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_t.G(1,idx,:)),'k-')
    hold on
    plot(freq/1e9,squeeze(rlgc_full_wave.G(1,idx,:)),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('G(1,%u) (S/m)',idx));
    title(sprintf('G(1,%u)',idx));
    legend({'W-element S','full-wave S'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% C
figure('Name','C matrix (cmp. full-wave-s and W-elem-s)')
% sgtitle({'Comparison Between Proposed Algorithm and';'HFSS: C Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_t.C(1,idx,:)),'k-')
    hold on
    plot(freq/1e9,squeeze(rlgc_full_wave.C(1,idx,:)),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('C(1,%u) (F/m)',idx));
    title(sprintf('C(1,%u)',idx));
    legend({'W-element S','full-wave S'},'Location','best','NumColumns',1)
    legend('boxoff')
end

%% cmp. rebuild s_rlgc_t and s_rlgc_full_wave

% rebuild S
[s_rlgc_t,~] = rlgc2s_t(rlgc_t.R,rlgc_t.L,rlgc_t.G,rlgc_t.C,lineLength,freq,z0);
[s_rlgc_full_wave,~] = rlgc2s_t(rlgc_full_wave.R,rlgc_full_wave.L,rlgc_full_wave.G,rlgc_full_wave.C,lineLength,freq,z0);

% plot
% external<-external
figure('Name','cmp. rebuild s_rlgc_t and s_rlgc_full_wave See') 
% sgtitle({'Cmp. between s-HFSSW using prop. method and ADS: See'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    % left axis: dB(S(i,j))
    yyaxis left
    plot(freq/1e9,db(squeeze(s_rlgc_t(1,idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(s_rlgc_full_wave(1,idx,:)),'voltage'),'g--')
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('dB(S(1,%u))',idx));
    title(sprintf('S(1,%u)',idx));
    
    % right axis: arg(S(i,j))
    yyaxis right
    plot(freq/1e9,angle(squeeze(s_rlgc_t(1,idx,:))),'m-')
    hold on
    plot(freq/1e9,angle(squeeze(s_rlgc_full_wave(1,idx,:))),'c--')
    hold off
    ylabel(sprintf('arg(S(1,%u))',idx));
    legend({'W-element(dB)','full-wave(dB)','W-element(arg)','full-wave(arg)'},'Location','best','NumColumns',2)
    legend('boxoff')
end

% external<-internal
figure('Name','cmp. rebuild s_rlgc_t and s_rlgc_full_wave Sei')
% sgtitle({'Cmp. between s-HFSSW using prop. method and ADS: Sei'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    % left axis: dB(S(i,j))
    yyaxis left
    plot(freq/1e9,db(squeeze(s_rlgc_t(1,idx+numOfLines,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(s_rlgc_full_wave(1,idx+numOfLines,:)), 'voltage'),'g--')
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('dB(S(1,%u))',idx+numOfLines));
    title(sprintf('S(1,%u)',idx+numOfLines));
    
    % right axis: arg(S(i,j))
    yyaxis right
    plot(freq/1e9,angle(squeeze(s_rlgc_t(1,idx+numOfLines,:))),'m-')
    hold on
    plot(freq/1e9,angle(squeeze(s_rlgc_full_wave(1,idx+numOfLines,:))),'c--')
    hold off
    ylabel(sprintf('arg(S(1,%u))',idx+numOfLines));
    legend({'W-element(dB)','full-wave(dB)','W-element(arg)','full-wave(arg)'},'Location','best','NumColumns',2)
    legend('boxoff')
end

% internal<-external
figure('Name','cmp. rebuild s_rlgc_t and s_rlgc_full_wave Sie') 
% sgtitle({'Cmp. between s-HFSSW using prop. method and ADS: Sie'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    % left axis: dB(S(i,j))
    yyaxis left
    plot(freq/1e9,db(squeeze(s_rlgc_t(numOfLines+1,idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(s_rlgc_full_wave(numOfLines+1,idx,:)),'voltage'),'g--')
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('dB(S(%u,%u))',numOfLines+1,idx));
    title(sprintf('S(%u,%u)',numOfLines+1,idx));
    
    % right axis: arg(S(i,j))
    yyaxis right
    plot(freq/1e9,angle(squeeze(s_rlgc_t(numOfLines+1,idx,:))),'m-')
    hold on
    plot(freq/1e9,angle(squeeze(s_rlgc_full_wave(numOfLines+1,idx,:))),'c--')
    hold off
    ylabel(sprintf('arg(S(%u,%u))',numOfLines+1,idx));
    legend({'W-element(dB)','full-wave(dB)','W-element(arg)','full-wave(arg)'},'Location','best','NumColumns',2)
    legend('boxoff')
end

% internal<-internal
figure('Name','cmp. rebuild s_rlgc_t and s_rlgc_full_wave Sii') 
% sgtitle({'Cmp. between s-HFSSW using prop. method and ADS: Sii'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    % left axis: dB(S(i,j))
    yyaxis left
    plot(freq/1e9,db(squeeze(s_rlgc_t(numOfLines+1,numOfLines+idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(s_rlgc_full_wave(numOfLines+1,numOfLines+idx,:)),'voltage'),'g--')
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('dB(S(%u,%u))',numOfLines+1,numOfLines+idx));
    title(sprintf('S(%u,%u)',numOfLines+1,numOfLines+idx));
    
    % right axis: arg(S(i,j))
    yyaxis right
    plot(freq/1e9,angle(squeeze(s_rlgc_t(numOfLines+1,numOfLines+idx,:))),'m-')
    hold on
    plot(freq/1e9,angle(squeeze(s_rlgc_full_wave(numOfLines+1,numOfLines+idx,:))),'c--')
    hold off
    ylabel(sprintf('arg(S(%u,%u))',numOfLines+1,numOfLines+idx));
    legend({'W-element(dB)','full-wave(dB)','W-element(arg)','full-wave(arg)'},'Location','best','NumColumns',2)
    legend('boxoff')
end
