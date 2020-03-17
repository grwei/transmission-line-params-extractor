%% Basic Information
% Overview
% Transmission-line parameters extractor
% MATLAB implementation of Patent US8892414B1
% Author Name: Guorui Wei
% Instructor Name: Bin Xia
% Created in: 2020-03-15 12:45

clc; clear; close all;

%% Import data

linelength = 0.0254; % Line Length(meters)

% Simulate data
filename_2line = 'data/CPL_1in_20200202.s4p';
SingleEnded4PortData = read(rfdata.data,filename_2line);
freq = SingleEnded4PortData.Freq;
freqPts = length(freq);
z0 = SingleEnded4PortData.Z0; % Reference Impedance
SingleEnded4PortData.S_Parameters = snp2smp(SingleEnded4PortData.S_Parameters,...
    z0,[1 2 3 4]); % Classic style
numOfLines = size(SingleEnded4PortData.S_Parameters,1)/2;
S_params_SI = sparameters(SingleEnded4PortData);

% Cadence-PowerSI-extracted params
filename_PowerSI = 'data/CPL_1in_20200202_PowerSI.csv';
opts = detectImportOptions(filename_PowerSI);
rlgc_PowerSI_mat = readtable(filename_PowerSI);
for freqIdx = 1:freqPts
    rlgc_PowerSI.R(1,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,3}/linelength;
    rlgc_PowerSI.R(1,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,4}/linelength;
    rlgc_PowerSI.R(2,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,4}/linelength;
    rlgc_PowerSI.R(2,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,5}/linelength;
    rlgc_PowerSI.L(1,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,3}/linelength;
    rlgc_PowerSI.L(1,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,4}/linelength;
    rlgc_PowerSI.L(2,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,4}/linelength;
    rlgc_PowerSI.L(2,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,5}/linelength;    
    rlgc_PowerSI.G(1,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,3}/linelength;
    rlgc_PowerSI.G(1,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,4}/linelength;
    rlgc_PowerSI.G(2,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,4}/linelength;
    rlgc_PowerSI.G(2,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,5}/linelength;
    rlgc_PowerSI.C(1,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,3}/linelength;
    rlgc_PowerSI.C(1,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,4}/linelength;
    rlgc_PowerSI.C(2,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,4}/linelength;
    rlgc_PowerSI.C(2,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,5}/linelength;
end


%%
rlgc_t = s2rlgc_t(SingleEnded4PortData.S_Parameters,linelength,freq,z0);

figure('Name','R matrix')
sgtitle('p.u.l Resistence')


