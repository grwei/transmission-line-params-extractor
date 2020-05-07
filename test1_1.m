%% Basic Information
%%% Overview
% Transmission-line parameters extractor
% MATLAB implementation of Patent US8892414B1
% Author Name: Guorui Wei
% Created in: 2020-05-06 10:45
% Example: Edge-Coupled Surface Microstrip 
% Note: Using larger frequency intervals

clc; clear; close all;

%% Import data
lineLength = 0.0254; % Line Length(meters)

% Import simulated data
filename_2line = 'data/CPL_1in_20200202.s4p';
SingleEnded4PortData = read(rfdata.data,filename_2line);
freq = SingleEnded4PortData.Freq;
freqPts = length(freq);
z0 = SingleEnded4PortData.Z0; % Reference Impedance
SingleEnded4PortData.S_Parameters = snp2smp(SingleEnded4PortData.S_Parameters,...
    z0,[1 2 3 4]); % Classic style
numOfLines = size(SingleEnded4PortData.S_Parameters,1)/2;

% Import Cadence-PowerSI-extracted params
filename_PowerSI = 'data/CPL_1in_20200202_PowerSI.csv';
opts = detectImportOptions(filename_PowerSI);
rlgc_PowerSI_mat = readtable(filename_PowerSI);
% Allocate memory
rlgc_PowerSI.R = zeros(numOfLines,numOfLines,freqPts);
rlgc_PowerSI.L = rlgc_PowerSI.R;
rlgc_PowerSI.C = rlgc_PowerSI.R;
rlgc_PowerSI.G = rlgc_PowerSI.R;
% Load data
for freqIdx = 1:freqPts
    rlgc_PowerSI.R(1,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,3}/lineLength;
    rlgc_PowerSI.R(1,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,4}/lineLength;
    rlgc_PowerSI.R(2,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,4}/lineLength;
    rlgc_PowerSI.R(2,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,5}/lineLength;
    rlgc_PowerSI.L(1,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,3}/lineLength;
    rlgc_PowerSI.L(1,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,4}/lineLength;
    rlgc_PowerSI.L(2,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,4}/lineLength;
    rlgc_PowerSI.L(2,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,5}/lineLength;    
    rlgc_PowerSI.G(1,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,3}/lineLength;
    rlgc_PowerSI.G(1,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,4}/lineLength;
    rlgc_PowerSI.G(2,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,4}/lineLength;
    rlgc_PowerSI.G(2,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,5}/lineLength;
    rlgc_PowerSI.C(1,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,3}/lineLength;
    rlgc_PowerSI.C(1,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,4}/lineLength;
    rlgc_PowerSI.C(2,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,4}/lineLength;
    rlgc_PowerSI.C(2,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,5}/lineLength;
end

%% Extract RLGC params using proposed method

idx_selected = ceil(1:5:freqPts);
rlgc_t = s2rlgc_t(SingleEnded4PortData.S_Parameters(:,:,idx_selected),lineLength,freq(idx_selected),z0,[],true);
check_consistence(rlgc_t.R, rlgc_t.L, rlgc_t.G, rlgc_t.C, lineLength, freq(idx_selected), z0);
