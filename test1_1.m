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

idx_selected = ceil(1:5:freqPts);
rlgc_t = s2rlgc_t(SingleEnded4PortData.S_Parameters(:,:,idx_selected),lineLength,freq(idx_selected),z0,[],true);
% check_consistence(rlgc_t.R, rlgc_t.L, rlgc_t.G, rlgc_t.C, lineLength, freq(idx_selected), z0);
