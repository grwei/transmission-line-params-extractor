%% Basic Information
%%% Overview
% Transmission-line parameters extractor
% MATLAB implementation of Patent US8892414B1
% Author Name: Guorui Wei
% Created in: 2020-03-15 12:45
% Example: Edge-Coupled Surface Microstrip 

clc; clear; close all;

%% Initialize
% 此节只需在每次数据改变后执行一次。
% 读取仿真数据，然后以.mat文件存储到工程根目录，以缩短程序多次运行时读数据时间。

% %%% Import simulated data
% lineLength = 0.0127; % Line Length(meters)
% filename_2line = 'data/2line/2lines_Polar_500mil.s4p';
% SingleEnded4PortData = read(rfdata.data,filename_2line);
% freq = SingleEnded4PortData.Freq;
% freqPts = length(freq);
% z0 = SingleEnded4PortData.Z0; % Reference Impedance
% numOfLines = size(SingleEnded4PortData.S_Parameters,1)/2;
% SingleEnded4PortData.S_Parameters = snp2smp(SingleEnded4PortData.S_Parameters,...
%     z0,1:1:2*numOfLines); % Classic style
% 
% %%% Import Cadence-PowerSI-extracted params
% % Allocate memory
% rlgc_PowerSI.R = zeros(numOfLines,numOfLines,freqPts);
% rlgc_PowerSI.L = rlgc_PowerSI.R;
% rlgc_PowerSI.C = rlgc_PowerSI.R;
% rlgc_PowerSI.G = rlgc_PowerSI.R;
% % Load data
% filename_PowerSI = 'data/2line/Transmission_RLGC_res.csv';
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
% %%% import Polar-RLGC
% filename_Polar = 'data/2line/PolarSi.xlsb';
% [status,sheets_Polar,xlFormat] = xlsfinfo(filename_Polar);
% odd_Polar = readtable(filename_Polar,'sheet',string(sheets_Polar(1))); % odd mode
% even_Polar = readtable(filename_Polar,'sheet',string(sheets_Polar(2))); % even mode
% % 
% Zc_odd_Polar = complex(odd_Polar{:,2},odd_Polar{:,3});
% gamma_odd_Polar = complex(odd_Polar{:,16},odd_Polar{:,17});
% Zc_even_Polar = complex(even_Polar{:,2},even_Polar{:,3});
% gamma_even_Polar = complex(even_Polar{:,16},even_Polar{:,17});
% 
% save('test1')

%% load data
load('test1')

%% Extract RLGC params using proposed method

rlgc_t = s2rlgc_t(SingleEnded4PortData.S_Parameters,lineLength,freq,z0,[],false);
% check_consistence(rlgc_t.R, rlgc_t.L, rlgc_t.G, rlgc_t.C, lineLength, freq, z0);
[s_rlgc_PowerSI,~] = rlgc2s_t(rlgc_PowerSI.R,rlgc_PowerSI.L,rlgc_PowerSI.G,rlgc_PowerSI.C,lineLength,freq,z0);

%% Calculate even-/odd- mode Zc and gamma using rlgc_t (method 1, by RLGC definition)
% pre-allocate memory
z_s_prop = zeros(freqPts,1); % self impedance p.u.l., i.e. z11 or z12, p.u.l. 
z_m_prop = z_s_prop;              % mutual impedance p.u.l., i.e. z12 or z21, p.u.l. 
y_s_prop = z_s_prop;              % self admittance p.u.l., i.e. y11 or y22, p.u.l. 
y_m_prop = z_s_prop;              % mutual admittance p.u.l., i.e. y12 or y21, p.u.l.

% build self-/mutual- impedance/admittance
for freqIdx = 1:freqPts
    z_s_prop(freqIdx) = complex(rlgc_t.R(1,1,freqIdx),2*pi*freq(freqIdx)*rlgc_t.L(1,1,freqIdx));
    z_m_prop(freqIdx) = complex(rlgc_t.R(1,2,freqIdx),2*pi*freq(freqIdx)*rlgc_t.L(1,2,freqIdx));
    y_s_prop(freqIdx) = complex(rlgc_t.G(1,1,freqIdx),2*pi*freq(freqIdx)*rlgc_t.C(1,1,freqIdx));
    y_m_prop(freqIdx) = complex(rlgc_t.G(1,2,freqIdx),2*pi*freq(freqIdx)*rlgc_t.C(1,2,freqIdx));
end

% calculate odd/even mode Zc and gamma
Zc_odd_prop = sqrt((z_s_prop - z_m_prop)./(y_s_prop - y_m_prop));
Zc_even_prop = sqrt((z_s_prop + z_m_prop)./(y_s_prop + y_m_prop));
gamma_odd_prop = sqrt((z_s_prop - z_m_prop).*(y_s_prop - y_m_prop));
gamma_even_prop = sqrt((z_s_prop + z_m_prop).*(y_s_prop + y_m_prop));

%% Calculate even-/odd- mode Zc and gamma using rlgc_t (method 2, by Zc/gamma property)
Zc_odd_prop2 = squeeze(rlgc_t.Zc(1,1,:) - rlgc_t.Zc(1,2,:));
gamma_odd_prop2 = squeeze(complex(rlgc_t.alpha(1,1,:) - rlgc_t.alpha(1,2,:),rlgc_t.beta(1,1,:) - rlgc_t.beta(1,2,:)));
Zc_even_prop2 = squeeze(rlgc_t.Zc(1,1,:) + rlgc_t.Zc(1,2,:));
gamma_even_prop2 = squeeze(complex(rlgc_t.alpha(1,1,:) + rlgc_t.alpha(1,2,:),rlgc_t.beta(1,1,:) + rlgc_t.beta(1,2,:)));

%% Calculate even-/odd- mode Zc and gamma using rlgc_PowerSI (only method 1 available)
% because \beta is not available in PowerSI-extracted params

% pre-allocate memory
z_s_PowerSI = zeros(freqPts,1); % self impedance p.u.l., i.e. z11 or z12, p.u.l. 
z_m_PowerSI = z_s_PowerSI;              % mutual impedance p.u.l., i.e. z12 or z21, p.u.l. 
y_s_PowerSI = z_s_PowerSI;              % self admittance p.u.l., i.e. y11 or y22, p.u.l. 
y_m_PowerSI = z_s_PowerSI;              % mutual admittance p.u.l., i.e. y12 or y21, p.u.l.

% build self-/mutual- impedance/admittance
for freqIdx = 1:freqPts
    z_s_PowerSI(freqIdx) = complex(rlgc_t.R(1,1,freqIdx),2*pi*freq(freqIdx)*rlgc_t.L(1,1,freqIdx));
    z_m_PowerSI(freqIdx) = complex(rlgc_t.R(1,2,freqIdx),2*pi*freq(freqIdx)*rlgc_t.L(1,2,freqIdx));
    y_s_PowerSI(freqIdx) = complex(rlgc_t.G(1,1,freqIdx),2*pi*freq(freqIdx)*rlgc_t.C(1,1,freqIdx));
    y_m_PowerSI(freqIdx) = complex(rlgc_t.G(1,2,freqIdx),2*pi*freq(freqIdx)*rlgc_t.C(1,2,freqIdx));
end

% calculate odd/even mode Zc and gamma
Zc_odd_PowerSI = sqrt((z_s_PowerSI - z_m_PowerSI)./(y_s_PowerSI - y_m_PowerSI));
Zc_even_PowerSI = sqrt((z_s_PowerSI + z_m_PowerSI)./(y_s_PowerSI + y_m_PowerSI));
gamma_odd_PowerSI = sqrt((z_s_PowerSI - z_m_PowerSI).*(y_s_PowerSI - y_m_PowerSI));
gamma_even_PowerSI = sqrt((z_s_PowerSI + z_m_PowerSI).*(y_s_PowerSI + y_m_PowerSI));

%% Cmp. between Zc/gamma prop.extraction method 1 and 2

figure('Name','cmp. of Zc/gamma prop.ext. method 1 and 2')
%%% Zc_odd
subplot(2,2,1)
% left axis: real(Zc_odd_prop)
yyaxis left
plot(freq/1e9,real(Zc_odd_prop),'k-')
hold on
plot(freq/1e9,real(Zc_odd_prop2),'g--')
grid on
xlabel('Freq(GHz)');
ylabel('\Re(Zc,o) (Ohms)');
title('Zc,o');
% right axis: imag(Zc_odd_prop)
yyaxis right
plot(freq/1e9,imag(Zc_odd_prop),'m-')
plot(freq/1e9,imag(Zc_odd_prop2),'c--');
hold off
ylabel('\Im(Zc,o) (Ohms)');
legend({'\Re(Zc,o)-prop.1','\Re(Zc,o)-prop.2','\Im(Zc,o)-prop.1','\Im(Zc,o)-prop.2'},'Location','best','NumColumns',2)
legend('boxoff')

%%% Zc_even
subplot(2,2,3)
% left axis: real(Zc_even_prop)
yyaxis left
plot(freq/1e9,real(Zc_even_prop),'k-')
hold on
plot(freq/1e9,real(Zc_even_prop2),'g--')
grid on
xlabel('Freq(GHz)');
ylabel('\Re(Zc,e) (Ohms)');
title('Zc,e');
% right axis: imag(Zc_even_prop)
yyaxis right
plot(freq/1e9,imag(Zc_even_prop),'m-')
plot(freq/1e9,imag(Zc_even_prop2),'c--');
hold off
ylabel('\Im(Zc,e) (Ohms)');
legend({'\Re(Zc,e)-prop.1','\Re(Zc,e)-prop.2','\Im(Zc,e)-prop.1','\Im(Zc,e)-prop.2'},'Location','best','NumColumns',2)
legend('boxoff')

%%% gamma_odd
subplot(2,2,2)
% left axis: alpha_o
yyaxis left
plot(freq/1e9,real(gamma_odd_prop),'k-')
hold on
plot(freq/1e9,real(gamma_odd_prop2),'g--')
grid on
xlabel('Freq(GHz)');
ylabel('\alpha,o (Np/m)');
title('\gamma,o');
% right axis: beta_o
yyaxis right
plot(freq/1e9,imag(gamma_odd_prop),'m-')
plot(freq/1e9,imag(gamma_odd_prop2),'c--');
hold off
ylabel('\beta,o (rad/m)');
legend({'\alpha,o-prop.1','\alpha,o-prop.2','\beta,o-prop.1','\beta,o-prop.2'},'Location','best','NumColumns',2)
legend('boxoff')

%%% gamma_even
subplot(2,2,4)
% left axis: alpha_e
yyaxis left
plot(freq/1e9,real(gamma_even_prop),'k-')
hold on
plot(freq/1e9,real(gamma_even_prop2),'g--')
grid on
xlabel('Freq(GHz)');
ylabel('\alpha,e (Np/m)');
title('\gamma,e');
% right axis: beta_o
yyaxis right
plot(freq/1e9,imag(gamma_even_prop),'m-')
plot(freq/1e9,imag(gamma_even_prop2),'c--');
hold off
ylabel('\beta,e (rad/m)');
legend({'\alpha,e-prop.1','\alpha,e-prop.2','\beta,e-prop.1','\beta,e-prop.2'},'Location','best','NumColumns',2)
legend('boxoff')

%% Cmp. between Zc/gamma prop. and Polar

figure('Name','cmp. of Zc/gamma prop. and Polar')
%%% Zc_odd
subplot(2,2,1)
% left axis: real(Zc_odd_prop)
yyaxis left
plot(freq/1e9,real(Zc_odd_prop),'k-')
hold on
plot(freq/1e9,real(Zc_odd_Polar),'g--')
grid on
xlabel('Freq(GHz)');
ylabel('\Re(Zc,o) (Ohms)');
title('Zc,o');
% right axis: imag(Zc_odd_prop)
yyaxis right
plot(freq/1e9,imag(Zc_odd_prop),'m-')
plot(freq/1e9,imag(Zc_odd_Polar),'c--');
hold off
ylabel('\Im(Zc,o) (Ohms)');
legend({'\Re(Zc,o)-prop.','\Re(Zc,o)-Polar','\Im(Zc,o)-prop.','\Im(Zc,o)-Polar'},'Location','best','NumColumns',2)
legend('boxoff')

%%% Zc_even
subplot(2,2,3)
% left axis: real(Zc_even_prop)
yyaxis left
plot(freq/1e9,real(Zc_even_prop),'k-')
hold on
plot(freq/1e9,real(Zc_even_Polar),'g--')
grid on
xlabel('Freq(GHz)');
ylabel('\Re(Zc,e) (Ohms)');
title('Zc,e');
% right axis: imag(Zc_even_prop)
yyaxis right
plot(freq/1e9,imag(Zc_even_prop),'m-')
plot(freq/1e9,imag(Zc_even_Polar),'c--');
hold off
ylabel('\Im(Zc,e) (Ohms)');
legend({'\Re(Zc,e)-prop.','\Re(Zc,e)-Polar','\Im(Zc,e)-prop.','\Im(Zc,e)-Polar'},'Location','best','NumColumns',2)
legend('boxoff')

%%% gamma_odd
subplot(2,2,2)
% left axis: alpha_o
yyaxis left
plot(freq/1e9,real(gamma_odd_prop),'k-')
hold on
plot(freq/1e9,real(gamma_odd_Polar),'g--')
grid on
xlabel('Freq(GHz)');
ylabel('\alpha,o (Np/m)');
title('\gamma,o');
% right axis: beta_o
yyaxis right
plot(freq/1e9,imag(gamma_odd_prop),'m-')
plot(freq/1e9,imag(gamma_odd_Polar),'c--');
hold off
ylabel('\beta,o (rad/m)');
legend({'\alpha,o-prop.','\alpha,o-Polar','\beta,o-prop.','\beta,o-Polar'},'Location','best','NumColumns',2)
legend('boxoff')

%%% gamma_even
subplot(2,2,4)
% left axis: alpha_e
yyaxis left
plot(freq/1e9,real(gamma_even_prop),'k-')
hold on
plot(freq/1e9,real(gamma_even_Polar),'g--')
grid on
xlabel('Freq(GHz)');
ylabel('\alpha,e (Np/m)');
title('\gamma,e');
% right axis: beta_o
yyaxis right
plot(freq/1e9,imag(gamma_even_prop),'m-')
plot(freq/1e9,imag(gamma_even_Polar),'c--');
hold off
ylabel('\beta,e (rad/m)');
legend({'\alpha,e-prop.','\alpha,e-Polar','\beta,e-prop.','\beta,e-Polar'},'Location','best','NumColumns',2)
legend('boxoff')

%% Cmp. between Zc/gamma prop. and PowerSI

figure('Name','cmp. of Zc/gamma prop. and PowerSI')
%%% Zc_odd
subplot(2,2,1)
% left axis: real(Zc_odd_prop)
yyaxis left
plot(freq/1e9,real(Zc_odd_prop),'k-')
hold on
plot(freq/1e9,real(Zc_odd_PowerSI),'g--')
grid on
xlabel('Freq(GHz)');
ylabel('\Re(Zc,o) (Ohms)');
title('Zc,o');
% right axis: imag(Zc_odd_prop)
yyaxis right
plot(freq/1e9,imag(Zc_odd_prop),'m-')
plot(freq/1e9,imag(Zc_odd_PowerSI),'c--');
hold off
ylabel('\Im(Zc,o) (Ohms)');
legend({'\Re(Zc,o)-prop.1','\Re(Zc,o)-PowerSI','\Im(Zc,o)-prop.1','\Im(Zc,o)-PowerSI'},'Location','best','NumColumns',2)
legend('boxoff')

%%% Zc_even
subplot(2,2,3)
% left axis: real(Zc_even_prop)
yyaxis left
plot(freq/1e9,real(Zc_even_prop),'k-')
hold on
plot(freq/1e9,real(Zc_even_PowerSI),'g--')
grid on
xlabel('Freq(GHz)');
ylabel('\Re(Zc,e) (Ohms)');
title('Zc,e');
% right axis: imag(Zc_even_prop)
yyaxis right
plot(freq/1e9,imag(Zc_even_prop),'m-')
plot(freq/1e9,imag(Zc_even_PowerSI),'c--');
hold off
ylabel('\Im(Zc,e) (Ohms)');
legend({'\Re(Zc,e)-prop.1','\Re(Zc,e)-PowerSI','\Im(Zc,e)-prop.1','\Im(Zc,e)-PowerSI'},'Location','best','NumColumns',2)
legend('boxoff')

%%% gamma_odd
subplot(2,2,2)
% left axis: alpha_o
yyaxis left
plot(freq/1e9,real(gamma_odd_prop),'k-')
hold on
plot(freq/1e9,real(gamma_odd_PowerSI),'g--')
grid on
xlabel('Freq(GHz)');
ylabel('\alpha,o (Np/m)');
title('\gamma,o');
% right axis: beta_o
yyaxis right
plot(freq/1e9,imag(gamma_odd_prop),'m-')
plot(freq/1e9,imag(gamma_odd_PowerSI),'c--');
hold off
ylabel('\beta,o (rad/m)');
legend({'\alpha,o-prop.1','\alpha,o-PowerSI','\beta,o-prop.1','\beta,o-PowerSI'},'Location','best','NumColumns',2)
legend('boxoff')

%%% gamma_even
subplot(2,2,4)
% left axis: alpha_e
yyaxis left
plot(freq/1e9,real(gamma_even_prop),'k-')
hold on
plot(freq/1e9,real(gamma_even_PowerSI),'g--')
grid on
xlabel('Freq(GHz)');
ylabel('\alpha,e (Np/m)');
title('\gamma,e');
% right axis: beta_o
yyaxis right
plot(freq/1e9,imag(gamma_even_prop),'m-')
plot(freq/1e9,imag(gamma_even_PowerSI),'c--');
hold off
ylabel('\beta,e (rad/m)');
legend({'\alpha,e-prop.1','\alpha,e-PowerSI','\beta,e-prop.1','\beta,e-PowerSI'},'Location','best','NumColumns',2)
legend('boxoff')

%% Extracted RLGC compared with PowerSI

% R, L
figure('Name','R, L matrix')
sgtitle({'Comparison Between prop. and';'PowerSI: R, L Matrix'})
subplot(221)
plot(freq/1e9,squeeze(rlgc_t.R(1,1,:)),'k-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.R(1,1,:)),'g--')
hold off
grid on
xlabel('Freq(GHz)');
ylabel('R11(Ohms/m)');
title('R11 Comparison');
legend({'prop.','PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

subplot(222)
plot(freq/1e9,squeeze(rlgc_t.R(1,2,:)),'k-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.R(1,2,:)),'g--')
hold off
grid on
xlabel('Freq(GHz)');
ylabel('R12(Ohms/m)');
title('R12 Comparison');
legend({'prop.','PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

subplot(223)
plot(freq/1e9,squeeze(rlgc_t.L(1,1,:)),'k-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.L(1,1,:)),'g--')
hold off
grid on
xlabel('Freq(GHz)');
ylabel('L11(H/m)');
title('L11 Comparison');
legend({'prop.','PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

subplot(224)
plot(freq/1e9,squeeze(rlgc_t.L(1,2,:)),'k-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.L(1,2,:)),'g--')
hold off
grid on
xlabel('Freq(GHz)');
ylabel('L12(H/m)');
title('L12 Comparison');
legend({'prop.','PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

% C, G
figure('Name','C, G matrix')
sgtitle({'Comparison Between prop. and';' PowerSI: C, G Matrix'})
subplot(221)
plot(freq/1e9,squeeze(rlgc_t.C(1,1,:)),'m-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.C(1,1,:)),'k--')
hold off
grid on
% ave_y = mean(squeeze(rlgc_t.C(1,1,:)));
% del_r = 4e-5;
% ylim([ave_y*(1-del_r), ave_y*(1+del_r)])
xlabel('Freq(GHz)');
ylabel('C11(F/m)');
title('C11 Comparison');
legend({'prop.','PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

subplot(222)
plot(freq/1e9,squeeze(rlgc_t.C(1,2,:)),'m-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.C(1,2,:)),'k--')
hold off
grid on
% ave_y = mean(squeeze(rlgc_t.C(1,2,:)));
% del_r = 5e-4;
% ylim([ave_y*(1+del_r), ave_y*(1-del_r)])
xlabel('Freq(GHz)');
ylabel('C12(F/m)');
title('C12 Comparison');
legend({'prop.','PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

subplot(223)
plot(freq/1e9,squeeze(rlgc_t.G(1,1,:)),'m-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.G(1,1,:)),'k--')
hold off
grid on
xlabel('Freq(GHz)');
ylabel('G11(S/m)');
title('G11 Comparison');
legend({'prop.','PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

subplot(224)
plot(freq/1e9,squeeze(rlgc_t.G(1,2,:)),'m-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.G(1,2,:)),'k--')
hold off
grid on
xlabel('Freq(GHz)');
ylabel('G12(S/m)');
title('G12 Comparison');
legend({'prop.','PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

%% Calculated S-parameters using powerSI-extracted RLGC
% external<-external
figure('Name','s_rlgc_PowerSI: See') 
sgtitle({'Comparison Between Calculated S-parameters using';'PowerSI-RLGC and Original S-parameters: See'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    % left axis: dB(S(i,j))
    yyaxis left
    plot(freq/1e9,db(squeeze(s_rlgc_PowerSI(1,idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(SingleEnded4PortData.S_Parameters(1,idx,:)),'voltage'),'g--')
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('dB(S(1,%u))',idx));
    title(sprintf('S(1,%u)',idx));
    
    % right axis: arg(S(i,j))
    yyaxis right
    plot(freq/1e9,angle(squeeze(s_rlgc_PowerSI(1,idx,:))),'m-')
    hold on
    plot(freq/1e9,angle(squeeze(SingleEnded4PortData.S_Parameters(1,idx,:))),'c--')
    hold off
    ylabel(sprintf('arg(S(1,%u))',idx));
    legend({'PowerSI(dB)','original(dB)','PowerSI(arg)','original(arg)'},'Location','best','NumColumns',2)
    legend('boxoff')
end

% external<-internal
figure('Name','s_rlgc_PowerSI: Sei')
sgtitle({'Comparison Between Calculated S-parameters using';'PowerSI-RLGC and Original S-parameters: Sei'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    % left axis: dB(S(i,j))
    yyaxis left
    plot(freq/1e9,db(squeeze(s_rlgc_PowerSI(1,idx+numOfLines,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(SingleEnded4PortData.S_Parameters(1,idx+numOfLines,:)), 'voltage'),'g--')
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('dB(S(1,%u))',idx+numOfLines));
    title(sprintf('S(1,%u)',idx+numOfLines));
    
    % right axis: arg(S(i,j))
    yyaxis right
    plot(freq/1e9,angle(squeeze(s_rlgc_PowerSI(1,idx+numOfLines,:))),'m-')
    hold on
    plot(freq/1e9,angle(squeeze(SingleEnded4PortData.S_Parameters(1,idx+numOfLines,:))),'c--')
    hold off
    ylabel(sprintf('arg(S(1,%u))',idx+numOfLines));
    legend({'PowerSI(dB)','original(dB)','PowerSI(arg)','original(arg)'},'Location','best','NumColumns',2)
    legend('boxoff')
end

% internal<-external
figure('Name','s_rlgc_PowerSI: Sie') 
sgtitle({'Comparison Between Calculated S-parameters using';'PowerSI-RLGC and Original S-parameters: Sie'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    % left axis: dB(S(i,j))
    yyaxis left
    plot(freq/1e9,db(squeeze(s_rlgc_PowerSI(numOfLines+1,idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(SingleEnded4PortData.S_Parameters(numOfLines+1,idx,:)),'voltage'),'g--')
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('dB(S(%u,%u))',numOfLines+1,idx));
    title(sprintf('S(%u,%u)',numOfLines+1,idx));
    
    % right axis: arg(S(i,j))
    yyaxis right
    plot(freq/1e9,angle(squeeze(s_rlgc_PowerSI(numOfLines+1,idx,:))),'m-')
    hold on
    plot(freq/1e9,angle(squeeze(SingleEnded4PortData.S_Parameters(numOfLines+1,idx,:))),'c--')
    hold off
    ylabel(sprintf('arg(S(%u,%u))',numOfLines+1,idx));
    legend({'PowerSI(dB)','original(dB)','PowerSI(arg)','original(arg)'},'Location','best','NumColumns',2)
    legend('boxoff')
end

% internal<-internal
figure('Name','s_rlgc_PowerSI: Sii') 
sgtitle({'Comparison Between Calculated S-parameters using';'PowerSI-RLGC and Original S-parameters: Sii'})
num_of_columes = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,num_of_columes,idx)
    % left axis: dB(S(i,j))
    yyaxis left
    plot(freq/1e9,db(squeeze(s_rlgc_PowerSI(numOfLines+1,numOfLines+idx,:)),'voltage'),'k-')
    hold on
    plot(freq/1e9,db(squeeze(SingleEnded4PortData.S_Parameters(numOfLines+1,numOfLines+idx,:)),'voltage'),'g--')
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('dB(S(%u,%u))',numOfLines+1,numOfLines+idx));
    title(sprintf('S(%u,%u)',numOfLines+1,numOfLines+idx));
    
    % right axis: arg(S(i,j))
    yyaxis right
    plot(freq/1e9,angle(squeeze(s_rlgc_PowerSI(numOfLines+1,numOfLines+idx,:))),'m-')
    hold on
    plot(freq/1e9,angle(squeeze(SingleEnded4PortData.S_Parameters(numOfLines+1,numOfLines+idx,:))),'c--')
    hold off
    ylabel(sprintf('arg(S(%u,%u))',numOfLines+1,numOfLines+idx));
    legend({'PowerSI(dB)','original(dB)','PowerSI(arg)','original(arg)'},'Location','best','NumColumns',2)
    legend('boxoff')
end
