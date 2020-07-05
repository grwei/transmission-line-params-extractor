function [] = check_consistence(resistance, inductance, conductance, capacitance, linelength, freq, z0)
%CHECK_CONSISTENCE 检查RLGC参数提取算法s2rlgc_t.m是否引入了非模型误差
%   首先调用rlgc2s_t.m，由形参RLGC计算S参数。然后调用s2rlgc_t.m，由S参数提取RLGC。
% 若两个RLGC参数完全一致，说明s2rlgc_t.m未引入非模型误差。

%%% 由原始RLGC参数求解S参数
[s_params_rebuilt,rlgc_original] = rlgc2s_t(resistance, inductance, conductance, capacitance, linelength, freq, z0);

%%% 由S参数反求RLGC参数
rlgc_re_extracted = s2rlgc_t(s_params_rebuilt,linelength,freq,z0,[],false);

%%% 评估反求的RLGC模型对原始RLGC模型的差距
% under development...

%% 一致性检查：用RLGC参数重建S，再提取新RLGC
% 目的：若用于RLGC参数提取的原始S参数不满足无源、因果、互易等条件，就可能产生RLGC模型误差。
% 本函数从接受一个RLGC模型，求取其S参数；然后反求RLGC模型，与原RLGC模型比较
% 若新RLGC参数与原RLGC参数完全一致，说明s2rlgc_m.t的参数提取过程未引入非模型误差。
% 这样，就可以考虑通过修正原始S参数，以得到更优的RLGC模型。

numOfLines = size(resistance,1);     % Number of transmission lines
% R
figure('Name','R matrix of rebuild S')
sgtitle({'Comparison Between original RLGC and';'re-extracted RLGC using rebuilt S: R Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_original.R(1,idx,:)),'k-')
    hold on
    plot(freq/1e9,squeeze(rlgc_re_extracted.R(1,idx,:)),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('R(1,%u)(Ohms/m)',idx));
    title(sprintf('R(1,%u)',idx));
    legend({'original','re-extracted'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% L
figure('Name','L matrix of rebuild S')
sgtitle({'Comparison Between original RLGC and';'re-extracted RLGC using rebuilt S: L Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_original.L(1,idx,:)),'k-')
    hold on
    plot(freq/1e9,squeeze(rlgc_re_extracted.L(1,idx,:)),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('L(1,%u)(H/m)',idx));
    title(sprintf('L(1,%u)',idx));
    legend({'original','re-extracted'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% G
figure('Name','G matrix of rebuild S')
sgtitle({'Comparison Between original RLGC and';'re-extracted RLGC using rebuilt S: G Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_original.G(1,idx,:)),'k-')
    hold on
    plot(freq/1e9,squeeze(rlgc_re_extracted.G(1,idx,:)),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('G(1,%u)(S/m)',idx));
    title(sprintf('G(1,%u)',idx));
    legend({'original','re-extracted'},'Location','best','NumColumns',1)
    legend('boxoff')
end

% C
figure('Name','C matrix of rebuild S')
sgtitle({'Comparison Between original RLGC and';'re-extracted RLGC using rebuilt S: C Matrix'})
total = ceil(numOfLines/2);
for idx = 1:numOfLines
    subplot(2,total,idx)
    plot(freq/1e9,squeeze(rlgc_original.C(1,idx,:)),'k-')
    hold on
    plot(freq/1e9,squeeze(rlgc_re_extracted.C(1,idx,:)),'g--')
    hold off
    grid on
    xlabel('Freq(GHz)');
    ylabel(sprintf('C(1,%u)(F/m)',idx));
    title(sprintf('C(1,%u)',idx));
    legend({'original','re-extracted'},'Location','best','NumColumns',1)
    legend('boxoff')
end

end
