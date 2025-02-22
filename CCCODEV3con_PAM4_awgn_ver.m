clear all;
close all;
clc;

% 定义常量
modu_order = 4;      % 使用 PAM4 调制
EbNo_all = 0:12; % SNR 范围
numSymPerFrame = 50000; % 每帧符号数
tbdepth = 16; 
bits_per_symbol = log2(modu_order); % 每个符号的比特数

% 定义五种卷积码结构以及对应的 K 和 N 值
trellis_structs = {
    struct('trellis', poly2trellis(3, [5 7]), 'K', 1, 'N', 2),           % (2,1,3)
    struct('trellis', poly2trellis(7, [171 133]), 'K', 1, 'N', 2),       % (2,1,7)
    struct('trellis', poly2trellis([5 ,4], [23,35,0;0,05,13]), 'K', 2, 'N', 3), % (3,2,4)
    struct('trellis', poly2trellis([4, 4 ,4], [11,13,0,0;0,05,11,0;0,0,07,13]), 'K', 3, 'N', 4), % (4,3,4)
    struct('trellis', poly2trellis([5, 4, 4], [23 35 0 0; 0 05 13 0; 0 0 05 13]), 'K', 3, 'N', 4) % (4,3,5)
};

% 生成固定的随机输入比特序列
%dataIn = randi([0 1], 3*numSymPerFrame * bits_per_symbol, 1);
%讀取固定檔案
dataIn = load('dataIn.asv', '-ascii');
% dataIn = randi([0 1], 1, 8191);
for t = 1:length(trellis_structs)
    trellis = trellis_structs{t}.trellis;
    K = trellis_structs{t}.K;
    N = trellis_structs{t}.N;
    for i = 1:length(EbNo_all)
            EbNo = EbNo_all(i);
            snr_all(t, i) = EbNo + 10*log10(bits_per_symbol * (K/N));
    end
end

% 初始化 BER 矩阵
berCoded = zeros(length(trellis_structs), length(EbNo_all));
berUncoded = zeros(length(trellis_structs), length(EbNo_all));
execution_times = zeros(1, length(trellis_structs)); % 初始化执行时间向量

for t = 1:length(trellis_structs)
    trellis = trellis_structs{t}.trellis;
    K = trellis_structs{t}.K;
    N = trellis_structs{t}.N;
    
    % 开始计时
    tic;
    
    for i = 1:length(EbNo_all)
        %EbNo = EbNo_all(i);
        %snr = EbNo + 10*log10(bits_per_symbol * (K/N));
        
        % 对固定的 dataIn 进行卷积编码
        coded_data = convenc(dataIn', trellis); % 卷积编码
        
        symbol_indices = bi2de(reshape(coded_data, bits_per_symbol, []).', 'left-msb');

        % 执行 PAM4 调制
        %txSig = pammod(symbol_indices, modu_order);
        txSig = pammod(symbol_indices,modu_order,0,'gray');
        % 添加 AWGN 噪声到 PAM4 符号
        noisy_symbols = awgn(txSig, snr_all(t, i), 'measured'); 

        % 硬判决解调（将符号映射回比特）
        demod_symbols = pamdemod(noisy_symbols, modu_order, 0, 'gray');

        % 将解调后的符号转换为比特流
        demod_data = de2bi(demod_symbols, log2(modu_order), 'left-msb')';
        demod_data = demod_data(:);

        % Viterbi 解码
        decode_matlab = vitdec(demod_data, trellis, tbdepth, 'cont', 'hard');
        decode_matlab_cleaned = decode_matlab(tbdepth * K + 1:end);
        
        % 修正长度不匹配问题
        dataIn_cut = dataIn(1:length(decode_matlab_cleaned))';

        % 确保 dataIn_cut 和 decode_matlab_cleaned 是列向量
        dataIn_cut = dataIn_cut(:);
        decode_matlab_cleaned = decode_matlab_cleaned(:);
        
        % 计算编码后的 BER
        berCoded(t, i) = biterr(dataIn_cut, decode_matlab_cleaned) / length(decode_matlab_cleaned);
        
        % 计算未编码的 BER（理论值）
        berUncoded(t, i) = berawgn(snr_all(t, i), 'pam', modu_order); % 对应的未编码 BER
    end
    
    % 结束计时并记录执行时间
    execution_times(t) = toc;
    
    % 显示当前架构的执行时间
    fprintf('Execution time for trellis structure %d: %.4f seconds\n', t, execution_times(t));
end

% 绘制图表
figure;
styles = {'b--o', 'r--s', 'g--d', 'k--x', 'm--^'};  % 增加样式
for t = 1:length(trellis_structs)
    % 绘制编码后的 BER 曲线
    semilogy(EbNo_all, berCoded(t,:), styles{t});
    hold on;
    % 绘制对应的未编码 BER 曲线
    semilogy(EbNo_all, berUncoded(t,:), [styles{t}(1) '-']);  % 使用相同颜色但不同线型
end
grid on;
xlabel('SNR/dB','FontSize', 20);
ylabel('BER','FontSize', 20);
set(gca, 'FontSize', 20); % 設置座標軸字體大小

legend('(2,1,3) Coded', '(2,1,3) Uncoded', '(2,1,7) Coded', '(2,1,7) Uncoded', '(3,2,4) Coded', '(3,2,4) Uncoded', '(4,3,4) Coded', '(4,3,4) Uncoded', '(4,3,5) Coded', '(4,3,5) Uncoded', 'FontSize', 20);
title('BER vs. SNR for Different Convolutional Codes + PAM4)');

% 显示所有架构的执行时间
for t = 1:length(trellis_structs)
    fprintf('Total execution time for trellis structure %d: %.4f seconds\n', t, execution_times(t));
end