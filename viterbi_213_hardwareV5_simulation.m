close all; clear; clc;
%% 參數設定
M = 4;                  % 4-PAM
k = log2(M);            % 每個符號所含的 bit 數 (4-PAM -> 2 bits/symbol)
EsN0_dB = -12:2:16;       % 模擬的 Es/N0 (dB) 範圍

trellis = poly2trellis(3, [5 7]);
traceback = [2, 5, 12, 16, 1000]; 

%% 產生隨機位元
txBits = load('dataIn.asv', '-ascii');
txBits = txBits(1:1000);
txBits = transpose(txBits);

%% 將 bits 映射 (2 bits -> 1 symbol index)
conv_data = conv_hardware_213(txBits);
symIdxTx_data = bi2de(reshape(conv_data, k, []).','left-msb');
symIdxTx_Uncoded_data = bi2de(reshape(txBits, k, []).','left-msb');
%% 4-PAM 調變 (帶灰度編碼)
% phaseOffset=0, 'gray' => 使用 Gray Mapping
txSymbols = pammod(symIdxTx_data, M, 0, 'gray');    % 回傳的調變結果為 4 個星座點（實數）
txSymbols_Uncoded = pammod(symIdxTx_Uncoded_data, M, 0, 'gray');    % 回傳的調變結果為 4 個星座點（實數）

%% 理論計算用：先量測已成形波形的「符號能量(平均)」
% 在「基帶實數 PAM」情況下，可估計每符號平均能量。
% E_s = mean(abs(txSymbols).^2)；對 4-PAM 的預設星座點，也可手動帶入
Es = mean(abs(txSymbols).^2);

num_tb = length(traceback);  
num_EsN0 = length(EsN0_dB);

BER_HARDV5 = zeros(num_tb, num_EsN0); 
BER_MATLAB = zeros(num_tb, num_EsN0); 
PAM4_CHANNEL = zeros(1,num_EsN0); 


PAM4_ERRORS = zeros(1, num_EsN0); % 每個 SNR 錯誤的位數


%計算Uncoded ber
for ii = 1:length(EsN0_dB)
    thisEsN0_dB = EsN0_dB(ii);
    thisEsN0    = 10^(thisEsN0_dB/10);      % 線性值
    sigma = sqrt( Es / thisEsN0 );

    % 產生 AWGN 並加到已成形訊號 txWave
    noise = sigma * randn(size(txSymbols_Uncoded));
    rxSymbolBlock_Uncoded = txSymbols_Uncoded+noise;
    % 解調 (帶灰度解調)
    symIdxRx_Uncoded = pamdemod(rxSymbolBlock_Uncoded, M, 0, 'gray');

    [bit_Uncoded_errors, ber_Uncoded] = biterr(symIdxTx_Uncoded_data, symIdxRx_Uncoded);
    BER_Uncoded(ii) = ber_Uncoded;
    BER_Uncoded_errors(ii) = bit_Uncoded_errors;
end

%% 算一次BER
% for tb_idx = 1:num_tb
%     k_tb = traceback(tb_idx);  
% 
%     for ii = 1:length(EsN0_dB)
%         thisEsN0_dB = EsN0_dB(ii);
%         thisEsN0    = 10^(thisEsN0_dB/10);      % 線性值
%         % thisEsN0  = EsN0_dB(ii);
%         % -----------------------------
%         % *對"實數基帶 PAM"而言*
%         % SNR = Es / (N0/2) ==> N0 = 2 * Es / SNR
%         % 噪聲(單邊)方差 sigma^2 = N0/2 (因 randn() 預設功率=1)
%         % => sigma^2 = (2 * Es / SNR)/2 = Es / SNR
%         % => sigma = sqrt( Es / SNR )
%         % -----------------------------
%         sigma = sqrt( Es / thisEsN0 );
%         % 產生 AWGN 並加到已成形訊號 txWave
% 
%         noise = sigma * randn(size(txSymbols));
% 
%         rxSymbolBlock = txSymbols + noise;
% 
%         % 解調 (帶灰度解調)
%         symIdxRx = pamdemod(rxSymbolBlock, M, 0, 'gray');
%         recovered_symIdxRx_data = de2bi(symIdxRx, k, 'left-msb');
%         % 轉換為 1D 向量 (column-wise 展開)
%         recovered_bits = reshape(recovered_symIdxRx_data.',1,[]);
% 
%         % decoder_hradware_codeword = viterbi_hardwareV5_213(rxBits, k_tb);
%         decoded_matlab = vitdec(recovered_bits,trellis,k_tb, 'term' , 'hard' );
%         % [~, ber_HARDV5] = biterr(encoded_input, decoder_hradware_codeword);
%         % BER_HARDV5(tb_idx, k_awgndata) = ber_HARDV5;
% 
%         [bit_errors, ber_PAM4_CHANNEL] = biterr(conv_data, recovered_bits);
%         PAM4_CHANNEL(1,ii) = ber_PAM4_CHANNEL;
%         PAM4_ERRORS(ii) = bit_errors;
% 
%         % [ber_MATLAB_errors, ber_MATLAB] = biterr(symIdxTx_data, symIdxRx);
%         % BER_MATLAB(tb_idx, ii) = ber_MATLAB;
%         % BER_MATLAB_ERROR(tb_idx, ii) = ber_MATLAB_errors;
% 
%         [ber_MATLAB_errors, ber_MATLAB] = biterr(txBits, decoded_matlab);
%         BER_MATLAB(tb_idx, ii) = ber_MATLAB;
%         BER_MATLAB_ERROR(tb_idx, ii) = ber_MATLAB_errors;
%     end
% end

%% 每個SNR 測試10次平均
num_trials = 10; % 每個 SNR 測試 10 次

% 預先分配變數
BER_PAM4_CHANNEL = zeros(1, length(EsN0_dB));
PAM4_ERRORS = zeros(1, length(EsN0_dB));
BER_MATLAB = zeros(num_tb, length(EsN0_dB));
BER_MATLAB_ERROR = zeros(num_tb, length(EsN0_dB));
BER_HARDWARE_BEHAVIOR = zeros(num_tb, length(EsN0_dB));
BER_HARDWARE_BEHAVIOR_ERROR = zeros(num_tb, length(EsN0_dB));
for tb_idx = 1:num_tb
    k_tb = traceback(tb_idx);  

    for ii = 1:length(EsN0_dB)
        thisEsN0_dB = EsN0_dB(ii);
        thisEsN0 = 10^(thisEsN0_dB / 10); % 線性值
        sigma = sqrt(Es / thisEsN0); % 計算噪聲標準差

        % 初始化累加變數
        pam4_ber_sum = 0;
        pam4_errors_sum = 0;
        matlab_ber_sum = 0;
        matlab_errors_sum = 0;
        hardware_ber_sum = 0;
        hardware_errors_sum = 0;
        for trial = 1:num_trials
            % 產生 AWGN 並加到已成形訊號 txWave
            noise = sigma * randn(size(txSymbols));
            rxSymbolBlock = txSymbols + noise;

            % 解調 (帶灰度解調)
            symIdxRx = pamdemod(rxSymbolBlock, M, 0, 'gray');
            recovered_symIdxRx_data = de2bi(symIdxRx, k, 'left-msb');
            recovered_bits = reshape(recovered_symIdxRx_data.', 1, []);

            % MATLAB Viterbi 解碼
            decoded_matlab = vitdec(recovered_bits, trellis, k_tb, 'term', 'hard');
            decoded_hardware_behavior = viterbi_hardwareV5(recovered_bits,k_tb);
            
            % PAM4 Channel BER 計算
            [bit_errors, ber_PAM4_CHANNEL_trial] = biterr(conv_data, recovered_bits);
            pam4_ber_sum = pam4_ber_sum + ber_PAM4_CHANNEL_trial;
            pam4_errors_sum = pam4_errors_sum + bit_errors;

            % MATLAB Viterbi BER 計算
            [ber_MATLAB_errors, ber_MATLAB_trial] = biterr(txBits, decoded_matlab);
            matlab_ber_sum = matlab_ber_sum + ber_MATLAB_trial;
            matlab_errors_sum = matlab_errors_sum + ber_MATLAB_errors;

            % MATLAB Hardware Behavior Viterbi BER 計算
            [ber_hardware_errors, ber_hardware_trial] = biterr(txBits, decoded_hardware_behavior);
            hardware_ber_sum = hardware_ber_sum + ber_hardware_trial;
            hardware_errors_sum = hardware_errors_sum + ber_hardware_errors;
        end

        % 計算 10 次的平均 BER 和錯誤數
        BER_PAM4_CHANNEL(1, ii) = pam4_ber_sum / num_trials;
        PAM4_ERRORS(ii) = pam4_errors_sum / num_trials;
        BER_MATLAB(tb_idx, ii) = matlab_ber_sum / num_trials;
        BER_MATLAB_ERROR(tb_idx, ii) = matlab_errors_sum / num_trials;
        BER_HARDWARE_BEHAVIOR(tb_idx, ii)  = hardware_ber_sum / num_trials;
        BER_HARDWARE_BEHAVIOR_ERROR(tb_idx, ii) = hardware_errors_sum / num_trials;
    end
end


%% **繪圖 (BER / SER)** 軟+硬
figure;
hold on; 
markers = {'-*', '-d', '-^', '-v', '-<', '->', '-p', '-h'}; % 8 個標記

% 繪製未編碼 (Uncoded) BER
semilogy(EsN0_dB, BER_Uncoded, 'bo-', 'LineWidth', 1.5, 'DisplayName', 'Uncoded BER');
% semilogy(EsN0_dB, PAM4_CHANNEL, 'go-', 'LineWidth', 1.5, 'DisplayName', 'Uncoded_BER_PAM4_CHANNEL');

for tb_idx = 1:num_tb
    % 繪製硬體 Viterbi BER
    semilogy(EsN0_dB, BER_HARDWARE_BEHAVIOR(tb_idx, :), ...
        markers{mod(tb_idx, length(markers)) + 1}, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Hardware behavior Viterbi TB=%d', traceback(tb_idx)));

    % 繪製 MATLAB Viterbi BER
    % semilogy(EsN0_dB, BER_MATLAB(tb_idx, :), ...
    %     markers{mod(tb_idx + 4, length(markers)) + 1}, 'LineWidth', 1.5, 'LineStyle', '--', ...
    %     'DisplayName', sprintf('MATLAB Viterbi TB=%d', traceback(tb_idx)));
end

grid on;
set(gca, 'XTick', EsN0_dB); 
xlabel('E_s/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
legend show; 
title('BER Performance Comparison (Viterbi & PAM4)');
hold off; 

%---------VITERBI_213_HARDWARE_V5-----------
function [decoded_msg] = viterbi_hardwareV5(conv_code,D)
    survivors = cell(4,1);       % 儲存各狀態的存活路徑
    new_survivors = cell(4,1);   % 暫存新生成的存活路徑
    path_metrics = [0;3;3;3];    % 各狀態的路徑度量值
    decoded_msg = [];            % 儲存逐步輸出的解碼結果
    
    for step = 1:length(conv_code)/2
        idx = 2*step - 1;
        received_bits = conv_code(idx:idx+1);
        new_metrics = inf(4,1);    % 初始化新度量為無窮大
        new_survivors = cell(4,1); % 重置暫存器
        
        for current_state = 0:3
            for input_bit = 0:1
    
                next_state = viterbi_next_state(current_state, input_bit);
                output_dec = viterbi_outputs(current_state, input_bit);
                expected_bits = de2bi(output_dec, 2, 'left-msb');
                
                hamming_dist = sum(received_bits ~= expected_bits);
                
                candidate_metric = path_metrics(current_state+1) + hamming_dist;
                
                if candidate_metric < new_metrics(next_state+1)
                    new_metrics(next_state+1) = candidate_metric;
                    new_survivors{next_state+1} = [survivors{current_state+1}, input_bit];
                end
            end
        end
        
        path_metrics = new_metrics;
        survivors = new_survivors;
    
        row_lengths = cellfun(@length, survivors);
        if any(row_lengths >= D)
            [~, best_state] = min(path_metrics);
            best_path = survivors{best_state};
            decoded_msg = [decoded_msg, best_path];
            survivors = cell(4,1);
        end
    end
    
    [~, final_state] = min(path_metrics);
    remaining_bits = survivors{final_state};
    decoded_msg = [decoded_msg, remaining_bits];
end
%---------VITERBI_HARDWARE_TABLE-----------
function nextState = viterbi_next_state(currentState,inptBits)
    switch currentState
        case 0 
            if(inptBits == 0)
                nextState = 0;
            else
                nextState = 2;
            end
         case 1 
            if(inptBits == 0)
                nextState = 0;
            else
                nextState = 2;
            end
         case 2 
            if(inptBits == 0)
                nextState = 1;
            else
                nextState = 3;
            end
         case 3 
            if(inptBits == 0)
                nextState = 1;
            else
                nextState = 3;
            end
    end
end

function outputs = viterbi_outputs(currentState,inptBits)
    switch currentState
        case 0 
            if(inptBits == 0)
                outputs = 0;
            else
                outputs = 3;
            end
         case 1 
            if(inptBits == 0)
                outputs = 3;
            else
                outputs = 0;
            end
         case 2 
            if(inptBits == 0)
                outputs = 1;
            else
                outputs = 2;
            end
         case 3 
            if(inptBits == 0)
                outputs = 2;
            else
                outputs = 1;
            end
    end
end

%---------hardware_conv213_function-----------
function codeword = conv_hardware_213(msg_source)
    s1 = 0;
    s2 = 0;
    bit_string_length  = length(msg_source);
    codeword = zeros(1, bit_string_length * 2);
    for i = 1:bit_string_length
        u0 = xor(msg_source(i), s2);
        u1 = xor(xor(msg_source(i), s1), s2);
        s2 = s1;
        s1 = msg_source(i);
        codeword(2*i-1) = u0;
        codeword(2*i) = u1;
    end
end
