close all; clear; clc;
load('data.mat');
load('awgn_data.mat');
load('uncoded_data.mat');

conv_code = conv_hardware_213(encoded_input);
EsN0_dB = 0:12;   
traceback = [2, 5, 12, 16, 1000]; 
num_tb = length(traceback);  
num_EsN0 = length(EsN0_dB);

BER_HARDV5 = zeros(num_tb, num_EsN0); 
BER_MATLAB = zeros(num_tb, num_EsN0); 
PAM4_CHANNEL = zeros(1,num_EsN0); 
trellis = poly2trellis(3, [5 7]); 

PAM4_ERRORS = zeros(1, num_EsN0); % 每個 SNR 錯誤的位數

for tb_idx = 1:num_tb
    k_tb = traceback(tb_idx);  
    
    for k_awgndata = 1:num_EsN0
        rxBits = awgn_data(1:end, k_awgndata);
        rxBits = transpose(rxBits);
        
        decoder_hradware_codeword = viterbi_hardwareV5_213(rxBits, k_tb);
        decoded_matlab = vitdec(rxBits,trellis,k_tb, 'trunc' , 'hard' );
        
        [~, ber_HARDV5] = biterr(encoded_input, decoder_hradware_codeword);
        BER_HARDV5(tb_idx, k_awgndata) = ber_HARDV5;

        [~, ber_MATLAB] = biterr(encoded_input, decoded_matlab);
        BER_MATLAB(tb_idx, k_awgndata) = ber_MATLAB;
        
    end
end

for k_awgndata = 1:num_EsN0
    rxBits = awgn_data(1:end, k_awgndata);
    rxBits = transpose(rxBits);
    [bit_errors, ber_PAM4_CHANNEL] = biterr(conv_code, rxBits);
    PAM4_CHANNEL(1,k_awgndata) = ber_PAM4_CHANNEL;
    PAM4_ERRORS(k_awgndata) = bit_errors;
    % 顯示錯誤數
    fprintf('SNR = %d dB: Bit Errors = %d\n', EsN0_dB(k_awgndata), bit_errors);
end
% 顯示錯誤數據表
disp('PAM4 Channel Bit Errors per SNR (dB):');
disp(array2table([EsN0_dB' PAM4_ERRORS'], 'VariableNames', {'SNR_dB', 'BitErrors'}));
%% **繪圖 (BER / SER)** 軟+硬
figure;
hold on; 
markers = {'-*', '-d', '-^', '-v', '-<', '->', '-p', '-h'}; % 8 個標記

% 繪製未編碼 (Uncoded) BER
semilogy(EsN0_dB, berUncoded, 'bo-', 'LineWidth', 1.5, 'DisplayName', 'Uncoded BER');
semilogy(EsN0_dB, PAM4_CHANNEL, 'go-', 'LineWidth', 1.5, 'DisplayName', 'Uncoded_BER_PAM4_CHANNEL');

for tb_idx = 1:num_tb
    % 繪製硬體 Viterbi BER
    semilogy(EsN0_dB, BER_HARDV5(tb_idx, :), ...
        markers{mod(tb_idx, length(markers)) + 1}, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Hardware Viterbi TB=%d', traceback(tb_idx)));

    % 繪製 MATLAB Viterbi BER
    semilogy(EsN0_dB, BER_MATLAB(tb_idx, :), ...
        markers{mod(tb_idx + 4, length(markers)) + 1}, 'LineWidth', 1.5, 'LineStyle', '--', ...
        'DisplayName', sprintf('MATLAB Viterbi TB=%d', traceback(tb_idx)));
end

grid on;
set(gca, 'XTick', EsN0_dB); 
xlabel('E_s/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
legend show; 
title('BER Performance Comparison (Viterbi & PAM4)');
hold off; 

N_total_assumed = 2000; % 假設總測試位數bits
BER_berUncoded = berUncoded;
N_errors_estimated = BER_berUncoded .* N_total_assumed; % 反算錯誤位數
disp('估計的錯誤位數:');
disp(N_errors_estimated);


%---------VITERBI_213_HARDWARE_V5-----------
function [decoded_msg] = viterbi_hardwareV5_213(conv_code,D)
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
                    new_path = [survivors{current_state+1}, input_bit];
                    if length(new_path) > D
                        new_path = new_path(end-D+1:end); 
                    end
                    new_survivors{next_state+1} = new_path;
                end
            end
        end

        path_metrics = new_metrics;
        survivors = new_survivors;

        if step >= D
            [~, best_state] = min(path_metrics);
            best_path = survivors{best_state};
            if length(best_path) >= D
                decoded_bit = best_path(1);
                decoded_msg = [decoded_msg, decoded_bit];

                for state = 1:4
                    if length(survivors{state}) >= 1
                        survivors{state} = survivors{state}(2:end);
                    end
                end
            end
        end
    end

    % 處理剩餘路徑
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