% MATLAB script to modify the first 392 lines of edge.txt

% --- 1. 生成第一部分内容 (替换原前392行) ---

% 生成序列 1: 1 3, 3 5, ..., 97 99
col1_seq1 = (1:2:97)';          % 第一列: 1, 3, 5, ..., 97
col2_seq1 = col1_seq1 + 2;      % 第二列: 3, 5, 7, ..., 99
new_block_1 = [col1_seq1, col2_seq1];

% 生成序列 2: 2 4, 4 6, ..., 98 100
col1_seq2 = (2:2:98)';          % 第一列: 2, 4, 6, ..., 98
col2_seq2 = col1_seq2 + 2;      % 第二列: 4, 6, 8, ..., 100
new_block_2 = [col1_seq2, col2_seq2];

% 合并第一部分内容
content_part1 = [new_block_1; new_block_2];


% --- 2. 生成第二部分内容 (替换原392行之后的内容) ---

% 生成序列 3: 101 102, 102 103, ..., 125 126
col1_seq3 = (101:125)';         % 第一列: 101, 102, ..., 125
col2_seq3 = col1_seq3 + 1;      % 第二列: 102, 103, ..., 126
content_part2 = [col1_seq3, col2_seq3];


% --- 3. 将所有新内容合并为一个矩阵 ---
final_matrix = [content_part1; content_part2];


% --- 5. Overwrite the original edge.txt file ---
output_filename = 'two_rod_edge.txt';
writematrix(final_matrix, output_filename, 'Delimiter', ' ');

disp(['Success! The file "' output_filename '" has been modified and overwritten.']);