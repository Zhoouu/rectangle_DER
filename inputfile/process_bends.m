% MATLAB script to completely regenerate the content for bend.txt

% --- 1. Define the sequences for the first column ---
% This creates the three separate ranges you specified.
col1_part1 = (1:48)';
col1_part2 = (50:97)';
col1_part3 = (99:122)';

% --- 2. Combine the sequences into a single first column ---
first_col = [col1_part1; col1_part2; col1_part3];

% --- 3. Generate the second column ---
% Following the pattern in the original file, we assume the second 
% column contains the adjacent number (first column + 1).
second_col = first_col + 1;

% --- 4. Create the final two-column matrix ---
final_matrix = [first_col, second_col];

% --- 5. Write the final matrix to bend.txt, overwriting it ---
output_filename = 'two_rod_bends.txt';
try
    writematrix(final_matrix, output_filename, 'Delimiter', ' ');
    disp(['Success! The file "' output_filename '" has been completely overwritten with the new data.']);
catch ME
    disp(['Error: Could not write to file "' output_filename '".']);
    disp('Please check file permissions or if the file is open in another program.');
    rethrow(ME);
end