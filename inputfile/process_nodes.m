% MATLAB script to process, merge, and save node data

% --- 1. Read file ---
try
    fileContent = fileread('node.txt');
catch ME
    if strcmp(ME.identifier, 'MATLAB:fileread:cannotOpenFile')
        disp("Error: 'node.txt' not found. Ensure it's in the current directory or MATLAB path.");
        return;
    else
        rethrow(ME);
    end
end

% --- 2. Clean and parse data ---
cleanedContent = regexprep(fileContent, '\', '');
numbers = str2double(strsplit(cleanedContent)); 
numbers = numbers(~isnan(numbers)); 

if mod(length(numbers), 3) ~= 0
    disp('Warning: Total number of data points is not a multiple of 3. Trimming excess data.');
    numbers = numbers(1:end-mod(length(numbers), 3));
end

% --- 3. Reshape all data into a single matrix ---
numRows = length(numbers) / 3;
dataMatrix = reshape(numbers, 3, numRows)';

% --- 4. Find the delimiter and split the data ---
delimiterRow = [0 0 0];
[isRowPresent, rowIndex] = ismember(delimiterRow, dataMatrix, 'rows');

mainData = [];
secondBlockData = [];

if isRowPresent
    % First block is data before the delimiter
    mainData = dataMatrix(1:rowIndex-1, :);
    % Second block is data from the delimiter to the end
    secondBlockData = dataMatrix(rowIndex:end, :);
else
    % If no delimiter, treat all data as the first block
    mainData = dataMatrix;
    disp('Warning: Data delimiter [0 0 0] not found. No data to merge.');
end

% --- 5. Process the first data block ---
% Ensure we are only working with the first 400 lines
if size(mainData, 1) > 400
    mainData = mainData(1:400, :);
end
totalRows = size(mainData, 1);

% Create indices to extract 1st and 5th line of every 8-line block
indices_first = 1:8:totalRows;
indices_fifth = 5:8:totalRows;
final_indices = sort([indices_first, indices_fifth]);

% This is the "new matrix" from your request
processedFirstBlock = mainData(final_indices, :);

% --- 6. Merge the new matrix with the second block ---
finalMatrix = [processedFirstBlock; secondBlockData];

% --- 7. Write the final matrix to a new file ---
outputFileName = 'two_rod_nodes.txt';
writematrix(finalMatrix, outputFileName, 'Delimiter', ' ');

disp(['Operation complete. Merged data has been saved to "' outputFileName '".']);