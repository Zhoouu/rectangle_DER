function modifyNodes()
    % Define the name of the target file.
    fileName = 'two_rod_nodes.txt';

    % --- Step 1: Read the original data from the file ---
    % Use a try-catch block to handle potential file reading errors.
    try
        original_data = readmatrix(fileName);
    catch ME
        error('Failed to read the file: %s. Please ensure the file exists and is correctly formatted. MATLAB error: %s', fileName, ME.message);
    end

    % --- Step 2: Define parameters for the new semi-circular arc ---
    num_pairs = 50;      % The first 100 rows consist of 50 pairs of nodes.
    x_start = 0.0;       % Starting x-coordinate.
    x_end = 0.5;         % Ending x-coordinate.
    y_at_edge = 0.021;   % The y-magnitude at the start and end points.

    % --- Step 3: Generate the new coordinates ---
    % Create a linearly spaced vector for the 50 new x-coordinates.
    x_new = linspace(x_start, x_end, num_pairs)';

    % The semi-circle is centered on the x-axis at the midpoint of the x-range.
    center_x = (x_start + x_end) / 2; % This will be 0.25

    % Calculate the radius squared (R^2) of the circular arc.
    % The arc must pass through (x_start, y_at_edge).
    % The equation of the circle is (x - center_x)^2 + y^2 = R^2.
    radius_sq = (x_start - center_x)^2 + y_at_edge^2;

    % Calculate the new positive y-coordinates using the circle equation.
    % y = sqrt(R^2 - (x - center_x)^2)
    y_new_positive = sqrt(radius_sq - (x_new - center_x).^2);

    % --- Step 4: Assemble the new 100x3 data block ---
    new_section = zeros(100, 3);

    % Assign the new values, alternating between positive and negative y for symmetry.
    % Odd rows (1, 3, 5, ..., 99) get the positive y-values.
    new_section(1:2:end, 1) = x_new;
    new_section(1:2:end, 2) = y_new_positive;

    % Even rows (2, 4, 6, ..., 100) get the negative y-values.
    new_section(2:2:end, 1) = x_new;
    new_section(2:2:end, 2) = -y_new_positive;

    % The third column remains zero as per the original file format.
    new_section(:, 3) = 0;


    % --- Step 5: Replace the first 100 rows of the original data ---
    % Check if the original data has at least 100 rows.
    if size(original_data, 1) < 100
        error('The input file "%s" has fewer than 100 rows. Cannot perform the modification as requested.', fileName);
    end
    
    modified_data = original_data; % Make a copy to modify
    modified_data(1:100, :) = new_section;


    % --- Step 6: Write the modified data back to the file ---
    % This will overwrite the original file.
    try
        writematrix(modified_data, fileName, 'Delimiter', ' ');
        fprintf('Successfully modified the first 100 rows of "%s" and saved the changes.\n', fileName);
    catch ME
        error('Failed to write to the file: %s. Please check file permissions. MATLAB error: %s', fileName, ME.message);
    end
end