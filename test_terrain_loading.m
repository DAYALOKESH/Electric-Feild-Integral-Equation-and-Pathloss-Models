% Simple test to verify terrain data loading
clear; clc;

fprintf('Testing terrain data loading...\n');

% Test loading X.04 data
try
    data = load('X.04');
    fprintf('Terrain data loaded successfully: %d points\n', size(data, 1));
    fprintf('Distance range: %.1f to %.1f m\n', min(data(:,1)), max(data(:,1)));
    fprintf('Height range: %.1f to %.1f m\n', min(data(:,2)), max(data(:,2)));
    
    % Test first 100m filtering
    max_dist = 100.0;
    valid_idx = data(:,1) <= max_dist;
    filtered_data = data(valid_idx, :);
    fprintf('First 100m data: %d points\n', size(filtered_data, 1));
    
catch ME
    fprintf('Error: %s\n', ME.message);
end

fprintf('Test completed.\n');