function [terrain_distance, terrain_height] = load_terrain(filename)
    % Load terrain data from tab/space-separated file
    % Input: filename - name of terrain data file (e.g., 'X.04')
    % Output: terrain_distance - array of distance values in meters
    %         terrain_height - array of height values in meters
    %
    % The terrain file should contain distance-height pairs separated by
    % whitespace (tabs or spaces)
    
    try
        % Check if file exists
        if ~exist(filename, 'file')
            error('Terrain file %s not found', filename);
        end
        
        % Load the terrain data using MATLAB's load function
        % This handles tab-separated or space-separated data automatically
        terrain_data = load(filename);
        
        % Extract distance and height columns
        if size(terrain_data, 2) < 2
            error('Terrain file must contain at least 2 columns (distance and height)');
        end
        
        terrain_distance = terrain_data(:, 1);  % First column: distance (m)
        terrain_height = terrain_data(:, 2);    % Second column: height (m)
        
        % Validate the data
        if length(terrain_distance) ~= length(terrain_height)
            error('Distance and height arrays must have the same length');
        end
        
        if any(isnan(terrain_distance)) || any(isnan(terrain_height))
            warning('NaN values found in terrain data, removing affected points');
            valid_indices = ~(isnan(terrain_distance) | isnan(terrain_height));
            terrain_distance = terrain_distance(valid_indices);
            terrain_height = terrain_height(valid_indices);
        end
        
        if length(terrain_distance) < 2
            error('Insufficient valid terrain data points');
        end
        
        % Sort by distance to ensure monotonic ordering
        [terrain_distance, sort_idx] = sort(terrain_distance);
        terrain_height = terrain_height(sort_idx);
        
        % Display terrain data summary
        fprintf('Terrain data loaded successfully:\n');
        fprintf('  Number of points: %d\n', length(terrain_distance));
        fprintf('  Distance range: %.2f to %.2f meters\n', min(terrain_distance), max(terrain_distance));
        fprintf('  Height range: %.2f to %.2f meters\n', min(terrain_height), max(terrain_height));
        
    catch ME
        error('Failed to load terrain data from %s: %s', filename, ME.message);
    end
end