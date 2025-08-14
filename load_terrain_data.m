function [x_terrain, y_terrain, n_points] = load_terrain_data(filename, max_distance, delta_x)
%% Load and Process Terrain Data
% Loads terrain data from file and interpolates to required spatial resolution
% 
% Inputs:
%   filename - Name of terrain data file (X.04)
%   max_distance - Maximum distance to analyze (m)
%   delta_x - Spatial discretization step (m)
%
% Outputs:
%   x_terrain - Distance points for analysis (m)
%   y_terrain - Interpolated terrain heights (m)
%   n_points - Number of analysis points

try
    % Load terrain data using MATLAB's load function
    terrain_data = load(filename);
    x_raw = terrain_data(:, 1);  % Distance values
    y_raw = terrain_data(:, 2);  % Height values
    
    fprintf('Original terrain data: %d points\n', length(x_raw));
    fprintf('Distance range: %.1f to %.1f m\n', min(x_raw), max(x_raw));
    fprintf('Height range: %.1f to %.1f m\n', min(y_raw), max(y_raw));
    
    % Filter data for first 100 meters (testing requirement)
    valid_indices = x_raw <= max_distance;
    x_filtered = x_raw(valid_indices);
    y_filtered = y_raw(valid_indices);
    
    % Create uniform grid for analysis based on electromagnetic wavelength
    x_terrain = 0:delta_x:max_distance;
    n_points = length(x_terrain);
    
    % Interpolate terrain heights to uniform grid
    % Handle edge cases for extrapolation
    if length(x_filtered) < 2
        error('Insufficient terrain data for analysis range');
    end
    
    % Use linear interpolation with extrapolation
    y_terrain = interp1(x_filtered, y_filtered, x_terrain, 'linear', 'extrap');
    
    % Validate interpolated data
    if any(isnan(y_terrain))
        error('Invalid terrain interpolation - check data quality');
    end
    
    fprintf('Analysis grid: %d points from 0 to %.1f m\n', n_points, max_distance);
    fprintf('Grid spacing: %.3f m (λ/4)\n', delta_x);
    
catch ME
    error('Failed to load terrain data from %s: %s', filename, ME.message);
end

end