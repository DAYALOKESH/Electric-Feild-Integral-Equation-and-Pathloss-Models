function [electric_field, field_magnitude] = calculate_electric_field(x_terrain, y_terrain, ...
    surface_current, x_source, y_source, beta_0, omega, epsilon_0, delta_x, n_points)
%% Calculate Electric Field Distribution
% Computes total electric field using surface current and incident field
% Based on equations 7 & 8 from efie.txt reference
%
% Inputs:
%   x_terrain, y_terrain - Terrain geometry
%   surface_current - Surface current distribution from FBSM
%   x_source, y_source - Source position
%   beta_0 - Wave number
%   omega - Angular frequency
%   epsilon_0 - Permittivity
%   delta_x - Spatial step
%   n_points - Number of analysis points
%
% Outputs:
%   electric_field - Complex electric field distribution
%   field_magnitude - Magnitude of electric field

% Initialize electric field array
electric_field = zeros(n_points, 1);

fprintf('Computing electric field distribution...\n');

%% Helper Functions (from efie.txt reference)

% Distance from source to observation point (with 2.4m offset as in C++ code)
R_source_obs = @(p) sqrt((x_source - x_terrain(p))^2 + (y_source - y_terrain(p) - 2.4)^2);

% Distance from surface point to observation point (with 2.4m offset)  
R_surf_obs = @(p, q) sqrt((x_terrain(q) - x_terrain(p))^2 + ((y_terrain(q) + 2.4) - y_terrain(p))^2);

% Distance between terrain points
R_p_q = @(p, q) sqrt((x_terrain(q) - x_terrain(p))^2 + (y_terrain(q) - y_terrain(p))^2);

% Incident field radiation
EiRad = @(dist) -((beta_0^2)/(4*omega*epsilon_0)) * hankel_h2(beta_0 * dist);

% Impedance function for field calculation
Z_field = @(R) ((beta_0^2)/(4*omega*epsilon_0)) * hankel_h2(beta_0 * R);

%% Electric Field Calculation (Equations 7 & 8 from efie.txt)
% Et[index] = EiRad(R_source_obs(index)) - sum(J[n] * R_p_q(n,n+1) * Z(R_surf_obs(n,index)))

for index = 1:n_points
    % Initialize scattered field contribution
    scattered_field = 0;
    
    % Sum contributions from all surface current elements (equation 7 & 8)
    for n = 1:index
        % Segment length for current element
        if n < n_points
            segment_length = R_p_q(n, n+1);
        else
            segment_length = delta_x; % Use grid spacing for last element
        end
        
        % Scattered field contribution from surface current element n
        R_obs = R_surf_obs(n, index);
        scattered_field = scattered_field + surface_current(n) * segment_length * Z_field(R_obs);
    end
    
    % Total field = Incident field - Scattered field (equation from C++ code)
    incident_field = EiRad(R_source_obs(index));
    electric_field(index) = incident_field - scattered_field;
    
    % Progress indicator
    if mod(index, max(1, floor(n_points/10))) == 0
        fprintf('Electric field calculation: %d%% complete\n', round(100*index/n_points));
    end
end

%% Calculate field magnitude
field_magnitude = abs(electric_field);

% Apply distance normalization as in C++ code (20*log10(|E|/sqrt(R_source_obs)))
% This gives field strength in dB relative to 1 V/m normalized by distance
field_magnitude_db = zeros(n_points, 1);
for index = 1:n_points
    R_norm = sqrt(R_source_obs(index));
    if R_norm > 0
        field_magnitude_db(index) = 20 * log10(field_magnitude(index) / R_norm);
    else
        field_magnitude_db(index) = -inf; % Handle division by zero
    end
end

fprintf('Electric field range: %.2e to %.2e V/m\n', min(field_magnitude), max(field_magnitude));
fprintf('Normalized field range: %.1f to %.1f dB\n', min(field_magnitude_db(~isinf(field_magnitude_db))), ...
    max(field_magnitude_db(~isinf(field_magnitude_db))));

% Store both linear and dB values for visualization
field_magnitude = struct('linear', field_magnitude, 'dB', field_magnitude_db);

end

function h = hankel_h2(z)
%% Second kind Hankel function of order 0
% H0^(2)(z) = J0(z) - i*Y0(z)
% Using MATLAB's besselh function

h = besselh(0, 2, z);

% Handle numerical issues for small arguments
small_arg_idx = abs(z) < 1e-10;
if any(small_arg_idx)
    % For very small arguments, use series expansion to avoid numerical issues
    h(small_arg_idx) = 1 - 2i/pi * (log(z(small_arg_idx)/2) + 0.5772156649015329); % Euler's gamma
end

end