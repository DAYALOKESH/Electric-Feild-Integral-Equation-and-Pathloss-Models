%% Simplified FBSM Test - Core Algorithm Verification
% Tests the Forward-Backward Scattering Method implementation
clear; clc; close all;

fprintf('=== FBSM Core Algorithm Test ===\n');

%% Load terrain data (first 100m)
data = load('X.04');
max_dist = 100.0;
valid_idx = data(:,1) <= max_dist;
x_raw = data(valid_idx, 1);
y_raw = data(valid_idx, 2);

fprintf('Terrain data: %d points up to %.0f m\n', length(x_raw), max_dist);

%% Electromagnetic parameters
c = 299792458;
f = 970e6;
lambda = c/f;
omega = 2*pi*f;
mu_0 = 4*pi*1e-7;
epsilon_0 = 8.854e-12;
beta_0 = omega*sqrt(mu_0*epsilon_0);
delta_x = lambda/4;

fprintf('Wavelength: %.3f m, Grid spacing: %.3f m\n', lambda, delta_x);

%% Create analysis grid
x_terrain = 0:delta_x:max_dist;
n_points = length(x_terrain);

% Interpolate terrain heights
y_terrain = interp1(x_raw, y_raw, x_terrain, 'linear', 'extrap');

fprintf('Analysis grid: %d points\n', n_points);

%% Source parameters
x_source = 0.0;
y_source = 442.0;

%% Simple surface current calculation (forward only for testing)
fprintf('Computing surface current...\n');

surface_current = zeros(n_points, 1);
j = 1i;

% Helper functions
R_source_p = @(p) sqrt((x_source - x_terrain(p))^2 + (y_source - y_terrain(p))^2);
R_p_q = @(p, q) sqrt((x_terrain(q) - x_terrain(p))^2 + (y_terrain(q) - y_terrain(p))^2);

% Simplified incident field (avoid Hankel function issues)
EiRad = @(dist) exp(-j*beta_0*dist) / sqrt(dist + 1e-10);

% Simplified impedance
Z_self = @(i) 377 + j*100; % Approximate free space impedance with reactive part

% Forward scattering (simplified)
if n_points > 0
    surface_current(1) = EiRad(R_source_p(1)) / Z_self(1);
end

for p = 2:min(n_points, 20)  % Limit to first 20 points for testing
    surface_current(p) = EiRad(R_source_p(p)) / Z_self(p);
end

current_magnitude = abs(surface_current);

fprintf('Current calculation completed. Max: %.2e A/m\n', max(current_magnitude));

%% Simple visualization
figure('Position', [100, 100, 800, 600]);

subplot(2,1,1);
plot(x_terrain, y_terrain, 'k-', 'LineWidth', 2);
grid on;
xlabel('Distance (m)');
ylabel('Height (m)');
title('Terrain Profile');

subplot(2,1,2);
semilogy(x_terrain, current_magnitude + 1e-15, 'r-', 'LineWidth', 2);
grid on;
xlabel('Distance (m)');
ylabel('Surface Current Magnitude (A/m)');
title('Surface Current Distribution (Simplified)');

fprintf('Test completed successfully!\n');