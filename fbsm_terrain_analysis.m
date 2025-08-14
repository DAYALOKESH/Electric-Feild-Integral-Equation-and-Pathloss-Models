%% Forward and Backward Scattering Method for Electromagnetic Terrain Analysis
% Implementation based on Electric Field Integral Equation (EFIE)
% Reference: Balanis "Advanced Engineering Electromagnetics" Page 680
% Analyzes electromagnetic scattering from terrain surface using X.04 data

clear; clc; close all;

%% Electromagnetic Constants and Parameters (from efie.txt reference)
c = 299792458;                    % Speed of light (m/s)
f = 970e6;                        % Frequency (Hz)
lambda = c/f;                     % Wavelength (m)
omega = 2*pi*f;                   % Angular frequency (rad/s)
mu_0 = 4*pi*1e-7;                 % Permeability of free space (H/m)
epsilon_0 = 8.854e-12;            % Permittivity of free space (F/m)
eta_0 = sqrt(mu_0/epsilon_0);     % Impedance of free space (ohm)
beta_0 = omega*sqrt(mu_0*epsilon_0); % Wave number (rad/m)

% Discretization parameters
delta_x = lambda/4;               % Spatial discretization step (m)
gross_step = 10.0;                % Original terrain sampling step (m)

% Source parameters
x_source = 0.0;                   % Source x-position (m)
y_source = 442.0;                 % Source y-position (m) 
I_source = 1.0;                   % Source current amplitude (A)

% Test range: first 100 meters only as specified
max_distance = 100.0;             % Maximum distance for analysis (m)

fprintf('=== Forward and Backward Scattering Method ===\n');
fprintf('Frequency: %.1f MHz\n', f/1e6);
fprintf('Wavelength: %.3f m\n', lambda);
fprintf('Analysis range: 0 to %.0f meters\n', max_distance);

%% Load and Process Terrain Data
fprintf('\nLoading terrain data from X.04...\n');
[x_terrain, y_terrain, n_points] = load_terrain_data('X.04', max_distance, delta_x);

fprintf('Terrain points loaded: %d\n', n_points);
fprintf('Spatial resolution: %.3f m\n', delta_x);

%% Calculate Surface Current using Forward-Backward Scattering Method
fprintf('\nCalculating surface current distribution...\n');
[surface_current, current_magnitude] = calculate_surface_current(x_terrain, y_terrain, ...
    x_source, y_source, beta_0, omega, epsilon_0, mu_0, delta_x, n_points);

fprintf('Surface current calculation completed.\n');

%% Calculate Electric Field Distribution  
fprintf('\nCalculating electric field distribution...\n');
[electric_field, field_magnitude] = calculate_electric_field(x_terrain, y_terrain, ...
    surface_current, x_source, y_source, beta_0, omega, epsilon_0, delta_x, n_points);

fprintf('Electric field calculation completed.\n');

%% Visualize Results
fprintf('\nGenerating visualizations...\n');
visualize_results(x_terrain, y_terrain, current_magnitude, field_magnitude, lambda);

fprintf('\nAnalysis completed successfully!\n');