% Forward and Backward Scattering Method (FBSM) for Electromagnetic Scattering Analysis
% Implementation based on Electric Field Integral Equation (EFIE) method
% Reference: efie.txt C++ implementation (Balanis page 680)
%
% This script implements electromagnetic scattering analysis on terrain data
% using the Forward and Backward Scattering Method with EFIE formulation.
%
% Author: MATLAB FBSM Implementation
% Date: 2024

function fbsm_main()
    % Clear workspace and set up environment
    clear all;
    close all;
    clc;
    
    fprintf('=== Forward and Backward Scattering Method (FBSM) ===\n');
    fprintf('Loading terrain data and initializing parameters...\n');
    
    try
        % Load terrain data from X.04 file
        [terrain_distance, terrain_height] = load_terrain('X.04');
        fprintf('Successfully loaded terrain data: %d points\n', length(terrain_distance));
        
        % Set electromagnetic parameters (from efie.txt reference)
        params = set_em_parameters();
        fprintf('Electromagnetic parameters initialized\n');
        
        % Calculate surface current using EFIE with forward/backward scattering
        fprintf('Computing surface current distribution...\n');
        [surface_current, current_positions] = calculate_scattering(terrain_distance, terrain_height, params);
        
        % Calculate electric field distribution
        fprintf('Computing electric field distribution...\n');
        [electric_field, field_positions] = calculate_electric_field(surface_current, current_positions, terrain_distance, terrain_height, params);
        
        % Create visualizations
        fprintf('Generating visualizations...\n');
        create_visualizations(terrain_distance, terrain_height, surface_current, current_positions, electric_field, field_positions, params);
        
        fprintf('FBSM analysis completed successfully!\n');
        
    catch ME
        fprintf('Error in FBSM analysis: %s\n', ME.message);
        fprintf('Error occurred in: %s\n', ME.stack(1).name);
        rethrow(ME);
    end
end

function params = set_em_parameters()
    % Set electromagnetic parameters based on efie.txt reference
    params = struct();
    
    % Physical constants
    params.c = 299792458;                    % Speed of light (m/s)
    params.epsilon_0 = 8.854e-12;           % Free space permittivity
    params.mu_0 = 4*pi*1e-7;                % Free space permeability
    params.eta_0 = sqrt(params.mu_0/params.epsilon_0); % Free space impedance
    
    % Operating parameters (from efie.txt)
    params.frequency = 970e6;                % Frequency: 970 MHz
    params.omega = 2*pi*params.frequency;    % Angular frequency
    params.lambda = params.c/params.frequency; % Wavelength
    params.beta_0 = params.omega*sqrt(params.mu_0*params.epsilon_0); % Phase constant
    
    % Discretization parameters
    params.delta_x = params.lambda/4.0;     % Spatial discretization (lambda/4)
    params.gross_step = 10.0;               % Gross step from terrain data
    
    % Source parameters (from efie.txt)
    params.x_source = 0.0;                  % Source x-position
    params.y_source = 442.0;                % Source y-position
    params.current_amplitude = 1.0;         % Source current amplitude
    params.field_amplitude = 1.0;           % Electric field amplitude
    
    % Observation height offset
    params.obs_height_offset = 2.4;         % Observation point height offset
    
    % Tolerance for numerical calculations
    params.tolerance = 1e-15;
end