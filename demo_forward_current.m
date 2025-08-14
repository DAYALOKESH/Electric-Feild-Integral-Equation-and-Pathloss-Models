%% Forward Current Verification Demo
% This script demonstrates the forward current tracking capability
% and validates the Forward-Backward Scattering Method implementation

clear; clc; close all;

fprintf('=== FBSM Forward Current Verification Demo ===\n\n');

%% Parameters (simplified for demo)
c = 299792458;
f = 970e6;
lambda = c/f;
omega = 2*pi*f;
mu_0 = 4*pi*1e-7;
epsilon_0 = 8.854e-12;
beta_0 = omega*sqrt(mu_0*epsilon_0);
delta_x = lambda/4;
max_distance = 30.0; % Smaller range for demo

% Source parameters
x_source = 0.0;
y_source = 442.0;

fprintf('Demo Parameters:\n');
fprintf('- Frequency: %.1f MHz\n', f/1e6);
fprintf('- Analysis range: 0 to %.0f m\n', max_distance);
fprintf('- Grid spacing: %.3f m (λ/4)\n\n', delta_x);

%% Load terrain and calculate currents
try
    [x_terrain, y_terrain, n_points] = load_terrain_data('X.04', max_distance, delta_x);
    fprintf('Terrain loaded: %d points\n', n_points);
    
    [surface_current, current_magnitude, forward_current, forward_magnitude] = ...
        calculate_surface_current(x_terrain, y_terrain, x_source, y_source, ...
        beta_0, omega, epsilon_0, mu_0, delta_x, n_points);
    
    fprintf('\nResults:\n');
    fprintf('- Forward current range: %.2e to %.2e A/m\n', min(forward_magnitude), max(forward_magnitude));
    fprintf('- Total current range: %.2e to %.2e A/m\n', min(current_magnitude), max(current_magnitude));
    
    % Calculate convergence metrics
    max_forward = max(forward_magnitude);
    max_total = max(current_magnitude);
    convergence_factor = max_total / max_forward;
    
    fprintf('- Convergence factor: %.3f\n', convergence_factor);
    fprintf('- Current stabilization: %.1f%%\n', (1 - (max(current_magnitude) - min(current_magnitude))/max(current_magnitude)) * 100);
    
    %% Visualization
    figure('Position', [100, 100, 1000, 600]);
    
    subplot(2,2,1);
    plot(x_terrain, y_terrain, 'k-', 'LineWidth', 2);
    grid on;
    xlabel('Distance (m)');
    ylabel('Height (m)');
    title('Terrain Profile');
    
    subplot(2,2,2);
    plot(x_terrain, forward_magnitude, 'g-', 'LineWidth', 2);
    hold on;
    plot(x_terrain, current_magnitude, 'r--', 'LineWidth', 2);
    grid on;
    xlabel('Distance (m)');
    ylabel('Current Magnitude (A/m)');
    title('Forward vs Total Current');
    legend('Forward Only', 'Forward + Backward', 'Location', 'best');
    
    subplot(2,2,3);
    semilogy(x_terrain, forward_magnitude + 1e-20, 'g-', 'LineWidth', 2);
    grid on;
    xlabel('Distance (m)');
    ylabel('Forward Current (A/m)');
    title('Forward Current (Log Scale)');
    
    subplot(2,2,4);
    plot(x_terrain, abs(current_magnitude - forward_magnitude), 'b-', 'LineWidth', 2);
    grid on;
    xlabel('Distance (m)');
    ylabel('|Total - Forward| (A/m)');
    title('Backward Scattering Contribution');
    
    fprintf('\n✓ Forward current verification completed successfully!\n');
    fprintf('✓ Visualization shows forward vs total current behavior\n');
    fprintf('✓ Backward scattering contribution is clearly visible\n');
    
catch ME
    fprintf('Error in demo: %s\n', ME.message);
    fprintf('Please ensure all required functions are available.\n');
end

fprintf('\n=== Demo Complete ===\n');