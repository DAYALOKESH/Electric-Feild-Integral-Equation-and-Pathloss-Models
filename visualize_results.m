function visualize_results(x_terrain, y_terrain, current_magnitude, field_magnitude, lambda)
%% Visualize Results with Professional Quality Plots
% Creates publication-quality visualizations of surface current and electric field
%
% Inputs:
%   x_terrain - Distance points (m)
%   y_terrain - Terrain height profile (m) 
%   current_magnitude - Surface current magnitude (A/m)
%   field_magnitude - Electric field magnitude (V/m and dB)
%   lambda - Electromagnetic wavelength (m)

%% Figure Setup
% Set default figure properties for publication quality
% Note: groot not available in Octave, using alternative approach
try
    set(groot, 'DefaultFigureColor', 'white');
    set(groot, 'DefaultAxesBox', 'on');
    set(groot, 'DefaultAxesLineWidth', 1.2);
    set(groot, 'DefaultAxesFontSize', 12);
    set(groot, 'DefaultAxesFontWeight', 'bold');
    set(groot, 'DefaultLineLineWidth', 2);
catch
    % Octave fallback - set properties individually for each plot
end

%% Figure 1: Terrain Profile and Surface Current Distribution
figure('Position', [100, 100, 1200, 800], 'Name', 'FBSM Analysis Results');

% Subplot 1: Terrain Profile
subplot(2, 2, 1);
plot(x_terrain, y_terrain, 'k-', 'LineWidth', 2);
grid on;
xlabel('Distance (m)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Height (m)', 'FontSize', 12, 'FontWeight', 'bold');
title('Terrain Profile (X.04)', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0, max(x_terrain)]);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');

% Add wavelength reference
hold on;
y_range = max(y_terrain) - min(y_terrain);
text(0.05*max(x_terrain), min(y_terrain) + 0.1*y_range, ...
    sprintf('λ = %.2f m', lambda), 'FontSize', 10, 'BackgroundColor', 'white');
hold off;

% Subplot 2: Surface Current Magnitude
subplot(2, 2, 2);
semilogy(x_terrain, current_magnitude, 'r-', 'LineWidth', 2);
grid on;
xlabel('Distance (m)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Surface Current Magnitude (A/m)', 'FontSize', 12, 'FontWeight', 'bold');
title('Surface Current Distribution (FBSM)', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0, max(x_terrain)]);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');

% Add statistics text box
current_max = max(current_magnitude);
current_min = min(current_magnitude);
current_mean = mean(current_magnitude);
text(0.6*max(x_terrain), 0.7*current_max, ...
    sprintf('Max: %.2e A/m\nMean: %.2e A/m\nMin: %.2e A/m', ...
    current_max, current_mean, current_min), ...
    'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', 'black');

% Subplot 3: Electric Field Magnitude (Linear Scale)
subplot(2, 2, 3);
semilogy(x_terrain, field_magnitude.linear, 'b-', 'LineWidth', 2);
grid on;
xlabel('Distance (m)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Electric Field Magnitude (V/m)', 'FontSize', 12, 'FontWeight', 'bold');
title('Electric Field Distribution (Linear)', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0, max(x_terrain)]);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');

% Subplot 4: Electric Field Magnitude (dB Scale)  
subplot(2, 2, 4);
% Remove infinite values for plotting
field_dB_plot = field_magnitude.dB;
field_dB_plot(isinf(field_dB_plot)) = min(field_dB_plot(~isinf(field_dB_plot))) - 20;

plot(x_terrain, field_dB_plot, 'b-', 'LineWidth', 2);
grid on;
xlabel('Distance (m)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Electric Field (dB relative to 1 V/m)', 'FontSize', 12, 'FontWeight', 'bold');
title('Electric Field Distribution (Normalized dB)', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0, max(x_terrain)]);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');

% Add dynamic range information
field_range = max(field_dB_plot) - min(field_dB_plot);
text(0.6*max(x_terrain), 0.8*max(field_dB_plot), ...
    sprintf('Dynamic Range: %.1f dB', field_range), ...
    'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', 'black');

%% Enhance overall figure appearance
try
    sgtitle('Forward and Backward Scattering Method - Electromagnetic Analysis', ...
        'FontSize', 16, 'FontWeight', 'bold');
catch
    % Octave doesn't have sgtitle, skip main title
end

% Adjust subplot spacing
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

%% Figure 2: Combined Analysis View
figure('Position', [200, 200, 1000, 600], 'Name', 'FBSM Combined Analysis');

% Combined surface current and electric field plot
try
    % Try MATLAB's yyaxis (newer versions)
    yyaxis left;
    semilogy(x_terrain, current_magnitude, 'r-', 'LineWidth', 2.5);
    ylabel('Surface Current Magnitude (A/m)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
    set(gca, 'YColor', 'r');

    yyaxis right;
    semilogy(x_terrain, field_magnitude.linear, 'b-', 'LineWidth', 2.5);
    ylabel('Electric Field Magnitude (V/m)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');
    set(gca, 'YColor', 'b');
catch
    % Fallback for Octave - use plotyy
    [ax, h1, h2] = plotyy(x_terrain, log10(current_magnitude), x_terrain, log10(field_magnitude.linear));
    set(h1, 'Color', 'r', 'LineWidth', 2.5);
    set(h2, 'Color', 'b', 'LineWidth', 2.5);
    set(ax(1), 'YColor', 'r');
    set(ax(2), 'YColor', 'b');
    ylabel(ax(1), 'Surface Current Magnitude (log scale)', 'FontSize', 12, 'Color', 'r');
    ylabel(ax(2), 'Electric Field Magnitude (log scale)', 'FontSize', 12, 'Color', 'b');
end

xlabel('Distance (m)', 'FontSize', 12, 'FontWeight', 'bold');
title('Surface Current and Electric Field vs Distance', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
xlim([0, max(x_terrain)]);
set(gca, 'XMinorTick', 'on');

% Add legend
legend({'Surface Current', 'Electric Field'}, 'Location', 'best', 'FontSize', 11);

% Add electromagnetic parameters text box
param_text = sprintf('Frequency: 970 MHz\nWavelength: %.2f m\nAnalysis Points: %d\nMax Distance: %.0f m', ...
    lambda, length(x_terrain), max(x_terrain));
text(0.02*max(x_terrain), 0.5*max(field_magnitude.linear), param_text, ...
    'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', 'black', 'VerticalAlignment', 'top');

%% Save figures (optional)
% Uncomment to save figures
% print(gcf, 'FBSM_Combined_Analysis', '-dpng', '-r300');

fprintf('Visualization completed successfully!\n');
fprintf('Generated professional-quality plots showing:\n');
fprintf('  - Terrain profile from X.04 data\n');
fprintf('  - Surface current distribution from FBSM\n'); 
fprintf('  - Electric field magnitude (linear and dB scales)\n');
fprintf('  - Combined analysis view\n');

end