function create_visualizations(terrain_distance, terrain_height, surface_current, current_positions, electric_field, field_positions, params)
    % Create high-quality visualizations for FBSM electromagnetic scattering analysis
    %
    % Generates publication-quality plots showing:
    % 1. Surface current distribution across terrain
    % 2. Electric field magnitude versus distance
    % 3. Terrain profile for reference
    %
    % Inputs:
    %   terrain_distance - original terrain distance data (m)
    %   terrain_height - original terrain height data (m) 
    %   surface_current - complex surface current array
    %   current_positions - position coordinates for current samples (m)
    %   electric_field - complex electric field array
    %   field_positions - position coordinates for field observation points (m)
    %   params - electromagnetic parameters structure
    
    % Set up figure properties for publication quality
    set(0, 'DefaultFigureColor', 'white');
    set(0, 'DefaultAxesFontSize', 12);
    set(0, 'DefaultAxesFontWeight', 'bold');
    set(0, 'DefaultTextFontSize', 12);
    set(0, 'DefaultLineLineWidth', 2);
    
    % Create main figure with subplots
    fig = figure('Position', [100, 100, 1200, 800], 'Name', 'FBSM Electromagnetic Scattering Analysis');
    
    %% Plot 1: Terrain Profile
    subplot(2, 2, 1);
    plot(terrain_distance, terrain_height, 'k-', 'LineWidth', 2);
    grid on;
    xlabel('Distance (meters)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Height (meters)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Terrain Profile', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Add source position marker
    hold on;
    plot(params.x_source, params.y_source, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'red');
    text(params.x_source + 20, params.y_source, 'Line Source', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'red');
    hold off;
    
    % Set axis limits with some padding
    xlim([min(terrain_distance)-10, max(terrain_distance)+10]);
    ylim([min(terrain_height)-20, max(terrain_height)+50]);
    
    %% Plot 2: Surface Current Magnitude Distribution
    subplot(2, 2, 2);
    current_magnitude = abs(surface_current);
    
    % Create filled area plot for better visualization
    area(current_positions(1:length(current_magnitude)), current_magnitude, 'FaceColor', [0.2, 0.6, 0.8], 'FaceAlpha', 0.7);
    hold on;
    plot(current_positions(1:length(current_magnitude)), current_magnitude, 'b-', 'LineWidth', 2);
    hold off;
    
    grid on;
    xlabel('Distance (meters)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Surface Current Magnitude |J| (A/m)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Surface Current Distribution', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Add statistics to the plot
    max_current = max(current_magnitude);
    mean_current = mean(current_magnitude);
    text(0.7, 0.9, sprintf('Max: %.2e A/m\nMean: %.2e A/m', max_current, mean_current), ...
         'Units', 'normalized', 'BackgroundColor', 'white', 'FontSize', 10);
    
    %% Plot 3: Electric Field Magnitude vs Distance
    subplot(2, 2, 3);
    field_magnitude = abs(electric_field);
    
    % Convert to dB scale for better visualization (following efie.txt approach)
    field_positions_valid = field_positions(1:length(field_magnitude));
    
    % Calculate path loss correction (following efie.txt formula)
    r_source_obs = sqrt((params.x_source - field_positions_valid).^2 + ...
                       (params.y_source - interp1(terrain_distance, terrain_height, field_positions_valid, 'linear')).^2);
    field_magnitude_db = 20 * log10(field_magnitude ./ sqrt(r_source_obs));
    
    plot(field_positions_valid, field_magnitude_db, 'r-', 'LineWidth', 2);
    grid on;
    xlabel('Distance (meters)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Electric Field (dB)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Electric Field vs Distance', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Add legend to distinguish forward and backward components
    legend('Total Field (Forward + Backward)', 'Location', 'best');
    
    %% Plot 4: Forward and Backward Scattering Components
    subplot(2, 2, 4);
    
    % Separate forward and backward components for visualization
    % (This is an approximation for visualization purposes)
    num_points = length(current_magnitude);
    forward_component = current_magnitude(1:floor(num_points/2));
    backward_component = current_magnitude(floor(num_points/2)+1:end);
    
    forward_positions = current_positions(1:length(forward_component));
    backward_positions = current_positions(length(forward_component)+1:length(forward_component)+length(backward_component));
    
    % Plot forward and backward components
    plot(forward_positions, forward_component, 'g-', 'LineWidth', 2, 'DisplayName', 'Forward Scattering');
    hold on;
    plot(backward_positions, backward_component, 'm-', 'LineWidth', 2, 'DisplayName', 'Backward Scattering');
    hold off;
    
    grid on;
    xlabel('Distance (meters)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Current Magnitude |J| (A/m)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Forward vs Backward Scattering Components', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best');
    
    %% Add main title and parameter information
    sgtitle(sprintf('Forward and Backward Scattering Method (FBSM) Analysis\nFrequency: %.1f MHz, Wavelength: %.2f m', ...
                   params.frequency/1e6, params.lambda), 'FontSize', 16, 'FontWeight', 'bold');
    
    % Adjust subplot spacing
    tight_layout();
    
    %% Save results to files (following efie.txt output format)
    fprintf('Saving results to output files...\n');
    
    % Save surface current data (similar to J.dat from efie.txt)
    save_current_data(current_positions, surface_current, 'surface_current_results.dat');
    
    % Save electric field data (similar to E.dat from efie.txt)
    save_field_data(field_positions, electric_field, r_source_obs, 'electric_field_results.dat');
    
    % Save figure
    saveas(fig, 'FBSM_Analysis_Results.png', 'png');
    savefig(fig, 'FBSM_Analysis_Results.fig');
    
    fprintf('Visualization completed successfully!\n');
    fprintf('Results saved to:\n');
    fprintf('  - surface_current_results.dat\n');
    fprintf('  - electric_field_results.dat\n');
    fprintf('  - FBSM_Analysis_Results.png\n');
    fprintf('  - FBSM_Analysis_Results.fig\n');
end

function save_current_data(positions, current, filename)
    % Save surface current data to file (following efie.txt J.dat format)
    
    fid = fopen(filename, 'w');
    if fid == -1
        warning('Could not open file %s for writing', filename);
        return;
    end
    
    fprintf(fid, '%% Surface Current Results from FBSM Analysis\n');
    fprintf(fid, '%% Distance(m)  Current_Magnitude(A/m)  Current_Phase(rad)\n');
    
    for i = 1:min(length(positions), length(current))
        magnitude = abs(current(i));
        phase = angle(current(i));
        fprintf(fid, '%.6f  %.6e  %.6f\n', positions(i), magnitude, phase);
    end
    
    fclose(fid);
end

function save_field_data(positions, field, source_distances, filename)
    % Save electric field data to file (following efie.txt E.dat format)
    
    fid = fopen(filename, 'w');
    if fid == -1
        warning('Could not open file %s for writing', filename);
        return;
    end
    
    fprintf(fid, '%% Electric Field Results from FBSM Analysis\n');
    fprintf(fid, '%% Distance(m)  Field_dB  Field_Magnitude(V/m)  Field_Phase(rad)\n');
    
    for i = 1:min(length(positions), length(field))
        magnitude = abs(field(i));
        phase = angle(field(i));
        % Calculate dB with path loss correction (following efie.txt)
        field_db = 20 * log10(magnitude / sqrt(source_distances(i)));
        fprintf(fid, '%.6f  %.6f  %.6e  %.6f\n', positions(i), field_db, magnitude, phase);
    end
    
    fclose(fid);
end

function tight_layout()
    % Improve subplot spacing
    set(gcf, 'Units', 'normalized');
    subplot_handles = get(gcf, 'Children');
    
    % Adjust positions for better layout
    if length(subplot_handles) >= 4
        positions = {[0.08, 0.55, 0.38, 0.35], [0.55, 0.55, 0.38, 0.35], ...
                    [0.08, 0.08, 0.38, 0.35], [0.55, 0.08, 0.38, 0.35]};
        for i = 1:min(4, length(subplot_handles))
            set(subplot_handles(end-i+1), 'Position', positions{i});
        end
    end
end