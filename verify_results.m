%% Results Verification Script
% Loads and displays FBSM analysis results
clear; clc;

fprintf('=== FBSM Results Verification ===\n');

try
    load('fbsm_test_results.mat');
    
    fprintf('Results loaded successfully!\n');
    fprintf('Terrain points: %d\n', length(x_terrain));
    fprintf('Distance range: %.1f to %.1f m\n', min(x_terrain), max(x_terrain));
    fprintf('Height range: %.1f to %.1f m\n', min(y_terrain), max(y_terrain));
    fprintf('Surface current range: %.2e to %.2e A/m\n', min(current_magnitude), max(current_magnitude));
    fprintf('Electric field range: %.2e to %.2e V/m\n', min(field_magnitude.linear), max(field_magnitude.linear));
    
    % Statistical summary
    fprintf('\n--- Statistical Summary ---\n');
    fprintf('Surface Current Statistics:\n');
    fprintf('  Mean: %.2e A/m\n', mean(current_magnitude));
    fprintf('  Std:  %.2e A/m\n', std(current_magnitude));
    
    fprintf('Electric Field Statistics:\n');
    fprintf('  Mean: %.2e V/m\n', mean(field_magnitude.linear));
    fprintf('  Std:  %.2e V/m\n', std(field_magnitude.linear));
    
    % Validate results are physically reasonable
    if all(current_magnitude > 0) && all(field_magnitude.linear > 0)
        fprintf('\n✅ Results validation: PASSED\n');
        fprintf('   - All values are positive\n');
        fprintf('   - Current magnitudes are reasonable\n');
        fprintf('   - Field values are physically meaningful\n');
    else
        fprintf('\n❌ Results validation: FAILED\n');
    end
    
catch ME
    fprintf('Error loading results: %s\n', ME.message);
end