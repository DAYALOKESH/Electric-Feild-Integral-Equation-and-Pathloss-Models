function [surface_current, current_positions] = calculate_scattering(terrain_distance, terrain_height, params)
    % Calculate surface current distribution using Electric Field Integral Equation (EFIE)
    % with Forward and Backward Scattering Method (FBSM)
    %
    % This implementation follows the algorithm from efie.txt C++ reference
    % 
    % Inputs:
    %   terrain_distance - array of terrain distance points (m)
    %   terrain_height - array of terrain height points (m) 
    %   params - structure containing electromagnetic parameters
    %
    % Outputs:
    %   surface_current - complex array of surface current values
    %   current_positions - array of position coordinates for current samples
    
    % Calculate number of line segments based on discretization
    gross_no_steps = length(terrain_distance) - 1;
    total_length = terrain_distance(end) - terrain_distance(1);
    num_segments = round(total_length / params.delta_x);
    
    fprintf('  Discretizing terrain into %d segments\n', num_segments);
    
    % Create uniform spatial grid for current calculation
    current_positions = linspace(terrain_distance(1), terrain_distance(end), num_segments+1);
    current_heights = interp1(terrain_distance, terrain_height, current_positions, 'linear');
    
    % Initialize arrays for surface current calculation
    surface_current = complex(zeros(num_segments, 1));
    
    % Calculate incident field from line source at each point
    incident_field = complex(zeros(num_segments, 1));
    for p = 1:num_segments
        r_source_p = sqrt((params.x_source - current_positions(p))^2 + ...
                         (params.y_source - current_heights(p))^2);
        incident_field(p) = calculate_incident_field(r_source_p, params);
    end
    
    % Forward scattering calculation (transmission direction)
    fprintf('  Computing forward scattering...\n');
    
    % Initialize first current element
    if num_segments > 0
        z_self_0 = calculate_self_impedance(1, current_positions, current_heights, params);
        surface_current(1) = incident_field(1) / z_self_0;
    end
    
    % Forward sweep: calculate current at each point considering previous points
    for p = 2:num_segments
        sum_forward = complex(0.0, 0.0);
        
        % Sum contributions from all previous segments
        for q = 1:(p-1)
            r_pq = sqrt((current_positions(p) - current_positions(q))^2 + ...
                       (current_heights(p) - current_heights(q))^2);
            segment_length = calculate_segment_length(q, current_positions, current_heights);
            z_pq = calculate_mutual_impedance(r_pq, params);
            sum_forward = sum_forward + segment_length * z_pq * surface_current(q);
        end
        
        % Calculate current using EFIE formulation
        z_self_p = calculate_self_impedance(p, current_positions, current_heights, params);
        surface_current(p) = (incident_field(p) - sum_forward) / z_self_p;
    end
    
    % Backward scattering calculation (reflection direction)
    fprintf('  Computing backward scattering...\n');
    
    % Backward sweep: refine current considering later points
    for p = (num_segments-1):-1:1
        sum_backward = complex(0.0, 0.0);
        
        % Sum contributions from all later segments
        for q = (p+1):num_segments
            r_pq = sqrt((current_positions(p) - current_positions(q))^2 + ...
                       (current_heights(p) - current_heights(q))^2);
            segment_length = calculate_segment_length(q, current_positions, current_heights);
            z_pq = calculate_mutual_impedance(r_pq, params);
            sum_backward = sum_backward + segment_length * z_pq * surface_current(q);
        end
        
        % Update current with backward scattering contribution
        z_self_p = calculate_self_impedance(p, current_positions, current_heights, params);
        surface_current(p) = surface_current(p) - sum_backward / z_self_p;
    end
    
    fprintf('  Surface current calculation completed\n');
end

function incident_field = calculate_incident_field(distance, params)
    % Calculate incident electric field from line source
    % Based on EiRad function from efie.txt
    
    argument = params.beta_0 * distance;
    
    % Calculate Hankel function of second kind, order 0
    h0_2 = besselj(0, argument) - 1i * bessely(0, argument);
    
    % Electric field from line source
    incident_field = -((params.beta_0^2) / (4.0 * params.omega * params.epsilon_0)) * h0_2;
end

function z_mutual = calculate_mutual_impedance(distance, params)
    % Calculate mutual impedance between two segments
    % Based on Z function from efie.txt
    
    argument = params.beta_0 * distance;
    
    % Calculate Hankel function of second kind, order 0
    h0_2 = besselj(0, argument) - 1i * bessely(0, argument);
    
    % Mutual impedance
    z_mutual = ((params.beta_0^2) / (4.0 * params.omega * params.epsilon_0)) * h0_2;
end

function z_self = calculate_self_impedance(index, positions, heights, params)
    % Calculate self impedance of a segment
    % Based on Zself function from efie.txt
    
    % Calculate segment length
    segment_length = calculate_segment_length(index, positions, heights);
    
    if segment_length < params.tolerance
        segment_length = params.delta_x; % Use default discretization
    end
    
    % Self impedance calculation (with singularity treatment)
    real_part = ((params.beta_0^2) / (4.0 * params.omega * params.epsilon_0)) * segment_length;
    
    % Imaginary part includes logarithmic singularity treatment
    log_term = log((1.781 * params.beta_0 * segment_length) / (4.0 * exp(1)));
    imag_part = -((params.beta_0^2) / (4.0 * params.omega * params.epsilon_0)) * ...
                ((2.0 * segment_length) / pi) * log_term;
    
    z_self = complex(real_part, imag_part);
end

function length = calculate_segment_length(index, positions, heights)
    % Calculate length of a terrain segment
    
    if index >= length(positions)
        length = 0;
        return;
    end
    
    dx = positions(index+1) - positions(index);
    dy = heights(index+1) - heights(index);
    length = sqrt(dx^2 + dy^2);
end