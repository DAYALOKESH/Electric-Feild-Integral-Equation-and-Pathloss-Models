function [electric_field, field_positions] = calculate_electric_field(surface_current, current_positions, terrain_distance, terrain_height, params)
    % Calculate electric field distribution above the terrain surface
    % Based on the electric field calculation from efie.txt
    %
    % Inputs:
    %   surface_current - complex array of surface current values
    %   current_positions - position coordinates for current samples
    %   terrain_distance - original terrain distance data
    %   terrain_height - original terrain height data
    %   params - electromagnetic parameters structure
    %
    % Outputs:
    %   electric_field - complex array of total electric field values
    %   field_positions - position coordinates for field observation points
    
    num_current_points = length(surface_current);
    field_positions = current_positions;
    electric_field = complex(zeros(size(field_positions)));
    
    fprintf('  Computing electric field at %d observation points...\n', length(field_positions));
    
    % Calculate electric field at each observation point
    for obs_idx = 1:length(field_positions)
        % Initialize scattered field component
        scattered_field = complex(0.0, 0.0);
        
        % Sum contributions from all current segments up to observation point
        for current_idx = 1:min(obs_idx, num_current_points)
            if current_idx <= length(surface_current)
                % Calculate distance from current segment to observation point
                % Observation point is at height + offset above terrain
                obs_height = interp1(terrain_distance, terrain_height, field_positions(obs_idx), 'linear') + params.obs_height_offset;
                current_height = interp1(terrain_distance, terrain_height, current_positions(current_idx), 'linear');
                
                r_surf_obs = sqrt((field_positions(obs_idx) - current_positions(current_idx))^2 + ...
                                 (obs_height - current_height)^2);
                
                % Calculate segment length for integration
                segment_length = calculate_field_segment_length(current_idx, current_positions, terrain_distance, terrain_height);
                
                % Calculate impedance for this distance
                z_distance = calculate_field_impedance(r_surf_obs, params);
                
                % Add contribution to scattered field
                scattered_field = scattered_field + surface_current(current_idx) * segment_length * z_distance;
            end
        end
        
        % Calculate incident field at observation point
        obs_height = interp1(terrain_distance, terrain_height, field_positions(obs_idx), 'linear') + params.obs_height_offset;
        r_source_obs = sqrt((params.x_source - field_positions(obs_idx))^2 + ...
                           (params.y_source - obs_height)^2);
        incident_field_obs = calculate_incident_field_obs(r_source_obs, params);
        
        % Total field = Incident field - Scattered field
        % (Note: Sign convention from efie.txt)
        electric_field(obs_idx) = incident_field_obs - scattered_field;
    end
    
    fprintf('  Electric field calculation completed\n');
end

function incident_field = calculate_incident_field_obs(distance, params)
    % Calculate incident electric field at observation point
    % Similar to incident field calculation but for observation points
    
    argument = params.beta_0 * distance;
    
    % Calculate Hankel function of second kind, order 0
    h0_2 = besselj(0, argument) - 1i * bessely(0, argument);
    
    % Electric field from line source
    incident_field = -((params.beta_0^2) / (4.0 * params.omega * params.epsilon_0)) * h0_2;
end

function z_field = calculate_field_impedance(distance, params)
    % Calculate impedance for electric field computation
    % Based on Z function from efie.txt for field calculations
    
    argument = params.beta_0 * distance;
    
    % Calculate Hankel function of second kind, order 0
    h0_2 = besselj(0, argument) - 1i * bessely(0, argument);
    
    % Impedance for field calculation
    z_field = ((params.beta_0^2) / (4.0 * params.omega * params.epsilon_0)) * h0_2;
end

function segment_length = calculate_field_segment_length(index, positions, terrain_distance, terrain_height)
    % Calculate segment length for field integration
    % Uses terrain-following discretization
    
    if index >= length(positions)
        segment_length = positions(2) - positions(1); % Use uniform spacing
        return;
    end
    
    % Get heights at segment endpoints
    height1 = interp1(terrain_distance, terrain_height, positions(index), 'linear');
    height2 = interp1(terrain_distance, terrain_height, positions(index+1), 'linear');
    
    % Calculate 3D segment length along terrain surface
    dx = positions(index+1) - positions(index);
    dy = height2 - height1;
    segment_length = sqrt(dx^2 + dy^2);
    
    % Ensure minimum segment length for numerical stability
    if segment_length < 1e-12
        segment_length = dx; % Use horizontal distance if heights are equal
    end
end