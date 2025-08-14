function [surface_current, current_magnitude] = calculate_surface_current(x_terrain, y_terrain, ...
    x_source, y_source, beta_0, omega, epsilon_0, mu_0, delta_x, n_points)
%% Calculate Surface Current using Forward-Backward Scattering Method
% Implements EFIE-based surface current calculation with forward and backward iterations
% Based on reference implementation in efie.txt
%
% Inputs:
%   x_terrain, y_terrain - Terrain geometry
%   x_source, y_source - Source position
%   beta_0 - Wave number
%   omega - Angular frequency
%   epsilon_0, mu_0 - Material constants
%   delta_x - Spatial step
%   n_points - Number of analysis points
%
% Outputs:
%   surface_current - Complex surface current distribution
%   current_magnitude - Magnitude of surface current

% Initialize surface current array
surface_current = zeros(n_points, 1);
j = 1i; % Imaginary unit

%% Helper Functions (from efie.txt reference)

% Distance from source to point p
R_source_p = @(p) sqrt((x_source - x_terrain(p))^2 + (y_source - y_terrain(p))^2);

% Distance between points p and q
R_p_q = @(p, q) sqrt((x_terrain(q) - x_terrain(p))^2 + (y_terrain(q) - y_terrain(p))^2);

% Incident field radiation (EiRad function from C++ code)
EiRad = @(dist) -((beta_0^2)/(4*omega*epsilon_0)) * hankel_h2(beta_0 * dist);

% Impedance matrix element Z(p,q) 
Z_pq = @(p, q) ((beta_0^2)/(4*omega*epsilon_0)) * hankel_h2(beta_0 * R_p_q(p, q));

% Self impedance Z_self
Z_self = @(i) ((beta_0^2)/(4*omega*epsilon_0)) * ...
    (delta_x - j * ((2*delta_x)/pi) * log((1.781*beta_0*delta_x)/(4*exp(1))));

fprintf('Computing forward scattering iteration...\n');

%% Forward Scattering Calculation (Equation 6 from efie.txt)
% J[0] = EiRad(R_source_p(0),0)/Zself(0)
if n_points > 0
    surface_current(1) = EiRad(R_source_p(1)) / Z_self(1);
end

% Forward iteration: J[p] = (EiRad(R_source_p(p),p) - SUM) / Zself(p)
for p = 2:n_points
    SUM = 0;
    for q = 1:(p-1)
        if q < n_points
            R_segment = R_p_q(q, q+1);
        else
            R_segment = delta_x; % Use grid spacing for last segment
        end
        SUM = SUM + R_segment * Z_pq(p, q) * surface_current(q);
    end
    
    surface_current(p) = (EiRad(R_source_p(p)) - SUM) / Z_self(p);
    
    % Progress indicator
    if mod(p, max(1, floor(n_points/10))) == 0
        fprintf('Forward iteration: %d%% complete\n', round(100*p/n_points));
    end
end

fprintf('Computing backward scattering iteration...\n');

%% Backward Scattering Calculation (Backward iteration from efie.txt)
% J[p] += (-1.0 * SUM) / Zself(p)
for p = (n_points-1):-1:2
    SUM = 0;
    for q = n_points:-1:(p+1)
        if q <= n_points && q > 1
            R_segment = R_p_q(q-1, q);
        else
            R_segment = delta_x; % Use grid spacing for segment
        end
        SUM = SUM + R_segment * Z_pq(p, q) * surface_current(q);
    end
    
    surface_current(p) = surface_current(p) + (-1.0 * SUM) / Z_self(p);
    
    % Progress indicator  
    if mod(n_points-p, max(1, floor(n_points/10))) == 0
        fprintf('Backward iteration: %d%% complete\n', round(100*(n_points-p)/n_points));
    end
end

%% Calculate current magnitude
current_magnitude = abs(surface_current);

fprintf('Surface current range: %.2e to %.2e A/m\n', min(current_magnitude), max(current_magnitude));

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