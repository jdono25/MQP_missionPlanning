function [v_new, dv_flyby_mag, speed_gain, turn_angle_deg, v_final_mag, dv_correction] = perform_triton_flyby_and_correction(r_sc, v_sc, r_triton, v_triton, altitude_km, rp_desired_alt_km, mu_neptune, plot_flag)
    % PERFORM_TRITON_FLYBY_AND_CORRECTION Computes flyby velocity update and correction burn delta-V.
    %
    % Inputs:
    %   r_sc              - Spacecraft position [x, y, z] (km, but converted to m internally)
    %   v_sc              - Spacecraft velocity [vx, vy, vz] (km/s, to m/s)
    %   r_triton          - Triton position [x, y, z] (km to m)
    %   v_triton          - Triton velocity [vx, vy, vz] (km/s to m/s)
    %   altitude_km       - Flyby altitude above Triton (km)
    %   rp_desired_alt_km - Desired periapsis altitude above Neptune for correction (km)
    %   mu_neptune        - Neptune mu (km^3/s^2)
    %   plot_flag         - True to plot orbits (optional, default false)
    %
    % Outputs:
    %   v_new            - Updated velocity after flyby [vx, vy, vz] (km/s)
    %   dv_flyby_mag     - Flyby delta-V magnitude (m/s)
    %   speed_gain       - Speed gain from flyby (m/s)
    %   turn_angle_deg   - Turn angle (deg)
    %   v_final_mag      - Final speed after flyby (m/s)
    %   dv_correction    - Correction burn delta-V at apoapsis (m/s)

    if nargin < 8
        plot_flag = false;
    end
    
    % Constants (in m)
    mu_neptune_m = mu_neptune * 1e9;  % km^3/s^2 to m^3/s^2
    mu_triton = 1.4276e12;  % m^3/s^2
    R_triton_body = 1353e3;  % m
    rNeptune = 24764e3;  % m
    
    % Convert inputs to m, m/s
    r_sc_m = r_sc * 1e3;
    v_sc_m = v_sc * 1e3;
    r_triton_m = r_triton * 1e3;
    v_triton_m = v_triton * 1e3;
    
    % Rotate frame to align r_sc with x-axis
    r_mag = norm(r_sc_m);
    theta = atan2(r_sc_m(2), r_sc_m(1));
    R = [cos(theta), sin(theta); -sin(theta), cos(theta)];  % Rotation matrix (2D)
    v_sc_rot = R * v_sc_m(1:2);
    v_triton_rot = R * v_triton_m(1:2);
    
    % Flyby calculation (adapted from triton_gravity_assist)
    v_inf_in_vec = v_sc_rot - v_triton_rot;
    v_inf_mag = norm(v_inf_in_vec);
    if v_inf_mag < 1e-5
        warning('No gravity assist possible.');
        v_new = v_sc;
        dv_flyby_mag = 0;
        speed_gain = 0;
        turn_angle_deg = 0;
        v_final_mag = norm(v_sc_m);
        dv_correction = 0;
        return;
    end
    r_p = R_triton_body + altitude_km * 1000;  % m
    denominator = 1 + (r_p * v_inf_mag^2) / mu_triton;
    delta = 2 * asin(1 / denominator);
    turn_angle_deg = rad2deg(delta);
    
    % Rotation matrices
    R_pos = [cos(delta), -sin(delta); sin(delta), cos(delta)];
    R_neg = [cos(delta), sin(delta); -sin(delta), cos(delta)];
    v_inf_out_1 = R_pos * v_inf_in_vec;
    v_inf_out_2 = R_neg * v_inf_in_vec;
    
    v_final_vec_1 = v_triton_rot + v_inf_out_1;
    v_final_vec_2 = v_triton_rot + v_inf_out_2;
    mag_1 = norm(v_final_vec_1);
    mag_2 = norm(v_final_vec_2);
    
    if mag_1 > mag_2
        v_final_rot = v_final_vec_1;
        v_final_mag = mag_1;
    else
        v_final_rot = v_final_vec_2;
        v_final_mag = mag_2;
    end
    
    dv_flyby_mag = norm(v_final_rot - v_sc_rot);
    speed_gain = v_final_mag - norm(v_sc_m);
    
    % Rotate back to original frame
    R_inv = R';  % Transpose for inverse
    v_new_m = R_inv * v_final_rot;
    v_new = [v_new_m; 0] / 1000;  % Back to km/s, add z=0
    
    % Correction burn delta-V (adapted from get_correcton_brun, without flyby redo)
    r_mag_m = norm(r_sc_m);
    v_mag_m = v_final_mag;
    h_vec = cross(r_sc_m, [v_new_m; 0]);
    h_mag = norm(h_vec);
    e_vec = cross([v_new_m; 0], h_vec) / mu_neptune_m - r_sc_m / r_mag_m;
    e_mag = norm(e_vec);
    a_m = -mu_neptune_m / (2 * (v_mag_m^2 / 2 - mu_neptune_m / r_mag_m));
    rp_m = a_m * (1 - e_mag);
    ra_m = a_m * (1 + e_mag);
    v_apo_m = (mu_neptune_m / h_mag) * (1 + e_mag * cos(pi));  % Should match v_mag_m at apo, but computed for verification
    
    rp_desired_m = (rp_desired_alt_km + 24764) * 1000;
    e_desired = (ra_m - rp_desired_m) / (ra_m + rp_desired_m);
    h_desired = sqrt(mu_neptune_m * ra_m * (1 - e_desired));
    v_apo_desired = (mu_neptune_m / h_desired) * (1 + e_desired * cos(pi));
    dv_correction = abs(v_apo_m - v_apo_desired);  % m/s
    
    % Optional plotting (adapted from plot_triton_flyby)
    if plot_flag
        % Generate orbit points pre-flyby (using original v_sc_m)
        [x_in, y_in, ~, T_in] = get_orbit_points(r_sc_m, v_sc_m, mu_neptune_m);
        % Post-flyby
        [x_out, y_out, ~, T_out] = get_orbit_points(r_sc_m, [v_new_m; 0], mu_neptune_m);
        % Triton orbit
        theta_plot = linspace(0, 2*pi, 200);
        x_triton_orb = a_Triton * cos(theta_plot) * 1e3;  % m (assuming a_Triton in km)
        y_triton_orb = a_Triton * sin(theta_plot) * 1e3;
        
        figure;
        hold on; axis equal; grid on;
        % Draw Neptune (exaggerated)
        R_nep_draw = 24622e3;
        rectangle('Position', [-R_nep_draw, -R_nep_draw, 2*R_nep_draw, 2*R_nep_draw], 'Curvature', [1,1], 'FaceColor', [0.3, 0.3, 0.8], 'EdgeColor', 'none');
        % Triton orbit
        plot(x_triton_orb, y_triton_orb, 'k--', 'LineWidth', 1);
        % Triton position
        plot(r_triton_m(1), r_triton_m(2), 'ko', 'MarkerFaceColor', 'white', 'MarkerSize', 8);
        % Pre-flyby orbit
        plot(x_in, y_in, 'b-', 'LineWidth', 1.5);
        % Post-flyby orbit
        plot(x_out, y_out, 'r-', 'LineWidth', 1.5);
        xlabel('Distance (m)'); ylabel('Distance (m)');
        title(sprintf('Triton Flyby (Alt: %d km)', altitude_km));
        hold off;
    end
end

% Helper subfunction (from your code)
function [x, y, a, Period] = get_orbit_points(r_vec, v_vec, mu)
    r = norm(r_vec);
    v = norm(v_vec);
    E = (v^2)/2 - mu/r;
    a = -mu / (2*E);
    h_vec = cross(r_vec, v_vec);
    e_vec = cross(v_vec, h_vec)/mu - r_vec/r;
    e = norm(e_vec);
    omega = atan2(e_vec(2), e_vec(1));
    Period = 2*pi*sqrt(abs(a)^3/mu);
    nu = linspace(0, 2*pi, 360);
    dist = a * (1 - e^2) ./ (1 + e*cos(nu));
    x_per = dist .* cos(nu);
    y_per = dist .* sin(nu);
    x = x_per * cos(omega) - y_per * sin(omega);
    y = x_per * sin(omega) + y_per * cos(omega);
end