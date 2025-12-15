function [S_complete, t_complete] = run_simulation(S0, dt, totalTime, number_sc_orbits,neptuneAtmoStart, data, vehicle_properites, rp_desired_alt_km, gravity_assist_alt)
    % RUN_SIMULATION Performs the ODE simulation for the spacecraft trajectory.
    %
    % Inputs:
    %   S0                   - Initial state vector (19x1)
    %   dt                   - Time step (s)
    %   number_sc_orbits     - Number of spacecraft orbits to simulate
    %   sc_Period            - Spacecraft orbital period (s)
    %   neptuneAtmoStart     - Neptune atmosphere start radius (km)
    %   mu1, mu2, mu3        - Gravitational parameters (km^3/s^2)
    %   data                 - Neptune GRAM data table
    %   vehicle_properites   - Vehicle properties array [m, Cd, A, R_n]
    %   n_efficiency         - Efficiency for power calculation
    %   sc_A                 - Spacecraft area (m^2)
    %
    % Outputs:
    %   S_complete           - Complete state trajectory matrix
    %   t_complete           - Complete time vector
    %
    % Assumptions:
    %   Function handles @event_exit_atmosphere, @event_enter_atmosphere,
    %   @EOM_InAtmosphere, @EOM_NoAtmosphere, @coe_from_sv, @sv_from_coe
    %   are defined in the calling scope or on the MATLAB path.
        
    G = 6.67430e-20; % km^3/kg/s^2
    m1 = 1.024e26; % Neptune [kg]
    % m2 = 2.14e22; % Triton [kg]
    m2 = 0; % Triton [kg] % Note: Set to 0 as per your code, but mu2 uses original m2? Adjust if needed
    m3 = 160; % Spacecraft [kg]
    mu1 = G * m1; %km^3/s^2
    mu2 = G * m2; %km^3/s^2
    mu3 = G * m3; %km^3/s^2

    % ----------- Create Tisserand Plot
    tissFig = Nep_Triton_Tisserand(1,5,0.25);

    % Add the initial orbit to the tisserand plot. First have to find the
    % periapsis and apoapsis heights
    coe = coe_from_sv(S0(13:15), S0(16:18), mu1);
    coe(6) = deg2rad(180);
    [R0_apo, ~] = sv_from_coe(coe, mu1);
    apo0 = norm(R0_apo);

    coe = coe_from_sv(S0(13:15), S0(16:18), mu1);
    coe(6) = deg2rad(0);
    [R0_peri, ~] = sv_from_coe(coe, mu1);
    peri0 = norm(R0_peri);

    figure(tissFig)
    loglog(apo0,peri0,'ro',DisplayName='Starting Orbit',MarkerSize=8)
    text(apo0,peri0,'#0')

    % Triton Properties
    R_triton_body = 1353.4;

    currentTime = 0;
    S_complete = [];
    t_complete = [];
    
    % Check if already in Neptune's atmosphere
    if norm(S0(13:15)) < neptuneAtmoStart
        inAtmo = true;
        atmopass_num = 1; % Initialize atmosphere pass counter
    else
        inAtmo = false;
        S0(end) = [];
        atmopass_num = 0;
    end
    
    % Run ODE simulation
    while currentTime < totalTime && atmopass_num <= number_sc_orbits && (totalTime - currentTime) > dt
        if inAtmo
            fprintf("Atmosphere pass #%d\n", atmopass_num)
            tspan = currentTime:dt:totalTime;
            options = odeset('Events', @event_exit_atmosphere, 'RelTol',1e-8, 'AbsTol',1e-9, 'MaxStep', dt);
            % Perform the ODE integration for the spacecraft's trajectory
            [t, S, te, Se, ie] = ode89(@(t, S) EOM_InAtmosphere(t, S, mu1, mu2, mu3, data, vehicle_properites), tspan, S0, options);
            S_complete = [S_complete; S];
            t_complete = [t_complete; t];
            currentTime = t(end);
            inAtmo = false;
            atmopass_num = atmopass_num + 1;

            % Display the results
            fprintf(" - Time spent in atmosphere: %0.2f minutes\n", (te - currentTime)/60)
           
            % Old state vector used to find the apoapsis before the aerobraking
            coe = coe_from_sv(S0(13:15), S0(16:18), mu1);
            coe(6) = deg2rad(180);
            [R0_apo, V0_apo] = sv_from_coe(coe, mu1);
            apo0 = norm(R0_apo);
            apo_v0 = norm(V0_apo);
            coe = coe_from_sv(S0(13:15), S0(16:18), mu1);
            coe(6) = deg2rad(0);
            [R0_peri, V0_peri] = sv_from_coe(coe, mu1);
            peri0 = norm(R0_peri);
            peri_v0 = norm(V0_peri);

            % New exit state vector used to find the new apoapsis altitude
            coe_e = coe_from_sv(Se(13:15), Se(16:18), mu1);
            coe_e(6) = deg2rad(180);
            [Re_apo, Ve_apo] = sv_from_coe(coe_e, mu1);
            apo_alt_e = norm(Re_apo);
            apo_v_e = norm(Ve_apo);
            coe_e = coe_from_sv(Se(13:15), Se(16:18), mu1);
            coe_e(6) = deg2rad(0);
            [Re_peri, Ve_peri] = sv_from_coe(coe_e, mu1);
            peri_alt_e = norm(Re_peri);
            peri_v_e = norm(Ve_peri);

            % ---------------- ADD rp ra TO TISSERAND PLOT
            figure(tissFig)
            loglog(apo_alt_e,peri_alt_e,'r*',DisplayName=sprintf('Aerobrake#%d',atmopass_num))
            text(apo_alt_e+200,peri_alt_e+200,sprintf('#%d',atmopass_num))

            deltaV = (norm(S0(16:18)) - norm(Se(16:18))) * 1000;
            fprintf(" - Change in apoapsis altitude: %0.2f km\n", apo0 - apo_alt_e)
            % fprintf(" - Change in apoapsis veclotiy: %0.3f km/s\n", apo_v0 - apo_v_e)
            % fprintf(" - Change in periapsis altitude: %0.2f km\n", peri_alt0 - peri_alt_e)
            % fprintf(" - Change in periapsis veclotiy: %0.3f km/s\n", peri_v0 - peri_v_e)
            fprintf(" - Change in velocity is : %0.2f m/s\n", deltaV)
            % Power Generated
            % Calculate power generated during the atmosphere pass
            % powerGenerated = Se(19) * sc_A * n_efficiency;
            % dqdt = diff(S(:,19)) ./ diff(t);
            dqdt = gradient(S(:,19), t); % More accurate than the diff function
            % powerGenerated = max(dqdt) * sc_A * n_efficiency;
            powerGenerated = (max(S(:,19)) * vehicle_properites(3) * vehicle_properites(end)) / 3600;
            fprintf(" - Power generated during atmosphere pass: %0.2f W\n", powerGenerated);
            
        else
            tspan = currentTime:dt:totalTime;
            options = odeset('Events', @(t, S) event_enter_atmosphere_new(t, S, gravity_assist_alt), 'RelTol',1e-8, 'AbsTol',1e-9, 'MaxStep', dt);
            [t, S, te, Se, ie] = ode89(@(t, S) EOM_NoAtmosphere(t, S, mu1, mu2, mu3), tspan, S0, options);
            S = [S, zeros(size(S,1),1)];
            S_complete = [S_complete; S];
            t_complete = [t_complete; t];
            currentTime = t(end);
            if ie == 1
                inAtmo = true;
            end
        end
        if ie == 1
            if inAtmo
                S0 = [Se'; 0];
            else
                S0 = Se';
                S0(end) = [];
            end
        elseif ie == 2
            sprintf("Crashed into the surface of Neptune at t = %d", te)
            break
        elseif ie == 3
            % Triton flyby event
            fprintf("Triton flyby at t = %0.2f s\n", te)
            Se = Se';  % Column vector
            r_sc = Se(13:15);
            v_sc = Se(16:18);
            r_triton = Se(7:9);
            v_triton = Se(10:12);
            alt_km = norm(r_sc - r_triton) - R_triton_body;  % Actual alt at stop
            [v_new, dv_flyby_mag, speed_gain, turn_angle_deg, v_final_mag, ~] = perform_triton_flyby_and_correction(r_sc, v_sc, r_triton, v_triton, alt_km, rp_desired_alt_km, mu1, false);
            fprintf(" - Flyby delta-V: %0.2f m/s, speed gain: %0.2f m/s, turn angle: %0.2f deg\n", dv_flyby_mag, speed_gain, turn_angle_deg)
            Se(16:18) = v_new;
    
            % % Get the peri and apo from after the flyby
            % % Then add to Tisserand Plot
            % coe = coe_from_sv(Se(13:15), Se(16:18), mu1);
            % coe(6) = deg2rad(180);
            % [R0_apo, ~] = sv_from_coe(coe, mu1);
            % apof = norm(R0_apo);
            % 
            % coe = coe_from_sv(Se(13:15), Se(16:18), mu1);
            % coe(6) = deg2rad(0);
            % [R0_peri, ~] = sv_from_coe(coe, mu1);
            % perif = norm(R0_peri);
            % 
            % a_f = (apof + perif)/2;
            % e_f = 1 - perif / a_f;
            % p_f = a_f * (1 - e_f^2);
            % 
            % % Get the peri and apo from after the flyby
            % % Then add to Tisserand Plot
            % coe = coe_from_sv(r_sc, v_sc, mu1);
            % coe(6) = deg2rad(180);
            % [R0_apo, ~] = sv_from_coe(coe, mu1);
            % apoi = norm(R0_apo);
            % 
            % coe = coe_from_sv(r_sc, v_sc, mu1);
            % coe(6) = deg2rad(0);
            % [R0_peri, ~] = sv_from_coe(coe, mu1);
            % peri = norm(R0_peri);
            % 
            % a_i = (apoi + peri)/2;
            % e_i = 1 - peri / a_i;
            % p_i = a_f * (1 - e_i^2);
            
        
            % figure(tissFig)
            % loglog(apof,perif,'g*',DisplayName=sprintf('Triton Flyb#%d',atmopass_num))
            % text(apof+200,perif+200,sprintf('#%d',atmopass_num))

            S0 = Se;  % Update state, remain out of atmo
        elseif ie == 4
            % Apoapsis correction event
            fprintf("Apoapsis correction at t = %0.2f s\n", te)
            Se = Se';  % Column vector
            r_sc = Se(13:15);
            v_sc = Se(16:18);
            r_triton = Se(7:9);  % Unused for correction, but pass dummy if needed
            v_triton = Se(10:12);  % Unused
            [~, ~, ~, ~, ~, dv_correction] = perform_triton_flyby_and_correction(r_sc, v_sc, r_triton, v_triton, 0, rp_desired_alt_km, mu1, false);  % Alt=0 dummy
            fprintf(" - Correction burn delta-V: %0.2f m/s\n", dv_correction)
            v_dir = v_sc / norm(v_sc);
            v_new = v_sc + (dv_correction / 1000) * v_dir;  % dv in m/s, v in km/s
            Se(16:18) = v_new;
            
            % Get the peri and apo from after the flyby
            % Then add to Tisserand Plot
            coe = coe_from_sv(Se(13:15), Se(16:18), mu1);
            coe(6) = deg2rad(180);
            [R0_apo, ~] = sv_from_coe(coe, mu1);
            apo0 = norm(R0_apo);
        
            coe = coe_from_sv(Se(13:15), Se(16:18), mu1);
            coe(6) = deg2rad(0);
            [R0_peri, ~] = sv_from_coe(coe, mu1);
            peri0 = norm(R0_peri);
        
            figure(tissFig)
            loglog(apo0,peri0,'b*',DisplayName=sprintf('ΔV Correction#%d',atmopass_num))
            text(apo0+200,peri0+200,sprintf('#%d',atmopass_num))

            S0 = Se;  % Update state, remain out of atmo
        elseif isempty(ie)
            disp("Simulation completed successfully")
        end
        if isempty(ie)
            disp("Simulation completed successfully")
        end
    end
end