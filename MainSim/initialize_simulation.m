function S0 = initialize_simulation(rp, ra, TA_start_deg, sc_props,r_ga)
    % INITIALIZE_SIMULATION Computes the initial state vector S0 for the simulation.
    %
    % Inputs:
    %   rp           - Periapsis radius (km)
    %   ra           - Apoapsis radius (km)
    %   TA_start_deg - Starting true anomaly (degrees)
    %   sc_props     - Struct with spacecraft properties:
    %                  m            - Mass (kg)
    %                  Cd           - Drag coefficient
    %                  A            - Area (m^2)
    %                  R_n          - Nose radius (m)
    %                  n_efficiency - Efficiency (not used in ODE, but included)
    %
    % Output:
    %   S0 - Initial state vector [r1; v1; r2; v2; r3; v3; 0]

    % Hardcoded constants (as per main script)
    G = 6.67430e-20; % km^3/kg/s^2
    m1 = 1.024e26; % Neptune [kg]
    m2 = 0; % Triton [kg] (set to 0, gravity negligible)
    mu1 = G * m1; % km^3/s^2
    mu2 = G * m2; % km^3/s^2
    mu3 = G * 160; % km^3/s^2
    rNeptune = 24764; % km
    neptuneAtmoStart = rNeptune + 1000; % km
    r1_0 = [0; 0; 0]; % km
    v1_0 = [0; 0; 0]; % km/s
    data = readtable("neptune-GRAM-avg.txt");
    data.Properties.VariableNames = ["Alt_km", "TotNDm3", "H2NDm3", "HeNDm3", "CH4NDm3", "Denskg_m3", "PresN_m2", "TempK", "Zeta"];
    a_Triton = 354759; % km
    e_Triton = 0;
    h_Triton = sqrt(mu1 * a_Triton * (1 - e_Triton^2));

    % Vehicle properties array for ODE
    vehicle_properites = [sc_props.sc_m, sc_props.sc_Cd, sc_props.sc_A, sc_props.R_n]; % Package for ODE solver
    % n_efficiency = sc_props.n_efficiency; % Defined but not used in provided ODE call

    % Compute initial orbit parameters
    a_probe = (rp + ra) / 2;
    sc_Period = 2 * pi * sqrt(a_probe^3 / mu1);

    % Display input values
    fprintf("--------------------INPUT VALUES--------------------\n")
    fprintf(" SC PROPERTIES:\n")
    fprintf(" * SC Period = %0.2f hours\n", sc_Period / 3600)
    fprintf(" * SC Periapsis Altitude = %0.2f km\n", rp - rNeptune)
    fprintf(" * SC Apoapsis Altitude = %0.2f km\n", ra - rNeptune)
    fprintf(" * SC Cd = %0.2f \n", sc_props.sc_Cd)
    fprintf(" * SC Area = %0.2f m^2\n", sc_props.sc_A)
    fprintf(" * SC mass = %0.2f kg\n", sc_props.sc_m)
    fprintf("----------------------------------------------------\n")

    e_probe = (ra - rp) / (ra + rp);
    h_probe = sqrt(mu1 * a_probe * (1 - e_probe^2));

    % ---- Calculation for Triton's starting TA ----
    % Constants for calc
    r_atm = neptuneAtmoStart; % km
    r_triton = a_Triton; % km

    % -------------------------------------------------------------------------
    % FIND TIME TO GET TO EDGE OF NEPTUNE ATMO
    % -------------------------------------------------------------------------
    TA_start = deg2rad(TA_start_deg); % rad
    % Initial orbit parameters
    p1 = h_probe^2 / mu1;
    % True anomaly at entry (r = r_atm, on descent after apoapsis)
    cos_TA_entry = (p1 / r_atm - 1) / e_probe;
    if abs(cos_TA_entry) > 1
        error('No atmosphere entry: Check rp < r_atm < ra');
    end
    TA_entry = 2 * pi - acos(cos_TA_entry); % > pi rad
    % TOF1: From TA_start to TA_entry
    TOF1 = time_from_TA_to_TA(TA_start, TA_entry, a_probe, e_probe, mu1);

    % -------------------------------------------------------------------------
    % FIND TIME SPENT IN NEPTUNE ATMO
    % -------------------------------------------------------------------------
    [r3_entry, v3_entry] = sv_from_coe([h_probe, e_probe, 0, 0, 0, TA_entry], mu1);
    TA0_T_placeholder = 0; % Arbitrary placeholder (negligible effect since mu2=0)
    coe_Triton_placeholder = [h_Triton, e_Triton, 0, 0, 0, deg2rad(TA0_T_placeholder)];
    [r2_placeholder, v2_placeholder] = sv_from_coe(coe_Triton_placeholder, mu1);
    S0_entry = [r1_0; v1_0; r2_placeholder; v2_placeholder; r3_entry; v3_entry; 0];
    % Aerobraking ODE setup
    dt = 10; % s (step size; adjust as needed)
    currentTime = 0; % Start ODE at entry time (sim starts at t=0)
    totalTime = 10000; % Max time in atm; increase if aerobrake takes longer
    tspan = currentTime:dt:totalTime;
    options = odeset('Events', @event_exit_atmosphere, 'RelTol',1e-8, 'AbsTol',1e-9, 'MaxStep', dt);
    %--- Perform the ODE integration for the spacecraft's trajectory
    [~, ~, te, Se, ~] = ode89(@(t, S) EOM_InAtmosphere(t,S,mu1,mu2,mu3,data,vehicle_properites), tspan, S0_entry, options);
    TOF2 = te; % Time in atmosphere (s)
    fprintf(' - Time spent in atmosphere: %0.2f minutes\n', TOF2 / 60);

    % -------------------------------------------------------------------------
    % FIND TIME FROM EXIT TO INTERCEPT
    % -------------------------------------------------------------------------
    R_exit = Se(13:15);
    V_exit = Se(16:18);
    coe_new = coe_from_sv(R_exit, V_exit, mu1); % [h, e, Omega, i, omega, TA] (TA in rad)
    h_new = coe_new(1);
    e_new = coe_new(2);
    TA_exit = coe_new(6); % rad
    a_new = h_new^2 / mu1 / (1 - e_new^2);
    % True anomaly at intercept (r = r_triton)

    TA_int1 = acos((1/e_new) * ((h_new^2 / ((r_triton+r_ga)*mu1))-1));

    % TOF3: From TA_exit to TA_int1 (as per working code, using TA_int1)
    TOF3 = time_from_TA_to_TA(TA_exit, TA_int1, a_new, e_new, mu1);
    % Total time from start to intercept
    total_time = TOF1 + TOF2 + TOF3;

    % -------------------------------------------------------------------------
    % FIND TRITON START TA
    % -------------------------------------------------------------------------
    Triton_t_to_int = TA_int1 * sqrt((r_triton^3) / mu1) - total_time;
    starting_TA_triton_rad = Triton_t_to_int / sqrt((r_triton^3) / mu1);
    TA0_T = rad2deg(starting_TA_triton_rad);

    % ---- End of calculation ----
    % Debug prints (remove after verification)
    fprintf('Debug: TOF1 = %0.2f s, TOF2 = %0.2f s, TOF3 = %0.2f s\n', TOF1, TOF2, TOF3);
    fprintf('Debug: total_time = %0.2f s (~%0.2f hours)\n', total_time, total_time/3600);
    fprintf('Debug: starting_TA_triton_rad = %0.4f rad (~%0.2f deg)\n', starting_TA_triton_rad, TA0_T);

    % Set Triton initial COE and state
    coe_Triton = [h_Triton, e_Triton, 0, 0, 0, deg2rad(TA0_T)]; % No extra negation needed
    [r2_0, v2_0] = sv_from_coe(coe_Triton, mu1);

    % Set spacecraft starting COE and state
    coe_probe = [h_probe, e_probe, 0, 0, 0, deg2rad(TA_start_deg)];
    [r3_0, v3_0] = sv_from_coe(coe_probe, mu1);

    % Initial state
    S0 = [r1_0; v1_0; r2_0; v2_0; r3_0; v3_0; 0];
end

% ---- Subfunctions ----
function E = ecc_anom(TA, e)
    if e >= 1 || e < 0
        error('Eccentricity must be 0 <= e < 1');
    end
    k = sqrt((1 - e) / (1 + e));
    half_TA = TA / 2;
    half_E = atan(k * tan(half_TA));
    E = 2 * half_E;
    if E < 0
        E = E + 2 * pi;
    end
end

function tof = time_from_TA_to_TA(TA1, TA2, a, e, mu)
    if TA2 < TA1
        TA2 = TA2 + 2 * pi;  % Assume one orbit max
    end
    E1 = ecc_anom(TA1, e);
    M1 = E1 - e * sin(E1);
    E2 = ecc_anom(TA2, e);
    M2 = E2 - e * sin(E2);
    n = sqrt(mu / a^3);
    tof = (M2 - M1) / n;
    if tof < 0
        error('Negative TOF: Check TA order');
    end
end