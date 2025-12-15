clear; close all; clc

% pool = parpool("Processes");
animate_orbits = true;

%% Create the Tisserand Plot
% Nep_Triton_Tisserand(1,5,0.25)

%% Sim code
% Constants
G = 6.67430e-20; % km^3/kg/s^2
m1 = 1.024e26; % Neptune [kg]
m2 = 2.14e22; % Triton [kg]
m2 = 0; % Triton [kg] % Note: Set to 0 as per your code, but mu2 uses original m2? Adjust if needed
m3 = 160; % Spacecraft [kg]
mu1 = G * m1; %km^3/s^2
mu2 = G * m2; %km^3/s^2
mu3 = G * m3; %km^3/s^2
%% Neptune Data
rNeptune = 24764; % km
neptuneAtmoStart = rNeptune + 1000; % km
% Starting state vector
r1_0 = [0; 0; 0]; % km
v1_0 = [0; 0; 0]; % km/s
% Neptune GRAM data
data = readtable("neptune-GRAM-avg.txt");
data.Properties.VariableNames = ["Alt_km", "TotNDm3", "H2NDm3", "HeNDm3", "CH4NDm3", "Denskg_m3", "PresN_m2", "TempK", "Zeta"];
%% Define Spacecraft properties
Triton_period = 507764.16;
% Physical Properties
% powerTrade = zeros(n_peri, n_periods);
sc_m = 160; % [kg]
sc_Cd = 1.024;
% sc_Cd = 0;
sc_A = 1.5; % [m^2]
R_n = 0.22; % [m]
n_efficiency = 0.02; % Effieceny of the power transfer through heat pumps
vehicle_properites = [sc_m, sc_Cd, sc_A, R_n, n_efficiency]; % Package them all up for OsDE solver
sc_props.sc_m = sc_m;
sc_props.sc_Cd = sc_Cd;
sc_props.sc_A = sc_A;
sc_props.R_n = R_n;
sc_props.n = n_efficiency;

rp = 25217.5;
ra = 730846;
a_probe = (rp+ra)/2;
sc_Period = 2*pi*sqrt(a_probe^3/mu1);
e_probe = (ra - rp)/(ra + rp);
h_probe = sqrt(mu1 * a_probe * (1 - e_probe^2));

%% Set up triton properties
Triton_period = 507764.16; % seconds
a_Triton = 354759;
e_Triton = 0;
h_Triton = sqrt(mu1 * a_Triton * (1 - e_Triton^2));

% ---- Inserted calculation for Triton's starting TA ----
% Constants for calc
r_atm = neptuneAtmoStart; % km
r_triton = a_Triton; % km
omega_triton = -2 * pi / Triton_period; % rad/s, negative for retrograde

% -------------------------------------------------------------------------
%                FIND TIME TO GET TO EDGE OF NEPTUNE ATMO
% -------------------------------------------------------------------------
TA_start_deg = 180; % <----------------------------------------------------- Starting TA for spacecraft [deg]
r_ga = 1000;        % <----------------------------------------------------- gravity assist distance to Triton [km]
TA_start = deg2rad(TA_start_deg); % rad

S0 = initialize_simulation(rp,ra,TA_start_deg,sc_props,r_ga);

%% ODE simulations
dt = 10;
number_sc_orbits = 3;
totalTime = number_sc_orbits*sc_Period;
currentTime = 0;

[S_complete, t_complete] = run_simulation(S0,dt,totalTime, number_sc_orbits,1000, data, vehicle_properites, rp-rNeptune, r_ga);
%% Plotting orbits
incl = 0; % degrees
yRot = [cos(deg2rad(incl)), 0, sin(deg2rad(incl));
                0         , 1,          0        ;
       -sin(deg2rad(incl)), 0, cos(deg2rad(incl))];

figure(1);
hold on; grid on; axis equal;
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title('Neptune–Triton–Spacecraft System');
view(3);
view([0 90])
xline(-ra,'k--')

% Neptune plot
plot3(S_complete(:,1),S_complete(:,2),S_complete(:,3),'k')
plot3(S_complete(1,1),S_complete(1,2),S_complete(1,3),'k*')
hold on

% Neptune sphere
R_neptune = 24764; % km
[u,v] = meshgrid(linspace(0,2*pi,50), linspace(0,pi,50));
x_sphere = R_neptune * cos(u) .* sin(v);
y_sphere = R_neptune * sin(u) .* sin(v);
z_sphere = R_neptune * cos(v);

surf(x_sphere, y_sphere, z_sphere, ...
     'FaceColor','b','EdgeColor','none','FaceAlpha',0.8);

% Triton plot
S_complete(:,7:9) = S_complete(:,7:9) * yRot;
plot3(S_complete(:,7), S_complete(:,8), S_complete(:,9), 'red--');
plot3(S_complete(1,7), S_complete(1,8), S_complete(1,9), 'red*');

% % Spacecraft orbit
% rotate the spacecrafts orbit: 
S_complete(:,13:15) = S_complete(:,13:15) * yRot;
plot3(S_complete(:,13),S_complete(:,14),S_complete(:,15),'blue')
plot3(S_complete(1,13),S_complete(1,14),S_complete(1,15),'blue*')
plot3(S_complete(end,13),S_complete(end,14),S_complete(end,15),'red*')

%% Plotting Sutton Graves correlation

% figure(2)
% plot(t_complete, S_complete(:,19))

%% Animations
% Initialize moving markers
figure(1)
sc_marker = plot3(S_complete(1,13),S_complete(1,14), S_complete(1,15), 'bo', 'MarkerFaceColor','b', 'MarkerSize',6);
triton_marker  = plot3(S_complete(1,7), S_complete(1,8), S_complete(1,9), 'ro', 'MarkerFaceColor','r', 'MarkerSize',6);
% neptune_marker = plot3(x_N(1), y_N(1), z_N(1), 'ko', 'MarkerFaceColor','c', 'MarkerSize',10);

if animate_orbits
    for i = 1:length(t_complete)
        % Update positions
        set(sc_marker, 'XData', S_complete(i,13), ...
                       'YData', S_complete(i,14), ...
                       'ZData', S_complete(i,15));

        set(triton_marker, 'XData', S_complete(i,7), ...
                           'YData', S_complete(i,8), ...
                           'ZData', S_complete(i,9));
        % 
        % set(neptune_marker, 'XData', x_N(i), ...
        %                     'YData', y_N(i), ...
        %                     'ZData', z_N(i));

        drawnow limitrate nocallbacks
        % pause(0.0001);
    end
end




%% ---- Helper functions (add these to your script or as subfunctions) ----
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