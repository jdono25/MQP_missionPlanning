function Sdot = EOM_InAtmosphere(t, S, mu1, mu2, mu3, data, vehicle_properties) 
% 3BP vector drag % THREE-BODY EQUATIONS OF MOTION (full Newtonian) 
% S = [r1; v1; r2; v2; r3; v3], with r,v in km and km/s 
% mu = G * m, in km^3/s^2 

k = 6.28e-5; % [sqrt(kg)/m] - value taken from chatGPT when asked to convert, 
% converted from k = 0.042 kg/(s·m^{3/2}·atm^{1/2}) from pg 65 in Sutton-Graves paper 

%% Spacecraft Properties 
m_sc = vehicle_properties(1); % kg
Cd = vehicle_properties(2);
A = vehicle_properties(3);
R_n = vehicle_properties(4); % vehicle nose radius [m]

%% --- Unpack states --- 
r1 = S(1:3); v1 = S(4:6); 
r2 = S(7:9); v2 = S(10:12); 
r3 = S(13:15); v3 = S(16:18); 
% --- Relative position vectors --- 
r12 = r1 - r2; 
r13 = r1 - r3; 
r21 = r2 - r1; 
r23 = r2 - r3; 
r31 = r3 - r1; 
r32 = r3 - r2; 
% --- Magnitudes --- 
r12_mag = norm(r12); 
r21_mag = norm(r21); 
r13_mag = norm(r13); 
r31_mag = norm(r31); 
r23_mag = norm(r23); 
r32_mag = norm(r32);

R_NEPTUNE_KM = 24764;
sc_alt_km = r13_mag - R_NEPTUNE_KM;
% Neptune_atmo_start = 1000; % km 
v_km_s = [v3(1); v3(2); v3(3)];
v_m_s  = v_km_s * 1000;           % m/s
v_mag  = norm(v_m_s);
% if sc_alt_km < Neptune_atmo_start && (sc_alt_km >= 0) && drag
% Atmospheric density lookup: **assumes neptuneDensity expects altitude in km**
rho = neptuneDensity(sc_alt_km, data);   % returns kg/m^3 per your table

% if v_mag > 0
    % drag accel (m/s^2), vector opposite v
    a_drag_m_s2 = -0.5 * rho * Cd * A * v_mag * v_m_s / m_sc;

    % convert to km/s^2
    a_drag_km = a_drag_m_s2 / 1000;
% end
% else
%     a_drag_km = zeros(3, 1); % Initialize drag acceleration to zero if no drag is applied
%     rho = 0;
% end


%% --- Accelerations --- 
a1 = -mu2 * (r12 / r12_mag^3) - mu3 * (r13 / r13_mag^3);
a2 = -mu3 * (r23 / r23_mag^3) - mu1 * (r21 / r21_mag^3);
a3 = -mu1 * (r31 / r31_mag^3) - mu2 * (r32 / r32_mag^3) + a_drag_km;

sg = k*sqrt(rho/R_n)*(v_mag)^3;

%% --- Repack --- 
Sdot = [v1; a1; v2; a2; v3; a3; sg];
% Sdot = [v1; a1; v2; a2; v3; a3];
