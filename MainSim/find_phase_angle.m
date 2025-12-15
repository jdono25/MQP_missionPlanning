function phase_anlge = find_phase_angle(sc_rp,sc_ra,sc_ta_1)
%find_phase_angle: finds the angle that triton should be at 
%                  given a starting angle of the spacecraft
%   The original math only really works for 

r_Neptune = 24764; % radius of Neptune [km]
mu = 6836529; % 

R_triton = 354759; % Orbital radius of triton, orbit assuemd circular [km]
h_triton = sqrt(R_triton*mu);
n_triton = sqrt(mu/(R_triton^3));

% Properities of SPACECRAFT
a_sc = (sc_rp + sc_ra)/2;
e_sc = (sc_ra - sc_rp)/(sc_rp + sc_ra);
h_sc = sqrt(sc_ra*mu*(1-e_sc));
n_sc = sqrt(mu/(a_sc^3));

ang_int = acos((1/e_sc)*((h_sc^2)/(mu*R_triton)-1));

if sc_ta_1 > 180
    % first find the time it would take to get from periapsis to apoapsis
    t_180 = pi*sqrt((a_sc^3)/mu);

    % Then find the time it takes to get to SC position from apoapsis
    sc_ta_rad = sc_ta_1 - 180;
    u_A = acos((e_sc+cos(sc_ta_rad))/(1+e_sc*cos(sc_ta_rad)));
    M_Asc = u_A - e_sc*sin(sc_ta_rad);
    t_Asc = M_Asc/n_sc;

    tA_total = t_180 + t_Asc;

    % Now we need to find the endpoint when the spacecraft intersects the
    % atmoshpere of Neptune
    ang_int = 2*pi - acos((1/e_sc)*((h_sc^2)/(mu*(r_Neptune+1000))-1));
    u_A2 = acos((e_sc+cos(ang_int))/(1+e_sc*cos(ang_int)));
    M_Asc2 = u_A2 - e_sc*sin(ang_int);
    t_Asc2 = M_Asc2/n_sc;   

    time_to_atmo = t_Asc2 - tA_total;
end


% Find the the time it take the spacecraft to get to the first intercept
% from the periapsis
sc_ta_1 = deg2rad(sc_ta_1);
u_A2 = acos((e_sc+cos(sc_ta_1))/(1+e_sc*cos(sc_ta_1)));
M_Asc2 = u_A2 - e_sc*sin(u_A2);
t_Asc2 = M_Asc2/n_sc;

end