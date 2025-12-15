function [value, isterminal, direction] = event_enter_atmosphere_new(t, S, gravity_assist_alt)
    rNeptune = 24764;  % km
    r_triton = 1353;  % km
    alt = norm(S(13:15)) - rNeptune;
    r2 = S(7:9);
    r3 = S(13:15);
    r23 = r2 - r3;
    r23_mag = norm(r23);

    % Enter Neptune's atmosphere
    value(1) = alt - 1000;
    isterminal(1) = 1;
    direction(1) = -1;

    % Crash into Neptune
    value(2) = alt;
    isterminal(2) = 1;
    direction(2) = 0;
    
    % Close to Triton for flyby
    value(3) = r23_mag - (r_triton + gravity_assist_alt);
    isterminal(3) = 1;
    direction(3) = -1;

    % Apoapsis (vr crosses 0 from + to -) time for correction burn
    % r_sc = S(13:15);
    % v_sc = S(16:18);
    % vr = dot(r_sc, v_sc) / norm(r_sc);  % km/s
    % value(4) = vr;
    % isterminal(4) = 1;
    % direction(4) = 1;
end