function [value, isterminal, direction] = event_enter_atmosphere(t, S)
    rNeptune = 24764;  % km
    alt = norm(S(13:15)) - rNeptune;
    r2 = S(7:9);
    r3 = S(13:15);
    r23 = r2 - r3;
    r23_mag = norm(r23);

    % We entered Neptune's atmosphere
    value(1) = alt - 1000;   % boundary at 1000 km
    isterminal(1) = 1;    % stop integration
    direction(1) = -1;    % only trigger on downward crossing

    % We crashed into the surface of the planet
    value(2) = alt;  % We crashed into the surface of the planet
    isterminal(2) = 1;
    direction(2) = 0;
end