function [value, isterminal, direction] = event_exit_atmosphere(t, S)
    rNeptune = 24764;  % km
    alt = norm(S(13:15)) - rNeptune;

    % We entered Neptune's atmosphere
    value(1) = alt - 1000;
    isterminal(1) = 1;
    direction(1) = 1;     % only trigger on upward crossing

    % We crashed into the surface of the planet
    value(2) = alt;  % We crashed into the surface of the planet
    isterminal(2) = 1;
    direction(2) = 0;
end