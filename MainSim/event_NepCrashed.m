function [value, isterminal, direction] = event_NepCrashed(t, S)
    rNeptune = 24764;  % km
    alt = norm(S(13:15)) - rNeptune;

    % We crashed into the surface of the planet
    value = alt;  % We crashed into the surface of the planet
    isterminal = 1;
    direction = 0;
end
