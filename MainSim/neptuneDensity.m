function density = neptuneDensity(altitude, data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Interpolate the density based on the provided altitude
density = interp1(data.Alt_km, data.Denskg_m3, altitude, 'linear', 'extrap');

% Ensure altitude is within the range of the data for accurate interpolation
if altitude < min(data.Alt_km) || altitude > max(data.Alt_km)
    warning('Altitude is out of range. Extrapolation may lead to inaccurate results.');
end