function satEmpheris(timeVals,posVals)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


fileName = 's_ephemeris.e';
fid = fopen(fileName, 'w');
fprintf(fid, 'stk.v.12.0\n');
fprintf(fid, 'BEGIN Ephemeris\n');
fprintf(fid, sprintf('NumberOfEphemerisPoints %d\n',length(posVals)));
fprintf(fid, 'CoordinateSystem Inertial\n');
fprintf(fid, 'CentralBody Neptune\n');
fprintf(fid, 'InterpolationMethod Lagrage\n');
fprintf(fid, 'InterpolationOrder 5\n');
fprintf(fid, 'EphemerisTimePos\n');

for i = 1:length(posVals)
    fprintf(fid, '%0.4f %0.4f %0.4f %0.4f\n', timeVals(i), posVals(i,1), posVals(i,2), 0);
end

fprintf(fid, 'END Ephemeris\n');