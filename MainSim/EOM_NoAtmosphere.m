function Sdot = EOM_NoAtmosphere(t, S, mu1, mu2, mu3) 
% 3BP vector drag % THREE-BODY EQUATIONS OF MOTION (full Newtonian) 
% S = [r1; v1; r2; v2; r3; v3], with r,v in km and km/s 
% mu = G * m, in km^3/s^2 

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

%% --- Accelerations --- 
a1 = -mu2 * (r12 / r12_mag^3) - mu3 * (r13 / r13_mag^3);
a2 = -mu3 * (r23 / r23_mag^3) - mu1 * (r21 / r21_mag^3);
a3 = -mu1 * (r31 / r31_mag^3) - mu2 * (r32 / r32_mag^3);

%% --- Repack --- 
Sdot = [v1; a1; v2; a2; v3; a3];
