function figureHandle = Nep_Triton_Tisserand(min_Vinf,max_Vinf,step)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%%
%{
Tisserand graph for MQP
%}

%Gravitational parameter of Triton, m^3/s^2
mu_second_vec = 1427.6;

%Gravitational parameter of Neptune, m^3/s^2 
mu_prim = 6836529; % 
r_Neptune = 24764;

R_sec_vec = 354759; %Orbital radius of planets
V_sec_vec = sqrt(mu_prim./R_sec_vec); % km/s, orbiting velocity
R_planet = 1352.5; %km, Equatorial radius of planets

%Define minimum flyby radius
R_ga_vec = R_planet + 100;      % default min flyby radius

%Define maximum magnitude of V-infinity
% V_inf_vec = 1.5*V_sec_vec;
% max_Vinf = 8;
% sprintf("Maximum v_inf is currently %0.3f", max_Vinf)
% V_inf_vec(4:end) = 50;
hax=axes;

% Initial Orbit
% rp = r_Neptune + 450;
Triton_period = 507772.8;
sc_Period = Triton_period/2;
rp_initial = r_Neptune + 550;
a_probe = (mu_prim*((sc_Period/(2*pi))^2))^(1/3); % km
ra_initial = 2*a_probe - rp_initial;
e_probe = (ra_initial - rp_initial)/(ra_initial + rp_initial);
h_probe = sqrt(mu_prim * a_probe * (1 - e_probe^2));

energy_probe_initial = -0.5*((mu_prim^2)/(h_probe^2)*(1-e_probe^2));

c = colormap(turbo(13));

%Set up of V-infinity step
vInf_vec = max_Vinf;
% step = 0.1;
% min_Vinf = 5;
numLines = length(min_Vinf:step:vInf_vec);
numind = 1:numLines;

ra_rp_Data = cell(1, numLines);

for j = 1:length(R_sec_vec)%1:length(mu_second_vec)  %iterate over each planet
    vInfInd = 1;
    %Store planet's characteristics
    mu_second = mu_second_vec(j);
    R_sec = R_sec_vec(j);
    V_sec = V_sec_vec(j);
    R_ga = R_ga_vec(j);
    
    for V_inf_sec = min_Vinf:step:max_Vinf %Iterater over V-infinity values
        
        %Define maximum turn angle, based on min flyby radius
        delta_turn = 2*asin(1/(1+R_ga*V_inf_sec^2/mu_second));
        %Define number of points between each flyby marker
        delta_n = 100;
        

        %alph is angle between planet velocity and V-infinity
        if V_inf_sec>V_sec
            %If V-infinity is 
            alph = [0:delta_turn/delta_n:(pi-acos(V_sec/V_inf_sec))...
                        0.999*(pi-acos(V_sec/V_inf_sec))];
        else
            alph = [0:delta_turn/delta_n:pi pi];
        end
        
        ind = 1:delta_n:length(alph);
        
        %Compute Spacecraft velocity and flight-path angle wrt Sun
        V_sc = sqrt(V_sec^2+V_inf_sec^2-2*V_sec*V_inf_sec*cos(pi-alph));
        fpa = asin(V_inf_sec./V_sc.*sin(pi-alph));
        
        %Compute Spacecraft's orbit characteristics wrt Sun
        energy = V_sc.^2/2 - mu_prim/R_sec;
        h = R_sec*V_sc.*cos(fpa);
        eccen = sqrt(1+2*energy.*h.^2/mu_prim^2);
        a = -mu_prim/2./energy;
        rp = a.*(1-eccen);
        % rp_filter1 = rp < r_Neptune + 10000;
        % rp_filter2 = rp > r_Neptune;
        % rp_fitled = rp_filter1 == rp_filter2;
        ra = a.*(1+eccen);

        rp_filter1 = ra < 10e7;
        rp_filter2 = rp > r_Neptune+400 ;
        rp_fitled = rp_filter1 == rp_filter2;
        
        ra_rp_Data{j,numind(vInfInd)} = [ra',rp'];
        vInfInd = vInfInd + 1;
        ra(ra<0)=nan;
        
        % filtered_ras = ra(rp_fitled);
        % filtered_rps = rp(rp_fitled);

        figureHandle = figure(10);
        % loglog(ra,rp,'-k','Color',c(j,:),'DisplayName',num2str(V_inf_sec))
        loglog(ra(rp_fitled),rp(rp_fitled),'-k','Color',c(j,:),'HandleVisibility','off')
        loglog(ra(ind),rp(ind),'*','Color',c(j,:),'HandleVisibility','off')
        hold on
    end
end

% Plot the resonance
ra_rp_forRes = ra_rp_Data{1,3};
% The format is planet orbits/SC orbits
resonances =       [1/1  , 1/2  , 1/3 , 2/3  , 11/10  , 3/2];
resonances_names = strings(1,length(resonances));
for indx = 1:length(resonances)
    [num, denom] = rat(resonances(indx));
    resonances_names(indx) = sprintf('%g:%g', num, denom);
end

colors = cool(length(resonances));
for indx = 1:length(resonances)
    resonance = resonances(indx);
    ra_res = ra_rp_forRes(:,1);
    rp_res = 2*R_sec_vec*((resonance^2)^(1/3)) - ra_res;
    figure(10)
    loglog(ra_res,rp_res,'DisplayName',resonances_names(indx),'Color',colors(indx,:))
end

warning('off','MATLAB:xyzline:NegativeDataIgnored')
yline(24764+1000,'b','DisplayName','Start Neptune Atmo')
warning('on','MATLAB:xyzline:NegativeDataIgnored')
legend show

%%Figure 1 Label
figure(10)
% plot([0 max(rp)],[195 195],'r-','LineWidth',2)
hold on
% plot(6.8413e5,25314,'r*');
grid on
title("r_a vs r_p plot")
xlabel('r_a [km]')
ylabel('r_p [km]')
xlim([2.9541e+05, 1.4884e+06])
ylim([2.0331e+04, 6.7868e+05])

