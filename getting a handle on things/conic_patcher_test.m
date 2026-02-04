format compact
clear
clc
clf reset

%----------

error_tolerance = 1e-11;
unit_vec_size = 6e6;

G = 6.6743e-11;

lunar_mass = 7.35e22; %kg
earth_mass = 5.972e24;

lunar_radius = 1.7374e6;
earth_radius = 6.371e6; %m

mu = G*earth_mass;

lunar_orbital_radius = 	363300e3; %m
lunar_eccentricity = 0.0549;
longitude_ascending_lunar = 0;
argument_of_peri_lunar = 0;
inc_lunar = 5.145; 

r_initial = [
cosd(longitude_ascending_lunar), -sind(longitude_ascending_lunar), 0
sind(longitude_ascending_lunar), cosd(longitude_ascending_lunar), 0
0, 0, 1
]*[
1, 0, 0
0, cosd(inc_lunar), -sind(inc_lunar)
0, sind(inc_lunar), cosd(inc_lunar)
]*[
cosd(argument_of_peri_lunar), -sind(argument_of_peri_lunar), 0
sind(argument_of_peri_lunar), cosd(argument_of_peri_lunar), 0
0, 0, 1
]*[0; lunar_orbital_radius; 0];

semi_major_axis = norm(r_initial)/(1-lunar_eccentricity);
inital_velocity = sqrt((mu/semi_major_axis)*((1+lunar_eccentricity)/(1-lunar_eccentricity)));

r_dot_initial = [
cosd(longitude_ascending_lunar), -sind(longitude_ascending_lunar), 0
sind(longitude_ascending_lunar), cosd(longitude_ascending_lunar), 0
0, 0, 1
]*[
1, 0, 0
0, cosd(inc_lunar), -sind(inc_lunar)
0, sind(inc_lunar), cosd(inc_lunar)
]*[
cosd(argument_of_peri_lunar), -sind(argument_of_peri_lunar), 0
sind(argument_of_peri_lunar), cosd(argument_of_peri_lunar), 0
0, 0, 1
]*[-inital_velocity; 0; 0];

state_initial = [r_initial(1:3).',r_dot_initial(1:3).'];
orbital_period = sqrt(((semi_major_axis^3)/mu)*(2*pi)^2);

orbit_timespan = linspace(0,orbital_period*1.1,3e3);

[lunar_timerange, lunar_state_matrix] = ode45(@(timerange, state_matrix)g_dynamics_twobody(timerange, state_matrix,earth_mass),orbit_timespan,state_initial,odeset('Reltol',error_tolerance));
lunar_orbit_statematrix = [lunar_timerange,lunar_state_matrix];

lunar_SOI_radius = semi_major_axis*(1-lunar_eccentricity)*(lunar_mass / (3*(earth_mass+lunar_mass)) ) ^ (1/3);


%% ejecting out of lunar SOI

SOI_ejection_lat = 0; %extra degrees, parameter but not shown
SOI_ejection_long = 0; %extra degrees, parameter

r_SOI_enter = [
1, 0, 0
0, cosd(SOI_ejection_lat), -sind(SOI_ejection_lat)
0, sind(SOI_ejection_lat), cosd(SOI_ejection_lat)
]*[
cosd(SOI_ejection_long), -sind(SOI_ejection_long), 0
sind(SOI_ejection_long), cosd(SOI_ejection_long), 0
0, 0, 1    
]*[0;-lunar_SOI_radius;0];

r_SOIluna_enter = r_SOI_enter.' + lunar_orbit_statematrix(1,2:4);

%% target MEO
MEO_radius = 2.7e8 + earth_radius %target altitude, parameter

MEO_target_inc = 5.145; %specified in sim
MEO_target_long_ascend = 0; %degrees parameter

MEO_target_TA = 100; %degrees, parameter

MEO_target_argument_peri = 0;

MEO_target_eccentricity = 0; 

target_semimajor = MEO_radius/(1-MEO_target_eccentricity);
target_orbital_period= sqrt(((target_semimajor^3)/mu)*(2*pi)^2);

MEO_target_velocity = sqrt((mu/target_semimajor)*((1+MEO_target_eccentricity)/(1-MEO_target_eccentricity)));

target_start_r = [
cosd(MEO_target_long_ascend), -sind(MEO_target_long_ascend), 0
sind(MEO_target_long_ascend), cosd(MEO_target_long_ascend), 0
0, 0, 1
]*[
1, 0, 0
0, cosd(MEO_target_inc), -sind(MEO_target_inc)
0, sind(MEO_target_inc), cosd(MEO_target_inc)
]*[
cosd(MEO_target_argument_peri), -sind(MEO_target_argument_peri), 0
sind(MEO_target_argument_peri), cosd(MEO_target_argument_peri), 0
0, 0, 1
]*[0; -MEO_radius; 0];
target_start_r = target_start_r.';

target_start_v = [
cosd(MEO_target_long_ascend), -sind(MEO_target_long_ascend), 0
sind(MEO_target_long_ascend), cosd(MEO_target_long_ascend), 0
0, 0, 1
]*[
1, 0, 0
0, cosd(MEO_target_inc), -sind(MEO_target_inc)
0, sind(MEO_target_inc), cosd(MEO_target_inc)
]*[
cosd(MEO_target_argument_peri), -sind(MEO_target_argument_peri), 0
sind(MEO_target_argument_peri), cosd(MEO_target_argument_peri), 0
0, 0, 1
]*[MEO_target_velocity; 0; 0];
target_start_v = target_start_v.';

target_n = [ %target orbital plane normal vector
cosd(MEO_target_long_ascend), -sind(MEO_target_long_ascend), 0
sind(MEO_target_long_ascend), cosd(MEO_target_long_ascend), 0
0, 0, 1
]*[
1, 0, 0
0, cosd(MEO_target_inc), -sind(MEO_target_inc)
0, sind(MEO_target_inc), cosd(MEO_target_inc)
]*[
cosd(MEO_target_argument_peri), -sind(MEO_target_argument_peri), 0
sind(MEO_target_argument_peri), cosd(MEO_target_argument_peri), 0
0, 0, 1
]*[0; 0; 1];
target_n = target_n.';

target_state_initial = [target_start_r,target_start_v];
target_orbit_timespan = linspace(0,target_orbital_period,2e3);

[target_timerange, target_state_matrix] = ode45(@(timerange, state_matrix)g_dynamics_twobody(timerange, state_matrix,earth_mass),target_orbit_timespan,target_state_initial,odeset('Reltol',error_tolerance));
target_orbit_statematrix = [target_timerange,target_state_matrix];

TA_list = TA_from_statematrix(target_orbit_statematrix,target_n);
target_orbit_statematrix = [target_orbit_statematrix, TA_list.'];

target_start_interp = interp_statematrix_TA(target_orbit_statematrix, MEO_target_TA);
target_r = target_start_interp(2:4);
target_v = target_start_interp(5:7);

transfer_dt = 0.5; %parameter, normalised

[v_1, v_2, delta_t, dv_match] = construct_transfer(r_SOIluna_enter, target_r, transfer_dt, mu, target_v); %dt will be a param

state_initial_earthtransit = [r_SOIluna_enter,v_1];
orbit_timespan_earthtransit = linspace(0,delta_t*1.1,3e3); %overshoot dt on purpose

[timerange_earthtransit, state_matrix_earthtransit] = ode45(@(timerange_earthtransit, state_matrix_earthtransit) g_dynamics_twobody(timerange_earthtransit, state_matrix_earthtransit, earth_mass), orbit_timespan_earthtransit, state_initial_earthtransit, odeset('Reltol',error_tolerance));
earthtransit_statematrix = [timerange_earthtransit,state_matrix_earthtransit];

%prune earthtransit statematrix using bisection
inds_sample = [1,height(earthtransit_statematrix)];
for n=1:50
    s_1 = norm(target_r - earthtransit_statematrix(inds_sample(1), 2:4));
    s_2 = norm(target_r - earthtransit_statematrix(inds_sample(2), 2:4));
    ind_mean = floor(mean(inds_sample));
    if s_1 >= s_2
        inds_sample(1) = ind_mean;
    else
        inds_sample(2) = ind_mean;
    end
    if abs(inds_sample(1) - inds_sample(2)) < 2
        break
    end
end
ind_cut = max(inds_sample);
earthtransit_statematrix(ind_cut:end,:)=[];


%find thrust req
impulse_assumption_angle = 10;
r_final = earthtransit_statematrix(end,2:4);
impulse_assumption_time = nan;
for n=height(earthtransit_statematrix):-4:1
    r_spec = earthtransit_statematrix(n,2:4);
    theta_r = acosd(dot(r_final,r_spec) / (norm(r_final)*norm(r_spec)));
    if theta_r > impulse_assumption_angle
        impulse_assumption_time = earthtransit_statematrix(end,1) - earthtransit_statematrix(n,1);
        break
    end
end

%making sure the earth orbit doesn't cross the lunar SOI
no_SOI_reentry = true;
for n=1:height(earthtransit_statematrix)
    orbital_r = norm(earthtransit_statematrix(n,2:4));
    if orbital_r > (lunar_orbital_radius - lunar_SOI_radius) || n==1
        lunar_check_statematrix(n,:) = interp_statematrix_timebetween(lunar_orbit_statematrix, lunar_orbit_statematrix(:,1), earthtransit_statematrix(n,1), earthtransit_statematrix(n,1));
        if norm(lunar_check_statematrix(2:4) - earthtransit_statematrix(n,2:4)) < lunar_SOI_radius
            error("re-enters lunar SOI")
            no_SOI_reentry = false;
            break
        end
    end

end

%% running back in lunar SOI to surface

SOI_enter_stateinital = earthtransit_statematrix(1,2:7) - lunar_check_statematrix(1,2:7);
SOI_backtrack_max_dt = sqrt(( ( (lunar_SOI_radius/2)^3 ) / (lunar_mass*G) ) * (2*pi)^2); %max dt allowed for backsolve back down to lunar surface
SOI_backtrack_timespan = linspace(SOI_backtrack_max_dt,0,2e3);

[timerange_lunarbacktrack, state_matrix_lunarbacktrack] = ode45(@(timerange_lunarbacktrack, state_matrix_lunarbacktrack) g_dynamics_twobody(timerange_lunarbacktrack, state_matrix_lunarbacktrack, lunar_mass), SOI_backtrack_timespan, SOI_enter_stateinital, odeset('Reltol',error_tolerance));
lunarbacktrack_statematrix = [timerange_lunarbacktrack,state_matrix_lunarbacktrack];

%do we hit the moon?
ind_eject = 0;
reaches_surface = false;
orbital_r_prev = 0;
for n=1:height(lunarbacktrack_statematrix)
    orbital_r_lunar = norm(lunarbacktrack_statematrix(n,2:4));
    lunar_altitude_series(n) = orbital_r_lunar;
    if orbital_r_lunar < lunar_radius
        ind_eject = n;
        reaches_surface = true;
        break
    end
    orbital_r_prev = orbital_r_lunar;
end
if reaches_surface
    closest_altitude = 0;
else
    closest_altitude = min(lunar_altitude_series) %if we miss the surface, how close do we get anyway?
end

stays_in_SOI = true;
ejection_possible = false;
avoids_collision = true;
if reaches_surface
    lunarbacktrack_statematrix(ind_eject+1:end,:) = [];
    surface_interp = interp1([orbital_r_lunar, orbital_r_prev],[0,1],lunar_radius);
    surface_eject_time = interp1([0,1],[lunarbacktrack_statematrix(end,1),lunarbacktrack_statematrix(end-1,1)],surface_interp);
    lunar_eject_statematrix = interp_statematrix_timebetween(lunarbacktrack_statematrix, lunarbacktrack_statematrix(:,1), surface_eject_time, surface_eject_time);
    lunarbacktrack_statematrix(end,1:7) = lunar_eject_statematrix;
    
    %finding exact revolution time
    [~,ind_rev] = mink(abs(lunar_orbit_statematrix(:,2)),5);
    ind_rev = sort(ind_rev);
    lunar_rev_time = interp1([lunar_orbit_statematrix(ind_rev(end-1),2), lunar_orbit_statematrix(ind_rev(end),2)], [lunar_orbit_statematrix(ind_rev(end-1),1), lunar_orbit_statematrix(ind_rev(end),1)], 0);
    
    lunar_backtrack_time_lunanorm(:,1) = repelem(lunar_rev_time, height(height(lunar_orbit_statematrix))) + (lunarbacktrack_statematrix(:,1)-lunarbacktrack_statematrix(1,1));
    lunarbacktrack_statematrix(:,1) = lunarbacktrack_statematrix(:,1)-lunarbacktrack_statematrix(1,1);
    
    %putting backtrack statematrix in geo coord system
    lunar_backtrack_statematrix_geocoords = zeros(height(lunarbacktrack_statematrix),6);
    for n=1:height(lunarbacktrack_statematrix)
        lunar_statematrix_during_backtrack(n,:) = interp_statematrix_timebetween(lunar_orbit_statematrix, lunar_orbit_statematrix(:,1), lunar_backtrack_time_lunanorm(n,1), lunar_backtrack_time_lunanorm(n,1));
        lunar_backtrack_statematrix_geocoords(n,:) = lunarbacktrack_statematrix(n,2:7) + lunar_statematrix_during_backtrack(n,2:7);
        if norm(lunar_backtrack_statematrix_geocoords(n,2:4) - lunar_statematrix_during_backtrack(n,2:4)) < lunar_SOI_radius && n < height(lunarbacktrack_statematrix)-10
            stays_in_SOI = false;
            error("fails to reach back to surface")
            %break
        end

        if norm(lunarbacktrack_statematrix(n,2:4)) - 5e3 < lunar_radius && n ~= height(lunarbacktrack_statematrix)
            avoids_collision = false;
            %break
        end
    end

    lunar_backtrack_statematrix_geocoords = [lunarbacktrack_statematrix(:,1), lunar_backtrack_statematrix_geocoords];
    
    %finding ejection angles
    eject_r_geo = lunarbacktrack_statematrix(end,2:4)+lunar_orbit_statematrix(1,2:4); %ejection from MD
    eject_r_local = lunarbacktrack_statematrix(end,2:4);
    eject_v_local = lunarbacktrack_statematrix(end,5:7);

    [firing_az, firing_el] = find_surface_azel(eject_r_local,eject_v_local); %azimuth works clockwise with Z as reference, eg 90 degrees is lunar east

    eject_lunar_r_geo = lunar_orbit_statematrix(1,2:4);
    eject_lunar_v_geo = lunar_orbit_statematrix(1,5:7);
    eject_lunar_normal = cross(eject_lunar_r_geo./norm(eject_lunar_r_geo),eject_lunar_v_geo./norm(eject_lunar_v_geo));
    
    [MD_lat,MD_long] = find_lunar_latlong(eject_lunar_normal, -eject_lunar_r_geo, eject_r_local)

    if firing_el > 1e-4 && firing_el < 180
        ejection_possible = true;
    end

    total_transit_statematrix = [flipud(lunar_backtrack_statematrix_geocoords(2:end,:)), repelem(1,height(lunar_backtrack_statematrix_geocoords)-1).'; earthtransit_statematrix, repelem(2,height(earthtransit_statematrix)).'];
    
    dv_match = norm(total_transit_statematrix(end,5:7) - target_v); %dv to match, not using the one from the lambert solver
    
    min_acc_req = dv_match/impulse_assumption_time;

else
    MD_lat = nan;
    MD_long = nan;
    firing_az = nan;
    firing_el = nan;
    total_transit_statematrix = repelem(nan,8);
    dv_match = inf;
    min_acc_req = nan;
end

stays_within_agreeable_time = true;
if delta_t > 20*24*3600
    stays_within_agreeable_time = false;
end

max_altitude = 0;
for n=1:20:height(total_transit_statematrix)
    if norm(total_transit_statematrix(n,2:4)) > max_altitude
        max_altitude = norm(total_transit_statematrix(n,2:4));
    end
end
if max_altitude > 1.25e9 %earth SOI get a bit weird past this point
    stays_in_SOI = false;
end

check_matrix = [ %if any of these are low its not a valid traj
reaches_surface
stays_in_SOI
no_SOI_reentry
stays_within_agreeable_time
ejection_possible
avoids_collision
];

return_scalar = [
dv_match
delta_t
firing_az
firing_el
MD_lat
MD_long
];

total_transit_statematrix;

%% plots

%subplot(1,2,1)

hold on
grid on
axis equal padded
ax = gca;
ax.FontSize = 20;
%set(gcf, 'Color', [1,1,1])

set(gcf, 'Color', [0,0,0])
set(ax, 'Color', [0,0,0])

view([20,20])

%plot earth
earth_map = imread('earth_map.jpg');
earth_map = flipud(earth_map);
[surface_map_x,surface_map_y,surface_map_z] = ellipsoid(0,0,0,earth_radius,earth_radius,earth_radius,70);
earth_obj = surf(surface_map_x,surface_map_y,surface_map_z,'EdgeColor','none'); 
set(earth_obj,'cdata',earth_map,'facecolor','texturemap')


%plot lunar SOI
lunar_at_SOIchange = [lunar_statematrix_during_backtrack(1,2:7)];
[surface_map_x,surface_map_y,surface_map_z] = ellipsoid(lunar_at_SOIchange(1),lunar_at_SOIchange(2),lunar_at_SOIchange(3),lunar_SOI_radius,lunar_SOI_radius,lunar_SOI_radius,50);
SOI_obj = surf(surface_map_x,surface_map_y,surface_map_z,FaceColor=[repelem(1,3)],FaceAlpha=0.1,EdgeAlpha=0.1,EdgeColor=[1,1,1]); 

% luna at SOI entry
lunar_map = imread('lunar_map.jpg');
lunar_map = flipud(lunar_map);
[surface_map_x,surface_map_y,surface_map_z] = ellipsoid(lunar_statematrix_during_backtrack(end,2),lunar_statematrix_during_backtrack(end,3),lunar_statematrix_during_backtrack(end,4),lunar_radius,lunar_radius,lunar_radius,80);
lunar_obj_1 = surf(surface_map_x,surface_map_y,surface_map_z,'EdgeColor','none'); 
set(lunar_obj_1,'cdata',lunar_map,'facecolor','texturemap','FaceAlpha',0.5)

%lunar orbit
plot3(lunar_state_matrix(:,1),lunar_state_matrix(:,2),lunar_state_matrix(:,3),"w")


%target orbit
plot3(target_orbit_statematrix(:,2),target_orbit_statematrix(:,3),target_orbit_statematrix(:,4),"w")
scatter3(target_r(1),target_r(2),target_r(3),"w","filled")
target_v_unit = target_v./norm(target_v);
quiver3(target_r(1),target_r(2),target_r(3), target_v_unit(1)*unit_vec_size,target_v_unit(2)*unit_vec_size,target_v_unit(3)*unit_vec_size,"k")

%earth transit
plot3(earthtransit_statematrix(:,2),earthtransit_statematrix(:,3),earthtransit_statematrix(:,4),"r")

%lunar backtrack
plot3(lunar_backtrack_statematrix_geocoords(:,2),lunar_backtrack_statematrix_geocoords(:,3),lunar_backtrack_statematrix_geocoords(:,4),"b",LineWidth=1)
%scatter3(eject_r_true(1),eject_r_true(2),eject_r_true(3),10,"k","filled")


% % %subplot(1,2,2)
% hold on
% grid on
% axis equal
% plot3(lunar_state_matrix(:,1),lunar_state_matrix(:,2),lunar_state_matrix(:,3),"k")
% xlim([r_initial(1)-1e7,r_initial(1)+1e7])
% ylim([r_initial(2)-1e7,r_initial(2)+1e7])
% zlim([r_initial(3)-1e7,r_initial(3)+1e7])
% 
% lunar_map = imread('lunar_map.jpg');
% lunar_map = flipud(lunar_map);
% [surface_map_x,surface_map_y,surface_map_z] = ellipsoid(r_initial(1),r_initial(2),r_initial(3),lunar_radius,lunar_radius,lunar_radius,80);
% lunar_obj_1 = surf(surface_map_x,surface_map_y,surface_map_z,'EdgeColor','none'); 
% set(lunar_obj_1,'cdata',lunar_map,'facecolor','texturemap','FaceAlpha',0.5)
% 
% plot3(lunarbacktrack_statematrix(:,2) + lunar_orbit_statematrix(1,2),lunarbacktrack_statematrix(:,3)+lunar_orbit_statematrix(1,3),lunarbacktrack_statematrix(:,4)+lunar_orbit_statematrix(1,4),"b")
% 
% eject_r_norm = eject_r_local./norm(eject_r_local);
% eject_v_norm = eject_v_local./norm(eject_v_local);
% 
% scatter3(eject_r_geo(1),eject_r_geo(2),eject_r_geo(3),"filled","k")
% quiver3(eject_r_geo(1), eject_r_geo(2), eject_r_geo(3), eject_r_norm(1)*unit_vec_size, eject_r_norm(2)*unit_vec_size, eject_r_norm(3)*unit_vec_size,"k")
% quiver3(eject_r_geo(1), eject_r_geo(2), eject_r_geo(3), eject_v_norm(1)*unit_vec_size, eject_v_norm(2)*unit_vec_size, eject_v_norm(3)*unit_vec_size,"k")
% 
% eject_lunar_r_norm = - eject_lunar_r_geo./norm(eject_lunar_r_geo);
% quiver3(eject_lunar_r_geo(1), eject_lunar_r_geo(2), eject_lunar_r_geo(3), eject_lunar_r_norm(1)*unit_vec_size, eject_lunar_r_norm(2)*unit_vec_size, eject_lunar_r_norm(3)*unit_vec_size,"k")

%view([0,90-inc_lunar])
view([30,30])

set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')

ax.XColor = 'w';
ax.YColor = 'w';
ax.ZColor = 'w';
ax.Title.Color  = 'w';
ax.XLabel.Color = 'w';
ax.YLabel.Color = 'w';
ax.GridColor = [1 1 1];
ax.MinorGridColor = [1 1 1];
ax.GridAlpha = 0.2;


function [lat,long] = find_lunar_latlong(n,a,b)
    n = n./norm(n);
    a = a./norm(a);
    b = b./norm(b);

    long_ref_x = a - dot(a,n)*n;
    long_ref_y = cross(n,long_ref_x);

    if dot(b,n) == 1
        b = b+[repelem(1e-12,3)]; %exactly polar
    end

    lat = asind(dot(b,n));
    long_tmp = atan2d( dot(b,long_ref_y) , dot(b,long_ref_x) );

    long = long_tmp;
end


function [az,el] = find_surface_azel(a,b)

    az_ref = [0,0,1]; % z axis defines 0 degrees asimuth

    a = a./norm(a); %normalise
    b = b./norm(b);

    el = asind(dot(a,b));

    b_perp_a = b - dot(b,a)*a;
    b_perp_a = b_perp_a./norm(b_perp_a);

    x = az_ref - dot(az_ref,a)*a;
    x = x./norm(x);
    x_ref = cross(a,x);

    az = -atan2d(dot(b_perp_a,x_ref),dot(b_perp_a,x)); %clockwise (lunar east) is positive asimuth
    az = rem(az+360,360);
end


function TA_list = TA_from_statematrix(statematrix, normal_vec)
    peri_ref = statematrix(1,2:4);
    if normal_vec(3) < 0
        normal_vec = -normal_vec;
    end
    for n=1:height(statematrix)
        r_spec = statematrix(n,2:4);
        TA_spec = atan2d( normal_vec * cross(peri_ref,r_spec).', dot(peri_ref,r_spec) );
        TA_list(n) = rem(TA_spec+360,360);
    end
end

function closest_inds = find_closest_inds(a,b)
    [~,closest_inds] = mink( abs(a-b), 2);
    closest_inds = sort(closest_inds).';
end

function state_out = interp_statematrix_TA(statematrix, TA)
    between_inds = find_closest_inds(statematrix(:,8), TA);
    for n=1:width(statematrix)
        state_out(n) = interp1([statematrix(between_inds(1),8),statematrix(between_inds(2),8)], [statematrix(between_inds(1),n),statematrix(between_inds(2),n)], TA, "linear","extrap");
    end
end

function state_out = interp_statematrix_timebetween(statematrix, timeseries, time_close, time_target)
    between_inds = find_closest_inds(timeseries, time_close);
    for n=1:width(statematrix)
        state_out(n) = interp1([statematrix(between_inds(1),1),statematrix(between_inds(2),1)], [statematrix(between_inds(1),n),statematrix(between_inds(2),n)], time_target, "linear","extrap");
    end
end

function output = g_dynamics_twobody(t,state,primary_mass)
    G = 6.67e-11;
    velocity = state(4:6);
    r = state(1:3); 
    r_vector = norm(r);
    r_unit_vector = r/r_vector;
    force_gravity = r_unit_vector*(-G*primary_mass/(r_vector^2));
    acceleration = force_gravity;
    output = [velocity; acceleration];
end