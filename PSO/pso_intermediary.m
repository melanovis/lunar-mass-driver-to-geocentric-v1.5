function [fitness, check_matrix, return_scalar, return_vector, total_transit_statematrix,lunar_r_check] = pso_intermediary(input_vector, lunar_SOI_radius, lunar_orbit_statematrix, latlong_target, target_radius, constrain_long)

    error_tolerance = 1e-10;
    
    MEO_target_inc = 5.145; %specified in sim

    G = 6.6743e-11;
    lunar_mass = 7.35e22; %kg
    earth_mass = 5.972e24;
    lunar_radius = 1.7374e6;
    earth_radius = 6.371e6; %m
    mu = G*earth_mass;
    lunar_orbital_radius = 	363300e3; %m
    lunar_obliquity = 6.687; 
    
    %% params
    
    SOI_ejection_lat = interp1([0,1],[-90,90],input_vector(1)) + 10^interp1([0,1],[-10,-1],input_vector(3)) ; %extra degrees, parameter but not shown
    SOI_ejection_long = interp1([0,1],[-180,180],input_vector(2)) + 10^interp1([0,1],[-10,-1],input_vector(4)) ;
    transfer_dt = input_vector(5); %parameter, normalised
    MEO_target_TA = interp1([0,1],[0,360],input_vector(6)); %degrees, parameter

    if abs(SOI_ejection_lat) > 90
        SOI_ejection_lat = 90*sign(SOI_ejection_lat);
    end
    if abs(SOI_ejection_long) > 180
        SOI_ejection_long = 180*sign(SOI_ejection_long);
    end
    if input_vector(7) < 0.75
        SOI_ejection_lat = MEO_target_inc; %should help with orbital plane selection
    end

    %% ejecting out of lunar SOI
    
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
    
    MEO_target_long_ascend = 0; %degrees parameter
    
    MEO_target_argument_peri = 0;
    
    MEO_target_eccentricity = 0; 
    
    target_semimajor = target_radius/(1-MEO_target_eccentricity);
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
    ]*[0; -target_radius; 0];
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
    ]*[0; -target_radius; 0];
    target_start_r = target_start_r.';
    
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
    target_orbit_timespan = linspace(0,target_orbital_period,1e3);
    
    [target_timerange, target_state_matrix] = ode45(@(timerange, state_matrix)g_dynamics_twobody(timerange, state_matrix,earth_mass),target_orbit_timespan,target_state_initial,odeset('Reltol',error_tolerance));
    target_orbit_statematrix = [target_timerange,target_state_matrix];
    
    TA_list = TA_from_statematrix(target_orbit_statematrix,target_n);
    target_orbit_statematrix = [target_orbit_statematrix, TA_list.'];
    
    target_start_interp = interp_statematrix_TA(target_orbit_statematrix, MEO_target_TA);
    target_r = target_start_interp(2:4);
    target_v = target_start_interp(5:7);
    
    [v_1, v_2, delta_t, ~] = construct_transfer(r_SOIluna_enter, target_r, transfer_dt, mu, target_v); %dt will be a param

    state_initial_earthtransit = [r_SOIluna_enter,v_1];
    orbit_timespan_earthtransit = linspace(0,delta_t*1.1,3.5e3);
    
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
                % re-enters lunar SOI
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
        closest_altitude = min(lunar_altitude_series); %if we miss the surface, how close do we get anyway?
    end
    
    SOI_check_good = true;
    ejection_possible = false;
    avoids_collision = true;
    lat_check = true;
    lunar_r_check = repelem(nan,6);

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
                SOI_check_good = false; % fails to reach back to surface
                break
            end
            if norm(lunarbacktrack_statematrix(n,2:4)) - 100e3 < lunar_radius && n < height(lunarbacktrack_statematrix)-10
                avoids_collision = false; %risk of lunar collision
                break
            end
        end
        
        lunar_r_check = lunar_statematrix_during_backtrack;
        lunar_backtrack_statematrix_geocoords = [lunarbacktrack_statematrix(:,1), lunar_backtrack_statematrix_geocoords];
        
        %finding ejection angles
        eject_r_geo = lunarbacktrack_statematrix(end,2:4)+lunar_orbit_statematrix(1,2:4); %ejection from MD
        eject_r_local = lunarbacktrack_statematrix(end,2:4);
        eject_v_local = lunarbacktrack_statematrix(end,5:7);
    
        [firing_az, firing_el] = find_surface_azel(eject_r_local,eject_v_local); %azimuth works clockwise with Z as reference, eg 90 degrees is lunar east
    
        eject_lunar_r_geo = lunar_orbit_statematrix(1,2:4);
        eject_lunar_v_geo = lunar_orbit_statematrix(1,5:7);
        eject_lunar_normal = cross(eject_lunar_r_geo./norm(eject_lunar_r_geo),eject_lunar_v_geo./norm(eject_lunar_v_geo));
        
        [MD_lat,MD_long] = find_lunar_latlong(eject_lunar_normal, -eject_lunar_r_geo, eject_r_local);
    
        MD_lat = MD_lat + lunar_obliquity*cosd(MD_long); %accounting for lunar obliquity
        if abs(MD_lat) > 80
            lat_check = false;
        end

        if firing_el > 1e-4 && firing_el < 30
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
        eject_r_local = repelem(nan,3);
        eject_v_local = repelem(nan,3);
    end
    
    stays_within_agreeable_time = true;
    if delta_t > 27*24*3600
        stays_within_agreeable_time = false;
    end
    
    max_altitude = 0;
    for n=1:20:height(total_transit_statematrix)
        if norm(total_transit_statematrix(n,2:4)) > max_altitude
            max_altitude = norm(total_transit_statematrix(n,2:4));
        end
    end
    if max_altitude > 1.3e9 %earth SOI get a bit weird past this point
        SOI_check_good = false;
    end
    
    target_proximity = norm(latlong_target-[MD_lat,MD_long]);
    if isnan(target_proximity)
        target_proximity = 360;
    end

    if abs(MD_long) > 90 && constrain_long
        long_check = false;
    else
        long_check = true;
    end

    check_matrix = [ %if any of these are low its not a valid traj
    reaches_surface
    SOI_check_good
    no_SOI_reentry
    stays_within_agreeable_time
    ejection_possible
    lat_check
    long_check
    avoids_collision
    ~isnan(MD_lat)
    ~isnan(MD_long)
    ];
    
    return_scalar = [
    dv_match
    delta_t
    firing_az
    firing_el
    MD_lat
    MD_long
    min_acc_req
    target_proximity
    ];
    
    return_vector = [
    eject_r_local
    eject_v_local
    ];

    if isnan(MD_long)
        MD_long = 180;
    end
    acc_tmp = min_acc_req;
    if isnan(acc_tmp)
        acc_tmp = 10;
    end
    
    %regarding fitness
    fitness = 5e4 / min([dv_match,1e5]);
    fitness = fitness * exp( -0.1 * min([abs(log10(closest_altitude)),10]) );
    fitness = fitness * exp( -0.01 * acc_tmp);

    if ~all(isnan(latlong_target)) %if is nan then no target provided
        fitness = fitness * exp( -0.1 * target_proximity);
        if target_proximity < 0.1
            fitness = fitness*2; %should help lock in precise targeting 
        end
    end
    
    if any(~check_matrix)
        fitness = fitness/2;
    end
    
    %longitude constraint guide
    if abs(MD_long) > 90 && constrain_long
        fitness = fitness * exp( -0.01 * (abs(MD_long)-90) );
    end

    fprintf("-")
end


%-------------------------------------------------


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
    statematrix = real(statematrix);
    time_target = real(time_target);
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