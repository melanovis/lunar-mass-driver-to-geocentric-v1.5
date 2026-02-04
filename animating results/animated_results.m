format compact
clear
clc
clf reset

%----------

load("results_formatted.mat")

record_video = true;

luna_following = true;

animation_length = 64; %secs
fps = 60;

unit_vec_size = 1e7;

if record_video
    v = VideoWriter("record_results", 'MPEG-4');
    v.FrameRate = 60;
    open(v);
end

earth_mass = 5.972e24;
earth_radius = 6.371e6; %m
G = 6.6743e-11;
mu = G*earth_mass;
error_tolerance = 1e-12;
lunar_SOI_radius = 5.7922e7;
lunar_radius = 1.7374e6;

%making target orbit
LEO = 1.2e6 + earth_radius;
MEO_target_inc = 5.145;
MEO_target_long_ascend = 0; 
MEO_target_argument_peri = 0;
MEO_target_eccentricity = 0; 
target_semimajor = LEO/(1-MEO_target_eccentricity);
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
]*[0; -LEO; 0];
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

target_state_initial = [target_start_r,target_start_v];
target_orbit_timespan = linspace(0,target_orbital_period,2e3);
[target_timerange, target_state_matrix] = ode45(@(timerange, state_matrix)g_dynamics_twobody(timerange, state_matrix,earth_mass),target_orbit_timespan,target_state_initial,odeset('Reltol',error_tolerance));
target_orbit_statemat = [target_timerange,target_state_matrix];

min_bounds = repelem(inf,4);
max_bounds = repelem(-inf,4);

sites = string(fieldnames(pso_results_formatted));
for n=1:numel(fieldnames(pso_results_formatted))
    tmp_struct = getfield(pso_results_formatted,sites(n));
    if n==1
        lunar_statemat = tmp_struct.lunar_orbit_statematrix;
    end
    tmp_statemat = tmp_struct.statematrix;
    for m=1:4
        min_bounds(m) = min(min_bounds(m),min(tmp_statemat(:,m)));
        max_bounds(m) = max(max_bounds(m),max(tmp_statemat(:,m)));
    end
end
for m=1:4
    min_bounds(m) = min(min_bounds(m),min(lunar_statemat(:,m)));
    if m~=1
        max_bounds(m) = max(max_bounds(m),max(lunar_statemat(:,m)));
    end
end
min_bounds(2:end) = min_bounds(2:end)*1.025;
max_bounds(2:end) = max_bounds(2:end)*1.025;
min_bounds(4) = min_bounds(4)*1.5;
max_bounds(4) = max_bounds(4)*1.5;

timeseries = linspace(min_bounds(1),max_bounds(1),animation_length*fps);
timeseries = sort(unique([timeseries,0,min(timeseries)-3600,max(timeseries)+3600]));

%correcting lunar statematrix
[~,ind_rev] = mink(abs(lunar_statemat(:,2)),5);
ind_rev = sort(ind_rev);
lunar_rev_time = interp1([lunar_statemat(ind_rev(end-1),2), lunar_statemat(ind_rev(end),2)], [lunar_statemat(ind_rev(end-1),1), lunar_statemat(ind_rev(end),1)], 0);

%making lunar backtracked statematrix
[~,inds_cycle] = mink(abs(lunar_statemat(:,1)-lunar_rev_time),2);
[lunar_state_finalcycle,~] = interp_statemat_time(lunar_statemat, lunar_rev_time);
lunar_statemat(inds_cycle:end,:) = [];
lunar_statemat = [lunar_statemat; lunar_state_finalcycle];
lunar_orbit_backtracked = [-(lunar_rev_time - lunar_statemat(:,1)), lunar_statemat(:,2:7)];
lunar_statemat = [lunar_orbit_backtracked(1:end-1,:); lunar_statemat];
lunar_statemat_original = lunar_statemat;
for n=1:2
    forward_time = max(lunar_statemat(:,1)) + 1e-10 + lunar_statemat_original(:,1) + abs(min(lunar_statemat_original(:,1)));
    additional_lunar_statematrix = [forward_time, lunar_statemat_original(:,2:7)];
    lunar_statemat = [lunar_statemat; additional_lunar_statematrix];
end


hold on
grid on
axis tight equal
xlim([min_bounds(2),max_bounds(2)])
ylim([min_bounds(3),max_bounds(3)])
zlim([min_bounds(4),max_bounds(4)])
%view([30,15])
view([0,90])
ax = gca;
ax.FontSize = 20;
set(gcf, 'Color', [0,0,0])
set(ax, 'Color', [0,0,0])
ax.XColor = 'w';
ax.YColor = 'w';
ax.ZColor = 'w';
ax.Title.Color  = 'w';
ax.XLabel.Color = 'w';
ax.YLabel.Color = 'w';
ax.GridColor = [1 1 1];
ax.MinorGridColor = [1 1 1];
ax.GridAlpha = 0.2;

xlabel("x (m)", Interpreter="latex", FontSize=20)
ylabel("y (m)", Interpreter="latex", FontSize=20)
zlabel("z (m)", Interpreter="latex", FontSize=20)

%plot earth
for n=1:2
    scatter(nan,nan,"k");
end
[surface_map_x,surface_map_y,surface_map_z] = ellipsoid(0,0,0,earth_radius,earth_radius,earth_radius,64);
lunar_obj = surf(surface_map_x,surface_map_y,surface_map_z,EdgeColor="none",FaceColor=[0,0,1]);

plot3(lunar_statemat(:,2),lunar_statemat(:,3),lunar_statemat(:,4),"w")
plot3(target_orbit_statemat(:,2),target_orbit_statemat(:,3),target_orbit_statemat(:,4),"w")

MD_fired = logical(zeros([1,numel(fieldnames(pso_results_formatted)) ]));

for ind_t = 1:numel(timeseries)/4

    ind_t
    
    transform_obj = hgtransform;

    time_spec = timeseries(ind_t);

    %where is luna?
    [lunar_statespec,~] = interp_statemat_time(lunar_statemat, time_spec);

    %lunar SOI
    %plot luna and SOI surfs

    lunar_r_tmp = [lunar_statespec(1,2:4)];
    if ~luna_following
        [surface_map_x,surface_map_y,surface_map_z] = ellipsoid(lunar_r_tmp(1),lunar_r_tmp(2),lunar_r_tmp(3),lunar_SOI_radius,lunar_SOI_radius,lunar_SOI_radius,64);
        SOI_obj = surf(surface_map_x,surface_map_y,surface_map_z,FaceColor=[repelem(1,3)],FaceAlpha=0.1,EdgeAlpha=0.01,EdgeColor=[1,1,1],parent=transform_obj); 
    else
        SOI_obj = plot(nsidedpoly(1000, 'Center', [lunar_r_tmp(1),lunar_r_tmp(2)], 'Radius', lunar_SOI_radius),facecolor = [1,1,1],FaceAlpha=0.1,EdgeAlpha=0.01,parent = transform_obj);
    end

    [surface_map_x,surface_map_y,surface_map_z] = ellipsoid(lunar_r_tmp(1),lunar_r_tmp(2),lunar_r_tmp(3),lunar_radius,lunar_radius,lunar_radius,64);
    lunar_obj = surf(surface_map_x,surface_map_y,surface_map_z,EdgeColor="none",FaceColor=[repelem(0.5,3)],parent=transform_obj);
    
    if ~luna_following
        p_1 = text(lunar_statespec(2),lunar_statespec(3),lunar_statespec(4)," luna",Color="w",FontSize=5,parent=transform_obj);
    else
        p_1 = [];
    end

    p_2 = [];
    p_3 = [];
    p_4 = [];
    p_5 = [];
    burn_vector = [];
    for ind_site = 1:numel(fieldnames(pso_results_formatted))
        tmp_struct = getfield(pso_results_formatted,sites(ind_site));
        tmp_statematrix = tmp_struct.statematrix;
        tmp_name = tmp_struct.name;
        if time_spec >= min(tmp_statematrix(:,1)) && time_spec <= max(tmp_statematrix(:,1))
            %good to plot this object
            [craft_statemat, ~] = interp_statemat_time(tmp_statematrix, time_spec);
            
            [~, between_inds_next] = interp_statemat_time(tmp_statematrix, timeseries(min([ind_t+1,numel(timeseries)])) );

            plot_burnvector = false;
            plot_MDvector = false;

            v_vector = repelem(0,3);
            if max(between_inds_next) >= height(tmp_statematrix)
                plot_burnvector = true;

                %need to find station recieving vector
                inds_sample = findclosest_statemat_bisection(target_orbit_statemat, craft_statemat(2:4));
                station_v = target_orbit_statemat(max(inds_sample),5:7);
                v_vector = station_v-craft_statemat(5:7);
                v_vector = v_vector./norm(v_vector);
            end
            if ~MD_fired(ind_site)
                plot_MDvector = true;
                v_vector = craft_statemat(5:7)./norm(craft_statemat(5:7));
                MD_fired(ind_site)=true;
            end

            %paylod
            p_2(ind_site) = scatter3(craft_statemat(2),craft_statemat(3),craft_statemat(4),5,"w","filled",parent=transform_obj);

            [~,ind_tail_end] = min(abs(tmp_statematrix(:,1) - (time_spec - 24*3600 ) ));
            [~,inds_tail_start] = mink( abs(tmp_statematrix(:,1)-time_spec), 2);
            fade_tail_statematrix = [tmp_statematrix(ind_tail_end:min(inds_tail_start),:); craft_statemat];

            %payload fade tail
            p_3(ind_site) =  patch([fade_tail_statematrix(:,2); nan],[fade_tail_statematrix(:,3); nan],[fade_tail_statematrix(:,4); nan],'white','EdgeColor','white',...
                'FaceVertexAlphaData',linspace(0,0.5,height(fade_tail_statematrix)+1).','AlphaDataMapping','none',...
                'EdgeAlpha','interp',parent=transform_obj);

            %elevation line
            p_4(ind_site) = plot3([craft_statemat(2),craft_statemat(2)],[craft_statemat(3),craft_statemat(3)],[craft_statemat(4),min_bounds(4)],'--',LineWidth=0.1,color=[1,1,1,0.1],parent=transform_obj);
        
            %labels
            p_5(ind_site) = text(craft_statemat(2),craft_statemat(3),craft_statemat(4)," MD "+ tmp_name,Color="w",FontSize=5,parent=transform_obj);

            if plot_MDvector || plot_burnvector
                burn_vector(ind_site) = quiver3(craft_statemat(2),craft_statemat(3),craft_statemat(4), v_vector(1)*unit_vec_size, v_vector(2)*unit_vec_size, v_vector(3)*unit_vec_size,"r",parent=transform_obj);
            end
        end
    end

    if time_spec >= 0
        string_additional = "t+";
    else
        string_additional = "t";
    end
    legend("- 1200km altitude rendezvous target (lunar OP) - single burn", "- "+string_additional + string( round( time_spec/3600 ,2)) + " hours SOI departure", Interpreter="latex", FontSize=16 ,location="northwest",textcolor="w")
    legend boxoff

    set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')

    xlim([lunar_r_tmp(1)-1.778e8,lunar_r_tmp(1)+1.778e8])
    ylim([lunar_r_tmp(2)-1e8,lunar_r_tmp(2)+1e8])

    drawnow()

    if record_video
        frame = getframe(gcf);
        writeVideo(v,frame)
    end

    if ind_t ~= numel(timeseries)
        delete(p_1);
        delete(p_2);
        delete(p_3);
        delete(p_4);
        delete(p_5);
        delete(SOI_obj);
        delete(lunar_obj);
        delete(burn_vector);
    end
end

if record_video
    close(v);
end

sound(sin(2*pi*400*(0:1/14400:0.15)), 14400);

function inds_sample = findclosest_statemat_bisection(statematrix, target_r)
    inds_sample = [1,height(statematrix)];
    
    for n=1:100
        s_1 = norm(target_r - statematrix(inds_sample(1), 2:4));
        s_2 = norm(target_r - statematrix(inds_sample(2), 2:4));
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
end

function [statemat_out,between_inds] = interp_statemat_time(statemat_in,time)
    between_inds = find_closest_inds(statemat_in(:,1), time);
    for n=1:width(statemat_in)
        statemat_out(n) = interp1([statemat_in(between_inds(1),1),statemat_in(between_inds(2),1)], [statemat_in(between_inds(1),n),statemat_in(between_inds(2),n)], time, "linear","extrap");
    end
end

function closest_inds = find_closest_inds(a,b)
    [~,closest_inds] = mink( abs(a-b), 2);
    closest_inds = sort(closest_inds).';
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
