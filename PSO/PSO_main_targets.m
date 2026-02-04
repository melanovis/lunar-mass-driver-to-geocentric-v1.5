format compact
clear
clc
clf reset

%----------

load("targets_formatted.mat") %target_matrix

target_altitude = 1.2e6; %meters

launch_sites = height(target_matrix);
runs_per_target = 5;

filename = "PSO_result.mat";

if ~exist(filename, 'file')
    ind_run = 1;
    ind_site = 1;
    pso_result_struct = struct();
    save(filename,"ind_run","ind_site","pso_result_struct")
end
load(filename)


ind_run = 1;
while ind_site <= launch_sites
    
    pso_result_struct = setfield(pso_result_struct,"site_"+string(ind_site), "name", target_matrix(ind_site,1) );
    pso_result_struct = setfield(pso_result_struct,"site_"+string(ind_site), "desc", target_matrix(ind_site,4) );

    while ind_run <= runs_per_target
        load(filename);
        
        latlong_target = str2double([target_matrix(ind_site,2:3)]);
        if ind_site == 1
            constrain_long = false;
        else
            constrain_long = true;
        end

        [global_best_fitness, global_best_scalar, global_best_vector, global_best_statematrix, global_best_checkgood, lunar_orbit_statematrix] = PSO_handler(target_altitude, latlong_target, constrain_long);
        
        pso_result_struct = setfield(pso_result_struct,"site_"+string(ind_site),"run_"+string(ind_run),"scalar", global_best_scalar );
        pso_result_struct = setfield(pso_result_struct,"site_"+string(ind_site),"run_"+string(ind_run),"vector", global_best_vector );
        pso_result_struct = setfield(pso_result_struct,"site_"+string(ind_site),"run_"+string(ind_run),"statematrix", global_best_statematrix );
        pso_result_struct = setfield(pso_result_struct,"site_"+string(ind_site),"run_"+string(ind_run),"checkgood", global_best_checkgood );
        pso_result_struct = setfield(pso_result_struct,"site_"+string(ind_site),"run_"+string(ind_run),"lunar_orbit_statematrix", lunar_orbit_statematrix );

        fprintf("completed run %i of site %i.\n",ind_run,ind_site);

        ind_run = ind_run + 1;
        save(filename,"ind_run","ind_site","pso_result_struct", "launch_sites","runs_per_target");
    end
    ind_site = ind_site + 1;
    ind_run = 1;
    save(filename,"ind_run","ind_site","pso_result_struct","launch_sites","runs_per_target");
end

sound(sin(2*pi*400*(0:1/14400:0.15)), 14400);
