format compact
clear
clc
clf reset

%----------

%rejecting bad results

load("PSO_result.mat")

greeks = [
"Α"
"Β"
"Γ"
"Δ"
"Ε"
"Ζ"
"Η"
"Θ"
"Ι"
"Κ"
"Λ"
"Μ"
"Ν"
"Ξ"
"Ο"
"Π"
"Ρ"
"Σ"
"Τ"
"Υ"
"Φ"
"Χ"
"Ψ"
"Ω"  
];

pso_results_formatted = struct();

sites = string(fieldnames(pso_result_struct));
for ind_site = 1:numel(fieldnames(pso_result_struct))
    runs = string( fieldnames( getfield(pso_result_struct,sites(ind_site)) ) );

    result_mat = [];
    for ind_run = 1:numel(runs)
        tmp_struct = getfield(pso_result_struct,sites(ind_site),runs(ind_run));
        run_dv = tmp_struct.scalar(1);
        run_long = tmp_struct.scalar(6);
        run_check = tmp_struct.checkgood;
        result_mat = [result_mat; ind_run,run_dv];
    end
    [~,ind_dv] = min(result_mat(:,2));

    pso_results_formatted = setfield(pso_results_formatted, sites(ind_site), getfield(pso_result_struct,sites(ind_site),runs(ind_dv)) );
    pso_results_formatted = setfield(pso_results_formatted, sites(ind_site), "name", "site "+greeks(ind_site));
end

save("results_formatted.mat","pso_results_formatted")