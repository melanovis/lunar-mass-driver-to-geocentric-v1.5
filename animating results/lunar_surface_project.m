format compact
clear
clc
clf reset

%----------

load("results_formatted.mat")
load("craters_formatted.mat")

unit_vec_size = 2;

hold on
grid on
axis tight equal
set(gcf, 'Color', [1,1,1])

ax = gca;
ax.FontSize = 20;

surfacemap = imread('lunar_map.png');
surfacemap = flipud(surfacemap);

wm = image(surfacemap,'xdata',[-180,180],'ydata',[-90 90]);

set(gca,'ydir','normal')
uistack(wm,'down')

sites = string(fieldnames(pso_results_formatted));

cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));
cmap = flip(cmap);

for n=1:numel(fieldnames(pso_results_formatted))
    tmp_struct = getfield(pso_results_formatted,sites(n));
    coord_list(n,:) = single( [tmp_struct.scalar(5),tmp_struct.scalar(6)] );
    firing_list(n,:) = single( [tmp_struct.scalar(3),tmp_struct.scalar(4)] );
    dv_list(n,1) = single( tmp_struct.scalar(1) );
    eject_v = tmp_struct.vector;
    eject_v_list(n,:) = norm(eject_v(2,:));
    TW_list(n,1) = single( tmp_struct.scalar(7)/9.81 );

    site_name(n,1) = erase(tmp_struct.name,"site ");
end


[~,ind_dvsort] = sort(dv_list);
ind_dvsort = single(ind_dvsort);

MD_centroid = [mean(coord_list(:,2)),mean(coord_list(:,1))];
scene_scale = 60;
x_bounds = [MD_centroid(1)-scene_scale,MD_centroid(1)+scene_scale];
y_bounds = [MD_centroid(2)-scene_scale/2,MD_centroid(2)+scene_scale/2];

xlim(x_bounds)
ylim(y_bounds)
% xlim([-180,180])
% ylim([-90,90])

labels_xlist = single( linspace(min(x_bounds),max(x_bounds)-10,numel(sites)+2) );
labels_xlist([1,end])=[];

labels_ydisp = linspace(min(y_bounds),max(y_bounds),20);
labels_ydisp = [labels_ydisp(5),labels_ydisp(end-4)];

labels_ylist = single( labels_ydisp(single(~iseven(1:numel(sites))+1)) );

[~,ind_longsort] = sort(coord_list(:,2));
ind_longsort = single(ind_longsort.');

coord_x_longsort = coord_list(ind_longsort,2);
coord_y_longsort = coord_list(ind_longsort,1);

dv_to_colour = linspace(min(dv_list),max(dv_list),height(cmap));

for n=1:numel(fieldnames(pso_results_formatted))

    ind_long_spec = n;

    name_spec = "site "+greek_retrofix(char(site_name(n)));
    string_label = string( name_spec + newline ...
        + "$\Delta$v - " + string( round(dv_list(ind_longsort(n))./1e3,3)) +" km/s" + newline ...
        + "T/W - " + round( TW_list(ind_longsort(n)),4)  + newline ...
        + "v$_{MD}$ - " + round(eject_v_list(ind_longsort(n))./1e3,3) + " km/s" + newline ...
        + "$\theta$ - " + round(firing_list(ind_longsort(n),1),1)+ "$^\circ$" + newline + "$\epsilon$ - " + sprintf("%.2e",firing_list(ind_longsort(n),2)) + "$^\circ$" +newline ...
        + "$\phi$ - " + round(coord_y_longsort(n),3)+ "$^\circ$" + newline + "$\lambda$ - " + round(coord_x_longsort(n),3) + "$^\circ$" );

    if labels_ylist(ind_long_spec) < 0
        align = "top";
    else
        align = "bottom";
    end

    plot([coord_list(ind_longsort(n),2), labels_xlist(ind_long_spec)],[coord_list(ind_longsort(n),1), labels_ylist(ind_long_spec)],"w",linewidth=0.1)

    text(labels_xlist(ind_long_spec),labels_ylist(ind_long_spec),string_label, 'Interpreter','latex', Color="w", fontsize=10, verticalalignment = align)
end

for n=1:numel(fieldnames(pso_results_formatted))
    theta = firing_list(n,1);
    eject_vector = [
    cosd(theta), sind(theta)
    -sind(theta), cosd(theta)
    ]*[0;1];
    quiver(coord_list(n,2),coord_list(n,1),eject_vector(1)*unit_vec_size, eject_vector(2)*unit_vec_size,"w",MaxHeadSize=1)
    scatter(coord_list(n,2),coord_list(n,1),120,"kx")
    %text(coord_list(n,2),coord_list(n,1), " "+round(dv_list(n)./1e3,2) ,color="w")
    [~,ind_colour] = min(abs(dv_to_colour-dv_list(n)));
    scatter(coord_list(n,2),coord_list(n,1),20,"filled",markeredgecolor="k",markerfacecolor=[cmap(ind_colour,:)]);
end

xlabel("longitude (degrees)", Interpreter="latex", FontSize=20)
ylabel("latitude (degrees)", Interpreter="latex", FontSize=20)

clim([min(dv_list),max(dv_list)]./1e3)
colormap(cmap)
h = colorbar;
set(get(h,'label'),'string','payload rendezvous $\Delta$v (km/s)', Interpreter="latex", FontSize=20);

set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')


function b = iseven(a)
    b = ~logical(rem(a,2));
end

function b = greek_retrofix(a)
    switch double(a)
        case 913  
            b="A";
        case 914  
            b="B";
        case 915  
            b="$\Gamma$";
        case 916  
            b="$\Delta$";
        case 917  
            b="E";
        case 918  
            b="Z";
        case 919  
            b="H";
        case 920  
            b="$\Theta$";
        case 921  
            b="I";
        case 922  
            b="K";
        case 923  
            b="$\Lambda$";
        case 924  
            b="M";
        case 925  
            b="N";
        case 926  
            b="$\Xi$";
        case 927  
            b="O";
        case 928  
            b="$\Pi$";
        case 929  
            b="P";
        case 931  
            b="$\Sigma$";
        case 932  
            b="T";
        case 933  
            b="$\Upsilon$";
        case 934  
            b="$\Phi$";
        case 935  
            b="X";
        case 936  
            b="$\Psi$";
        case 937  
            b="$\Omega$";
    otherwise
        b = "";
    end
end