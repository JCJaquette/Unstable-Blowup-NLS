clear
close all 
load('blowup_512_1m6_prep.mat') 


DotsPerInch =  400;
fig_placement_w = 500;
fig_placement_h = 300;
fig_width_1 = 930;% Norm
fig_width_2 = 930;% spacetime
fig_width_3 = 1000;% Proportion

fig_height_2=200;
fig_height_3=200;
fig_height_1=200;

figure(10)
set(gcf, 'Position',  [fig_placement_w, fig_placement_h+fig_height_3+fig_height_2, fig_width_1, fig_height_1])
figure(3)
set(gcf, 'Position',  [fig_placement_w, fig_placement_h+fig_height_3, fig_width_2, fig_height_2])
figure(6)
set(gcf, 'Position',  [fig_placement_w, fig_placement_h, fig_width_3, fig_height_3])

SAVE_FILE =1;


% Plots:
% Realtive_L2_modes
% Linfty norm
% space-time abs plot

max_mode_prob = 8;


% 
% Fig 6
figure(10)
set(gcf, 'Position',  [fig_placement_w, fig_placement_h, fig_width_1, fig_height_1])
t_max = 2.70145;

line_width=1;
plot(tspan,1./Linfty_norms,'LineWidth',line_width)


xlabel('$t$','Interpreter','latex')
% ylabel('$\left(\|u^{(N)}(t)\|_{L^\infty}\right)^{-1}$','Interpreter','latex')
% ylabel('$\|u^{(N)}(t)\|_{L^\infty}^{-1}$','Interpreter','latex')
ylabel('$1/\|u(t)\|_{L^\infty}$','Interpreter','latex')
% ylim([0,1.5])
xlim([0,t_max ])
% set(get(gca,'ylabel'),'rotation',0)

if (SAVE_FILE)
    cd('Figures')
    exportgraphics(gcf,'NLS_inv_Linfty.png','Resolution',DotsPerInch);
    cd ..;
end


%%%%%%%%%%%%%%%
% Plot Norm

figure(3)
set(gcf, 'Position',  [fig_placement_w+fig_width_1, fig_placement_h+fig_height_3, fig_width_2, fig_height_2])
x_grid = 0:0.0005:1;
t_max_space_time = 2.5;
for i=1:N_times
    if tspan(i)>t_max_space_time 
        time_index_max = i;
        break
    end
end
space_time = zeros(length(x_grid),time_index_max );
tic
parfor i=1:time_index_max 
    local_u = u{i};
    space_time(:,i)=abs(local_u(x_grid));
end
toc
s= surf(tspan(1:time_index_max),x_grid,space_time,'EdgeColor','none');
xlim([0,t_max_space_time ])
 clim([0,400])

xlabel('$t$','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
c=colorbar;
c.Label.String ='$|u(t,x)|$';
c.Label.Interpreter='latex';  
c.Label.Rotation=0;
view(2)

if (SAVE_FILE)
    cd('Figures')
    exportgraphics(gcf,'space_time_plot.png','Resolution',DotsPerInch);
    cd ..;
end



 

%%%%%%%%%%%%%

% Construct the L^2 proportions.


proportions_L2 = zeros(N_times,max_mode_prob);

for i =1:max_mode_prob
    data = ((abs(Fourier_List(:,i)).^2)./(L2_norms.^2))';
    if i >1
        data = 2*data;
    end
    proportions_L2(:,i)=data;
end


% Fig 6
figure(6)
set(gcf, 'Position',  [fig_placement_w+fig_width_1, fig_placement_h, fig_width_3, fig_height_3])
a=area(tspan,proportions_L2,'LineStyle','none');
a(max_mode_prob).FaceColor = [.25,.875,.8125];
ylim([0,1]);
xlim([0,tspan(end)])
% title([' $|c_i(t)|^2/\|c(t)\|^2$ with $L^2$ norm'],'Interpreter','latex')
% title('relative contribution of each Fourier mode to the $L^2$ norm','Interpreter','latex')
% legend_text = {'$c_0$','$c_1$','$c_2$','$c_3$','$c_4$','$c_5$','$c_6$','$c_7$'};
legend_text = {'$a_0$','$a_1$','$a_2$','$a_3$','$a_4$','$a_5$','$a_6$','$a_7$'};
try 
    l=legend(legend_text ,'Location','southeastoutside','Interpreter','latex',Direction='normal');
catch
    l=legend(legend_text ,'Location','southeastoutside','Interpreter','latex');
end

xlabel('$t$','Interpreter','latex'); 

if (SAVE_FILE)
    cd('Figures')
    exportgraphics(gcf,'Fig_6.png','Resolution',DotsPerInch);
    cd ..;
end


 