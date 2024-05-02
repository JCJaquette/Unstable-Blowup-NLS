clear
case_num =3; % A300 <-> 3;   E300 <-> 1;

if case_num == 3
    load('blowup_A300_4096_1m7_prep.mat')
else
    load('blowup_Exp300_4096_1m7_prep.mat')
end
    
SAVE_FILE=0;

DotsPerInch =  300;

fig_placement_w = 1200;
fig_placement_h = 600;
fig_width_1 = 400;% Norm
fig_height_1=325;


if case_num == 3% A =300
    T_blowup = .07443;  % Estimated blowup time ( 0.07443 ); 
    t_start = 0.07;
elseif case_num == 1
    T_blowup = 0.00455; 
    t_start = .000 ;
end

[J, J_max ] = find_index(t_start,T_blowup,tspan);

figure(1)
set(gcf, 'Position',  [fig_placement_w, fig_placement_h , fig_width_1, fig_height_1]);
plot(tspan,1./Linfty_norms)

xlabel('$t$','Interpreter','latex')
ylabel('$1/\|u(t)\|_{L^\infty}$','Interpreter','latex')


xlim([t_start,tspan(end)])
% set(get(gca,'ylabel'),'rotation',0)

if (SAVE_FILE)
    cd('Figures')
    if case_num == 1;
        exportgraphics(gcf,'E300_inv_Linfty.png','Resolution',DotsPerInch);
    else
        exportgraphics(gcf,'A300_inv_Linfty.png','Resolution',DotsPerInch);
    end
    cd ..
end