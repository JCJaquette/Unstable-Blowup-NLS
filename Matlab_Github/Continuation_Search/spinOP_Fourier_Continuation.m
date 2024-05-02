
clear
% close all
% mu = 1i+.02;
% mu = mu/abs(mu);
% mu=exp(pi*1i*(.49));
mu=1i;
n_modes = 256; 
tMax = 1;
delta = 2e-4; 
timeStep =1e-4;
tspan = 0:delta:tMax;  % Messing around

DotsPerInch=300;
SAVE_FIGURE = 1;

N_times = length(tspan);

tic, dom = [0 1]; x = chebfun('x',dom); 
% Define constant for trapping region.
Amp_cos = 30;
C2 = exp(-pi/2);
w =Amp_cos*cos(2*pi*x);
A_max = 150;

A_delta = 1;

% A_list = -440+( -A_max:A_delta:A_max);
% -1.117446238072070e+02
% A_list = -115:A_delta:-106;
% A_list = -1.892868406016300e+02

A_list = ( -A_max:A_delta:A_max);
% A_max=2;
% A_list = -5+( -A_max:.05:A_max);
a_length = length(A_list);
N_times = length(tspan);

solution_norms = zeros(a_length,N_times);
trapping_times =NaN+ zeros(a_length,1);
trapping_rat = zeros(a_length,N_times);


parfor i = 1:a_length
    i;
    A_local = A_list(i);
    w_local = A_local + w;
    
    S = spinop(dom,tspan);
    S.lin    = @(u) (mu)*diff(u,2);
    S.nonlin = @(u) (mu)*( u.^2) ; % spin cannot parse "u.*diff(u)"
    S.init =  w_local;
    try
        tic;
%         u = spin(S,n_modes,timeStep,'dataplot', 'abs','dealias','on');
        u = spin(S,n_modes,timeStep,'plot', 'off','dealias','on');
    %     u = spin(S,n_modes,timeStep);
        toc;
    catch
        disp('solution blew up')
        continue
    end
        for j = 1:N_times
            solution_norms(i,j) = sqrt(max(u{j}*conj(u{j})));
        end

        % Find the trapping region
        for j = 1: N_times
            trap_out = check_trapping( u{j} , n_modes);
            TRAPPED = trap_out(1);
            trapping_rat(i,j) = trap_out(2)/(trap_out(3)*C2);
            if TRAPPED ==1
                trapping_times(i) =tspan(j);
                break
            end
        end
        if TRAPPED ==0
            trapping_times(i) =2*tMax;
            % trapping_times(i) =tspan(j);
        end

end


normalized_sol =solution_norms;
for i = 1:a_length
    normalized_sol(i,:) = normalized_sol(i,:)/normalized_sol(i,1);
end


disp('Plotting')

fig_placement_w = 1000;
fig_placement_h = 300;
fig_width = 700;
fig_height=350;

figure(10)
set(gcf, 'Position',  [fig_placement_w, fig_placement_h, fig_width, fig_height])

s=surf(tspan,A_list,solution_norms)
view(2)
colorbar
s.EdgeColor = 'none';
c=colorbar;

xlabel('$t$','Interpreter','latex')
ylabel('$A$','Interpreter','latex')
c.Label.String ='$||u(t)||_{L^\infty}$';
c.Label.Interpreter='latex';  
c.Label.Rotation=0;
c.Label.FontSize=11;
z_max = max(max(solution_norms));
hold on
line_color = [.825,.825,.825];
plot3(trapping_times,A_list,z_max+0*A_list,"Color",line_color ,'LineWidth',3,'LineStyle','--')
hold off

xlim([0,tMax])
set(get(gca,'ylabel'),'rotation',0)

if SAVE_FIGURE 
    title_str='WideSearch_final.png';
    exportgraphics(gcf,title_str,'Resolution',DotsPerInch);
end


return



function  output = check_trapping( u_local ,n_modes)
    F_coeff = trigcoeffs(u_local);
    a0 = abs(F_coeff(n_modes/2+1));
    F_coeff(n_modes/2+1) =0;
    diff_norm = sum(abs(F_coeff));
%     
% \frac{\text{asum} e^{\frac{\pi  \text{a0}
%    r}{2}}}{\text{a0}^2}<r
% 
% \frac{\text{asum}}{\text{a0}}<e^{-\frac{\pi 
%    R}{2}} R
% 
    if diff_norm < a0 *(2/(exp(1)*pi))
        TRAPPED =1;
    else 
        TRAPPED =0;
    end
    output = [TRAPPED, diff_norm, a0];
end



% % % % % Log
% % % % figure
% % % % s=surf(tspan,A_list,log(normalized_sol)/log(10))
% % % % view(2)
% % % % s.EdgeColor = 'none';
% % % % c=colorbar;
% % % % 
% % % % xlabel('t')
% % % % ylabel('A')
% % % % c.Label.String ='||u(t)||_{L^\infty}';
% trapping_rat_2 = min(trapping_rat,10+0*trapping_rat);
% figure
% s2=surf(tspan,A_list,trapping_rat_2)
% s2.EdgeColor = 'none';
% colorbar
