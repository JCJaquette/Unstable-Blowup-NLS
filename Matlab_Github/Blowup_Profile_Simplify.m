% close
clear

case_num =1; % A300 <-> 3;   E300 <-> 1;

if case_num == 3
    load('blowup_A300_4096_1m7.mat')
    case_num =3;
else
    load('blowup_Exp300_4096_1m7.mat')
    case_num =1;
end



%%% TRACK THE MAXIMAL DERIVATIVE OF THE IMAGINARY PART! 
%%% TURN THE GRID ON FOR THE PLOTS

% J=920
% T_blowup = .07425;
% T_blowup = .07425;

%%%% Integrable case has blowup rate of (T-t)^-2 
%%%% Real initial data has blowup rate of (T-t)^-1 


% T_blowup = tspan(end)+10^-6;
close all 
SAVE_FILE = 0; 

DotsPerInch =  300;

fig_placement_w = 1200;
fig_placement_h = 600;
fig_width_1 = 400;% Norm
fig_height_1=325;


% There are three cases we consider

if case_num == 3 % A =300
    
    T_blowup = .07443;  % Estimated blowup time ( 0.07443 );
%     T_blowup = .0746;  % Estimated blowup time ( 0.07443 );
%     T_blowup = tspan(end); 
    t_plot = 0.0743;
    
elseif case_num == 1 % E =300
    T_blowup = .004603; % Estimated blowup time +/- .00002;
    T_blowup = .0047; % Estimated blowup time +/- .0001;
    %     T_blowup = tspan(end);
    J=5; % t= .003
    J_max = 224; % Index for blowup time. 
    
     t_plot = 0.0040;
end
[J , J_max ] = find_index(t_plot ,T_blowup,tspan);
 

% return
% J=4850;
figure
u_local=u{J};
plot(real(u_local))
hold on 
plot(imag((u_local)))
% plot((u_local*conj(u_local))^.5)
legend(['real(u)';'imag(u)'],'Location','north')
% title(['t = ', num2str( tspan(J))],'Interpreter','latex');
% xlim([0.,.25])
hold off
% 
%  

xlabel('$x$','Interpreter','latex')
ylabel('$u$','Interpreter','latex');
 
    set(gcf, 'Position',  [fig_placement_w, fig_placement_h , fig_width_1, fig_height_1])

if (SAVE_FILE)
    title_str = ['Self_similar_no', num2str(case_num),'_blowup_point',  '.png'];
%         exportgraphics(gcf,title_str,'Resolution',DotsPerInch);
    cd('Figures')
    exportgraphics(gcf,title_str,'Resolution',DotsPerInch);
    cd ..
end
 

