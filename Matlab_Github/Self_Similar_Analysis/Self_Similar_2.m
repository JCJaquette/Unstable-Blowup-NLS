clear
% close all 
if ~exist('N_times')
    % load('blowup_A300_4096_1m7.mat')
     % load('blowup_A200_4096_1m7.mat')
     load('blowup_Exp300_4096_1m7.mat')
     if imag(u{1}(.125)) > 0
         case_num = 1;
     elseif real(u{1}(0)) < 100
         case_num = 2;
     else
         case_num = 3;
     end
end 

%%%% Integrable case has blowup rate of ~(T-t)^-2 
%%%% Real initial data has blowup rate of ~(T-t)^-1 

 
% close all
format long
DotsPerInch =  300;
N_ss_plots = 500;
SAVE_FILE =1;

TRACK_BLOWUP =1;

RECOMPUTE_BU_POINT = 0;


fig_placement_w = 1200;
fig_placement_h = 600;
fig_width_1 = 450;% Norm
fig_height_1=325;



% There are three cases we consider



if case_num == 3% A =300
    T_blowup = .07443;  % Estimated blowup time ( 0.07443 ); 
    t_start = 0.07;
%     t_start = 0.0736;
    
    % Self-similar scaling
    alpha = 1;
    beta = .5;
    
    xi_max = 15;
    S_MAX = 11;
elseif case_num == 2% A =200
      
    T_blowup_1 =  0.084323116346734;
    T_blowup_2 =  0.084310553658291;
    T_blowup = (T_blowup_1+T_blowup_2)/2; % Estimated blowup time .08433;
    T_blowup = T_blowup_2; 
    t_start = .083;
    
    % Self-similar scaling
    alpha = 1;
    beta = .5;
    % Initial Blowup point to track
%     xi_local = .35;
    xi_max = 15;
    S_MAX = 12;
elseif case_num == 1

    T_blowup = 0.004572; 
    t_start = .000 ;
    
    % T_blowup = 0.00501; 
    
    % Self-similar scaling
    alpha = 2;  
    beta = alpha/2;
    
    % Initial Blowup point to track
    xi_local = .85;
    xi_max = 50;
        % Self-similar scaling
    S_MAX = 11; 
    fig_width_1=fig_width_1+20;
end

% % % Find the index for start of the time window
% % tspan_finding_index =  abs(tspan - t_start);
% %     [~,J_Index] = sort(tspan_finding_index);
% %     J_start = J_Index(1);
% %     J = J_start;
% % tspan_finding_index =  abs(tspan - T_blowup);
% %     [~,J_Index] = sort(tspan_finding_index);
% %     J_max = J_Index(1)-1;
    
[J, J_max ] = find_index(t_start,T_blowup,tspan);

if case_num ==3 % A =300
    r_local = .1591;
    i_local = .1591;
elseif  case_num==2% A =200
    r_local = .63897;
    i_local = .63897;
elseif  case_num==5% A =200
    r_local = .67;
    i_local = .67;
else % case_num ; exp
    r_local = .86692;
    i_local = .86692;
    
    J_max = min([J_max,902]); % Spectral blocking
end
    


tic




if ~exist('realMax_list') || RECOMPUTE_BU_POINT

    realMax_list     = zeros(N_times-J+1,1);
    imagInflect_list = zeros(N_times-J+1,1);

    tic
    for ii = fliplr(1:length(realMax_list))
        u_local=u{J+ii-1};
        disp(ii/length(realMax_list))

        realMax_list_local = roots(diff(real(u_local)));
        imagInflect_list_local = roots((diff(diff(imag(u_local)))));

        rDist =  abs(realMax_list_local - r_local);
        iDist =  abs(imagInflect_list_local - i_local);

        [~,rIndex] = sort(rDist);
        [~,iIndex] = sort(iDist);

        r_local = realMax_list_local(rIndex(1));
        i_local = imagInflect_list_local(iIndex(1));


        realMax_list(ii,1) = r_local;
        imagInflect_list(ii,1) = i_local;

    end
end
toc
disp('Plotting ...')


figure 
set(gcf, 'Position',  [200, 200, 450, 500])
hold on 
plot(tspan(J:N_times),realMax_list)
plot(tspan(J:N_times),imagInflect_list)
xline(T_blowup,'--')
xlabel('t')
hold off
ylabel('x_m(t)')



if (SAVE_FILE)
    title_str = ['Self_similar_no', num2str(case_num),'_blowup_point',  '.png'];
    cd('Figures')
    exportgraphics(gcf,title_str,'Resolution',DotsPerInch);
    cd ..
end
 

% Plot the self similar blowup stuff 


xi_grid_pts = n_modes+1;
if alpha ==1  % Real initial data
    xi = linspace(-xi_max,xi_max,xi_grid_pts);
else                % Exp Integrable
    xi = linspace(-xi_max,xi_max,xi_grid_pts);
end
x_rough = linspace(0,1,xi_grid_pts);


% Dealing with the time stuff
t_max = tspan(J_max-1);
t_min = tspan(J);
s_max = -log(T_blowup-(t_max));
s_min = -log(T_blowup-(t_min));

log_spaced_time_s = linspace(s_min,s_max,N_ss_plots);
log_spaced_time_t = -exp(-log_spaced_time_s)+T_blowup;

SS_times = log_spaced_time_t;
for i = 1:length(SS_times)
    tspan_finding_index =  abs(tspan - log_spaced_time_t(i));
    [~,J_Index] = sort(tspan_finding_index);
    SS_times(i) = J_Index(1);
end
SS_times=unique(SS_times);
 

rescaled_times = -log(T_blowup-tspan(SS_times));

time_list = 0*(1:length(x_rough));

X_list = zeros(length(x_rough),length(SS_times));
Y_list = X_list;
Z_list_real = X_list;
Z_list_imag = X_list;


for component = 0:1 % 0 is real; 1 is imag 
    
    counter = 0;

    for jj = SS_times
        counter = counter +1;
        if component == 0 
            u_local= real(u{jj});
        else
            u_local= imag(u{jj});
        end
        
        ii = jj-J+1;
        
        
%         mid = absMax_list(ii);

        if case_num == 3
            mid =0;
        elseif case_num == 2
            mid = .5;
        elseif case_num == 5
            mid = .5;
        else
            mid=0;
%             mid = (imagInflect_list(ii)+realMax_list(ii))/2;
        end
        if TRACK_BLOWUP 
            mid = (imagInflect_list(ii)+realMax_list(ii))/2;
        end
        
        s= -log(T_blowup-tspan(jj));
        scale_a=(T_blowup-tspan(jj))^alpha;
%         scale_b=(T_blowup-tspan(jj))^beta;
        
        scale_b=((T_blowup-tspan(jj))*s)^beta;
        
        
        time_list= s+0*time_list;
        

        x_grid = xi*scale_b;
        x_grid = x_grid + mid;
        
        
        
        u_rescale_r = scale_a*u_local(x_grid);
        
        X_list(:,counter) = xi;
        Y_list(:,counter) = time_list;
        if component == 0
            Z_list_real(:,counter) = u_rescale_r;
        else
            Z_list_imag(:,counter) = u_rescale_r;
        end
        
    end
   
end

rescaled_times = -log(T_blowup-tspan(SS_times));
s_limits = [rescaled_times(1),min([S_MAX,rescaled_times(end)])];


Z_lim_max = 1.05* max(max(([Z_list_real;Z_list_imag])));
Z_lim_min = 1.05* min(min(([Z_list_real;Z_list_imag])));
Z_limits = [Z_lim_min,Z_lim_max];

for component = 0:1
    
    figure
    set(gcf, 'Position',  [fig_placement_w, fig_placement_h-fig_height_1*component ...
        , fig_width_1, fig_height_1])
    
    if component == 0 
        Z_list = Z_list_real;
        title_str = 'real';
        title(title_str);
    else
        Z_list = Z_list_imag;
        title_str = 'imag';
        title(title_str);
    end
    s_fig = surf(Y_list,X_list,Z_list);
    % shading interp
    c = colorbar;
    % title( [title_str ,' (U)']);
    c.Label.String = [title_str ,'(U)'];
    c.Label.FontSize=11;
    s_fig.EdgeColor = 'none';
    xlim(s_limits); % s variable
    zlim(Z_limits);
    C_limits=max(abs(Z_limits));
    C_limits=[-C_limits,C_limits];
    caxis(C_limits)
    % caxis(Z_limits)
    ylabel('$y$','Interpreter','latex')
    xlabel('$s$','Interpreter','latex')
    
    view(2)
    
    if (SAVE_FILE)
        title_str = ['Self_similar_no', num2str(case_num),'_', title_str,  '.png'];
        cd('Figures')
        exportgraphics(gcf,title_str,'Resolution',DotsPerInch);
        cd ..
    end

end
