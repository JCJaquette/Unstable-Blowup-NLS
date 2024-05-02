clear;
% close all;

tic
mu = 1i;
T_Max = 120; 
delta_t = .005;
time_input = 0:delta_t :T_Max;

DotsPerInch =  300;

SAVE_FILE = 1;
% if SAVE_FILE
%     cd('Figures')
% end

N = chebop(0,T_Max );
% Initial Data: 

% Tuning
% A0 = -0.22301409257004942;
% A1 = 0.5;
% A2 = 0;
% A3 = 0;

A0 = -0.22007809004711737;
A1 = 0.47826583586743143;
A2 = 0.0917468218790068;
A3 = 0.012948006123315526;
% 
% A0 =-0.063744023333489;
% A1 =0;
% A2 =0.5079523381117685;
% A3 =0;

% A0 =-0.027885343504363874;
% A1 =0;
% A2 =0;
% A3 =0.5015485873004775;

N.op = @(t,c0,c1,c2,c3) [diff(c0)-  mu*(c0^2 +2*(c1^2+c2^2+c3^2) );
                        diff(c1)-   mu*(-c1         +2*(c0*c1+c1*c2+c3*c2));
                        diff(c2)-   mu*(-4*c2+2*c0*c2+c1^2+c3*c1);
                        diff(c3)-   mu*(-9*c3+2*c0*c3+2*c1*c2)];
N.lbc = [A0;A1;A2;A3];
[c0,c1,c2,c3] = N\0;

toc

disp('Graphing ... ')
tic

%  Plotting L2 Norms
norms = (c0*conj(c0))+ 2*((c1*conj(c1))+(c2*conj(c2))+(c3*conj(c3)));
figure


% 
% 
% if (SAVE_FILE)
%     exportgraphics(gcf,'ToyModel_L2.png');
% end

%  Plotting L\infty Norms
dom = [0 2*pi]; x = chebfun('x',dom); 

Linfty_norms=0*time_input;
disp('Computing L_{\infty} Norms')
tic
parfor i = 1:length(time_input)
    t=time_input(i);
    local_u=c0(t)+2*c1(t)*cos(x)+2*c2(t)*cos(2*x)+2*c3(t)*cos(3*x);
    Linfty_norms(i) = sqrt(max(local_u*conj(local_u)));
end
toc

hold on 
plot(time_input,Linfty_norms)
plot(time_input,sqrt(norms(time_input)))
xlabel('t','Interpreter','latex')
hold off 



disp('Computing Proportions')

%%% Relative Modes 


rel_0= abs(c0).^2./norms ;
rel_1= abs(c1).^2./norms ;
rel_2= abs(c2).^2./norms ;
rel_3= abs(c3).^2./norms ;


proportions_L2 = zeros(length(time_input),4);
proportions_L2(:,1) = rel_0(time_input);
proportions_L2(:,2) = 2*rel_1(time_input);
proportions_L2(:,3) = 2*rel_2(time_input);
proportions_L2(:,4) = 2*rel_3(time_input);


disp('Graphing')

fig_placement_w = 1000;
fig_placement_h = 300;
fig_width_1 = 400;
fig_width_2 = 500;
fig_height=350;

line_width=1;
figure(10)
plot(time_input,1./Linfty_norms,'LineWidth',line_width)
xlabel('$t$','Interpreter','latex')
% ylabel('$\left(\|u^{(N)}(t)\|_{L^\infty}\right)^{-1}$','Interpreter','latex')
% ylabel('$\|u^{(N)}(t)\|_{L^\infty}^{-1}$','Interpreter','latex')
% ylabel('$\|u(t)\|_{L^\infty}^{-1}$','Interpreter','latex')
ylabel('$1/\|u(t)\|_{L^\infty}$','Interpreter','latex')
ylim([0,1.5])
xlim([0,time_input(end)])
% set(get(gca,'ylabel'),'rotation',0)
set(gcf, 'Position',  [fig_placement_w, fig_placement_h, fig_width_1, fig_height])
if (SAVE_FILE)
    cd('Figures')
    exportgraphics(gcf,'ToyModel_inv_Linfty.png','Resolution',DotsPerInch);
    cd ..
end
% return

figure(11)
a=area(time_input,proportions_L2,'LineStyle','none');
set(gcf, 'Position',  [fig_placement_w+fig_width_1, fig_placement_h, fig_width_2, fig_height])
ylim([0,1]);
xlim([0,time_input(end)])
% title([' $|c_i(t)|^2/\|c(t)\|^2$ with $L^2$ norm'],'Interpreter','latex')
% title('relative contribution of each Fourier mode to the $L^2$ norm','Interpreter','latex')
legend_text = {'$a_0$','$a_1$','$a_2$','$a_3$'};
try 
    legend(legend_text ,'Location','eastoutside','Interpreter','latex','Direction','normal')
catch
    legend(legend_text ,'Location','eastoutside','Interpreter','latex')
end
xlabel('$t$','Interpreter','latex')

if (SAVE_FILE)
    cd('Figures')
    exportgraphics(gcf,'ToyModel_Relative_Modes.png','Resolution',DotsPerInch);
    cd ..
end

toc  

% 
% if SAVE_FILE
%     cd('..')
% end