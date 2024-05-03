
clear
close all
SAVE_FILE=1;

n_modes = 4096;  
timeStep =1e-7;


% tMax = .0044 ;% 300 exp(1i*2*pi*x)
tMax = .00466;% NLS blowup
 % computation time is about 50 seconds.
 % More time for processing data.

tMax_1 =.0000;


delta = 2.5e-6; 
delta_1 = 1e-3; 
delta_2 = 5e-6; 
delta_2 = max([delta_2 ,timeStep]);

  
tspan_1 = 0:delta_1:tMax_1;
tspan_2 = (tMax_1+delta_2):delta_2:tMax;
tspan = [tspan_1,tspan_2];
% tspan = [0 tMax];
% tspan = 0:delta:tMax;

N_times = length(tspan);

tic, dom = [0 1]; x = chebfun('x',dom); 

w =  300*exp(1i*2*pi*x);


S = spinop(dom,tspan);
S.lin    = @(u) (1i)*diff(u,2);
S.nonlin = @(u) (1i)*( u.^2) ; % spin cannot parse "u.*diff(u)"
S.init =  w;

% 'exprk5s8' pecec736
tic
u = spin(S,n_modes,timeStep,'dataplot', 'real','dealias','on','scheme','exprk5s8'); 
% tic, u = spin(S,n_modes,timeStep,'dataplot', 'abs'); 
toc

if (SAVE_FILE)
    save('blowup_Exp300_4096_1m7.mat')
end
cd .. 
addpath("Blowup_and_relative_Modes\")
process_data
cd("Self_Similar_Analysis\")

if (SAVE_FILE)
    save('blowup_Exp300_4096_1m7_prep.mat')
end

return

% % %  Hand cranked Lyapunov Perron method for finding pts on stbl manfiold
% % 
% % 
% % Zero_Coeff_List  = zeros(N_times,1); 
% % One_Coeff_List  = zeros(N_times,1); 
% % parfor i = 1: N_times
% %     local = trigcoeffs(u{i});
% %     Zero_Coeff_List(i) = transpose(local(n_modes/2+1));
% %     One_Coeff_List(i) = transpose(local(n_modes/2+2));
% % end
% % figure 
% % scatter(tspan(ceil(end/2):end),Zero_Coeff_List(ceil(end/2):end))
% % figure
% % plot(tspan ,log(abs(Zero_Coeff_List))/log(10))
% % Zero_Coeff_List(end)
% % 
% % slope = (( -0.009747398446094) - (-0.003751946505183)) /((-1.8929e+02)-(-1.89288e+02))*1.001;
% % 
% % 
% % estimate  = Zero_Coeff_List(1) -Zero_Coeff_List(end)/slope
% % 
