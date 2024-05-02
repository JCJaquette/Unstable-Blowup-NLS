
clear
% close all
SAVE_FILE=0;  

n_modes = 4096;  

% tMax_1 = 0.070;
tMax_1 = 0.073;
tMax = 0.0746;% NLS blowup A300; 0.07443
 

delta = 1e-4; 

delta_1 = 5e-5; 
delta_2 = 1e-6; 
% timeStep =(1/(2*pi))*5e-7;
% timeStep =(1/(2*pi))*5e-6;
% timeStep =1e-7;

timeStep =1e-7;
% timeStep =1e-6;

% timeStep =1e-8; % Messing around
  
tspan_1 = 0:delta_1:tMax_1;
tspan_2 = (tMax_1+delta_2):delta_2:tMax;
tspan = [tspan_1,tspan_2];
% tspan = [0 tMax];  % Messing around
% tspan = 0:delta:tMax;  % Messing around

N_times = length(tspan);

tic, dom = [0 1]; x = chebfun('x',dom); 
 
% 100 Cos  ( NLS Blowup at 0.2435 ) 
% w =         -41.590930582108129  +100*cos(2*pi*x);
 

% 200 Cos  ( NLS Blowup at 0.0834 ) 
% w = -1.117446238072070e+02 +200*cos(2*pi*x);

% 300 Cos  ( NLS Blowup at 0.38 ) 
% w =     -1.892868406016350e+02 +300*cos(2*pi*x);
  

% USED FOR DATA
w =     -1.892868406016350e+02 +300*cos(2*pi*x);
w=w-.5;

S = spinop(dom,tspan);
S.lin    = @(u) (1i)*diff(u,2);
S.nonlin = @(u) (1i)*( u.^2) ; % spin cannot parse "u.*diff(u)"
 
S.init =  w;

% 'exprk5s8' pecec736
figure
tic, u = spin(S,n_modes,timeStep,'dataplot', 'abs','dealias','off','scheme','exprk5s8'); 


if (SAVE_FILE)
    save('blowup_A300_4096_1m7.mat') 
end

return

% % %  Hand cranked Lyapunov Perron method for finding pts on stbl manfiold
% % 
% % 
% % Zero_Coeff_List  = zeros(N_times,1); 
% % One_Coeff_List  = zeros(N_times,1); 
% % for i = 1: N_times
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
% % % -112;;    -0.571877262992029
% % % -111.8;;  -0.133734777517694
% % 
% % %   ;;;  
% %   % ;; 
% % 
% % slope = (( -0.009747398446094) - (-0.003751946505183)) /((-1.8929e+02)-(-1.89288e+02))*1.001;
% % 
% % % slope = ((-2.489765580568670e-11) - (-1.758254772933510e-11)) /((-1.89286840601630e+02)-(-1.892868406016280e+02))*1.001;
% % 
% % estimate  = Zero_Coeff_List(1) -Zero_Coeff_List(end)/slope
 