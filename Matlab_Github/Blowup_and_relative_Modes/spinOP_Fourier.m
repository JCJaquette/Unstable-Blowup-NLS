
clear
close all
  
% Blowup probably around 2.70145 +/- 1e-4

% dt 1.6e-6 ; dof 512; blow up around 2.7010.
% dt 0.8e-6 ; dof 256; blow up around <2.7010.
% dt 8.0e-6 ; dof 128; blow up around 2.6900.

SAVE_FILE =1;  
n_modes = 512; 
 
tMax = 2.7;% NLS blowup
% tMax = 1;% NLS blowup

delta = 2e-4;  
timeStep =1e-6;
  
tspan = 0:delta:tMax;

N_times = length(tspan);

tic, dom = [0 1]; x = chebfun('x',dom); 

w = -5.3070235  +30*cos(2*pi*x);

S = spinop(dom,tspan);
S.lin    = @(u) (1i)*diff(u,2);
S.nonlin = @(u) (1i)*( u.^2) ; % spin cannot parse "u.*diff(u)"


S.init =  w;

% 'exprk5s8' pecec736
tic, u = spin(S,n_modes,timeStep,'dataplot', 'abs','dealias','on','scheme','exprk5s8'); 


Fourier_List  = zeros(N_times,n_modes/2);
Fourier_List2  = zeros(N_times,n_modes/2);
for i = 1: N_times
    local = trigcoeffs(u{i});
    Fourier_List(i,:) = transpose(local(n_modes/2+1:end));
end

for i = 1: N_times
    local = trigcoeffs(u{i});
    Fourier_List2(i,:) = transpose(local(1:n_modes/2));
end

figure 
hold on
scatter(tspan,abs(Fourier_List(:,1)),'.');
scatter(tspan,abs(Fourier_List(:,2)),'.');
scatter(tspan,abs(Fourier_List(:,3)),'.');
scatter(tspan,abs(Fourier_List(:,4)),'.');
scatter(tspan,abs(Fourier_List(:,5)),'.');
hold off

title('abs mode')



K = 0:(n_modes/2)-1;
L2_norms = zeros(N_times,1);
A1_norms=L2_norms;
H1_Norms = zeros(N_times,1);
H2_Norms = zeros(N_times,1);
for i = 1:N_times
    A1_norms(i) = abs(Fourier_List(i,1))+sum(2*abs(Fourier_List(i,2:end)));
    L2_norms(i) = sqrt(abs(Fourier_List(i,1))^2+2*sum(abs(Fourier_List(i,2:end).^2)));
    H1_Norms(i) = sqrt(2*sum(abs((K.*Fourier_List(i,:)).^2)));
    H2_Norms(i) = sqrt(2*sum(abs(((K.^2).*Fourier_List(i,:)).^2)));
end

figure 
hold on
scatter(tspan,L2_norms);
scatter(tspan,A1_norms);
scatter(tspan,H1_Norms);
scatter(tspan,H2_Norms);
hold off

title('L2 A1, and H1 and H2')

    
    figure 
hold on
scatter(tspan,abs(Fourier_List(:,6))./L2_norms,'.');
scatter(tspan,abs(Fourier_List(:,7))./L2_norms,'.');
scatter(tspan,abs(Fourier_List(:,8))./L2_norms,'.');
scatter(tspan,abs(Fourier_List(:,9))./L2_norms,'.');
scatter(tspan,abs(Fourier_List(:,10))./L2_norms,'.');
hold off

title('Relative abs mode 5-9')


    
    figure 
hold on
scatter(tspan,abs(Fourier_List(:,1))./L2_norms,'.');
scatter(tspan,abs(Fourier_List(:,2))./L2_norms,'.');
scatter(tspan,abs(Fourier_List(:,3))./L2_norms,'.');
scatter(tspan,abs(Fourier_List(:,4))./L2_norms,'.');
scatter(tspan,abs(Fourier_List(:,5))./L2_norms,'.');
hold off

title('Relative abs mode')

%%%%%%%%%% A1 Norms 


    figure 
hold on
scatter(tspan,abs(Fourier_List(:,6))./A1_norms,'.');
scatter(tspan,abs(Fourier_List(:,7))./A1_norms,'.');
scatter(tspan,abs(Fourier_List(:,8))./A1_norms,'.');
scatter(tspan,abs(Fourier_List(:,9))./A1_norms,'.');
scatter(tspan,abs(Fourier_List(:,10))./A1_norms,'.');
hold off

title('Relative abs mode 5-9')


figure
max_mode_prob = 20;
proportions = zeros(N_times,max_mode_prob+1);

other =1+0*A1_norms';
for i =1:max_mode_prob
    data = (abs(Fourier_List(:,i))./A1_norms)';
    if i >1
        data = 2*data;
    end
    proportions(:,i)=data;
    other = other -data;
end
% proportions(:,end)=other ;
bar(tspan,proportions,'stacked')
ylim([0,1]);
title('A1 proportions')



figure

proportions_L2 = zeros(N_times,max_mode_prob+1);

other =1+0*L2_norms.^2';
for i =1:max_mode_prob
    data = ((abs(Fourier_List(:,i)).^2)./(L2_norms.^2))';
%     keyboard
    if i >1
        data = 2*data;
    end
    proportions_L2(:,i)=data;
    other = other -data;
end
% proportions_L2(:,end)=other ;
bar(tspan,proportions_L2,'stacked')
ylim([0,1]);
title('l2 proportions')
    
    figure 
hold on
scatter(tspan,smoothdata(abs(Fourier_List(:,1)))./A1_norms,'.');
scatter(tspan,smoothdata(abs(Fourier_List(:,2)))./A1_norms,'.');
scatter(tspan,abs(Fourier_List(:,3))./A1_norms,'.');
scatter(tspan,abs(Fourier_List(:,4))./A1_norms,'.');
scatter(tspan,abs(Fourier_List(:,5))./A1_norms,'.');
hold off

title('Relative abs mode')


figure 
hold on
scatter(tspan,log(abs(Fourier_List(:,1)))/log(10),'.');
scatter(tspan,log(abs(Fourier_List(:,2)))/log(10),'.');
scatter(tspan,log(abs(L2_norms))/log(10),'.');
xlim([0,.5])
legend('|a_0|','|a_1|','||u||_{L^2}')
hold off

if (SAVE_FILE)
    save('blowup_512_1m6.mat') 
end
