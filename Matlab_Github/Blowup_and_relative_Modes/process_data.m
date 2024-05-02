% close all 
% clear
%%% Load the data file; ~~ 200MB to 1.1GB and depends on the Chebfun libraries
% load('blowup_512_1m6.mat') 
% load('blowup_A300_4096_1m7.mat')
% load('blowup_Exp300_4096_1m7.mat')



SAVE_FILE =1;


% The function is given by a cosine series where a_k = a_{-k}. 
Fourier_List  = zeros(N_times,n_modes/2); % Coefficients a_0, a_1,  ... a_N
parfor i = 1: N_times
    local = trigcoeffs(u{i});
    Fourier_List(i,:) = transpose(local(n_modes/2+1:end));
end


% Computes various norms at each time frame.

% Note: We only store what the solution looks like (ie its Fourier coefficients) 
%       once every 50 timesteps of the integrator. 


Linfty_norms = zeros(N_times,1);
L2_norms =  zeros(N_times,1);

parfor i = 1:N_times
    L2_norms(i) = sqrt(abs(Fourier_List(i,1))^2+2*sum(abs(Fourier_List(i,2:end).^2)));
    local_u = u{i};
    Linfty_norms(i) = sqrt(max(local_u*conj(local_u)));
end
    
invLinfty = (1./Linfty_norms);
invL2 = (1./L2_norms);
invZeroMode = 1./abs((Fourier_List(:,1)));

if (SAVE_FILE) && n_modes == 512
    save('blowup_512_1m6_prep.mat') 
end

