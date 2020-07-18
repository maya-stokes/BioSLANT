% Calculate derived timescales to use on axes in the paper. 
% Calls function "CalcDispersalTimescale.m"

ulist = [1 1.2 1.4 1.6 1.8 2 3:2:39]; nu = length(ulist); 
taulist = [1:10:200]; 

theta = 8;
rho = 5; 
ngen = 100; 
nx = 200; 
p_p = 2.5; 
niter = 100; 

% ts = (taulist+1)./theta; %speciation timescale
% td = zeros(nu,1); 
% 
% for i = 1:nu
%  td(i) = CalcDispersalTimescale(rho,theta,p_p,ulist(i),niter,ngen,nx); 
% end
% 
% save('DispersalTimescales.mat','td','ts')
%%

% Supp 1
rho = 4; theta = 4; p_p = 0.9; p_u = 2; tau = 0;  
td_1 = CalcDispersalTimescale(rho,theta,p_p,p_u,niter,ngen,nx); 
ts_1 = (tau+1)./theta; %speciation timescale

% Supp 4 
rho = 5; theta = 8; p_p = 0.4; p_u = 1.1; tau = 10;
td_4 = CalcDispersalTimescale(rho,theta,p_p,p_u,niter,ngen,nx); 
ts_4 = (tau+1)./theta;

% Supp 3
theta_list = [2 4 8 16 32 64]; rho = 8;  p_p = 0.9; p_u = 2; tau = 0; 
for i = 1:length(theta_list)
theta = theta_list(i);
td_3(i) = CalcDispersalTimescale(rho,theta,p_p,p_u,niter,ngen,nx); 
end
ts_3 = (tau+1)./theta_list; 
%%
% Supp 2a
theta_list = [2 4 8 16 32]; p_p = 2.5; 
for i = 1:length(theta_list)
 theta = theta_list(i); 
td_2(i) = CalcDispersalTimescale(rho,theta,p_p,p_u,niter,ngen,nx); 
end
%%
ts_2 = (tau+1)./theta_list;

save('Supp_Derived_Timescales.mat','td_1','td_4','td_3','td_2','ts_1','ts_2','ts_3','ts_4'); 