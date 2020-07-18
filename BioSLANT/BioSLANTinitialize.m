function [p,g,speciesmat,tmat] = BioSLANTinitialize(p,g,niter,save_name)
%BioSLANTinitialize creates initial conditions for the neutral model 
%
% Syntax
%
%     BioSLANTinitialize(p,g,niter,save_name)
%
%
% Description
%
%     BioSLANTinitialize takes a parameter structure, a landscape
%     attribute structure to create/spin-up an initial condition for the
%     model fully coupled to the changing landscapes. 
%
% Input arguments
%
%     p     parameter structure
%     g     landscape attributes created from function
%           "Bio_SLANT_Landscape_Initialize.m"
%     niter     number of random iterations
%     save_name basename of file where the resulting initial condition will
%               be saved
%
% Output arguments
%
%     speciesmat    matrix with species identity of every organism 
%
% Author: Maya Stokes (mstokes@mit.edu) and Taylor Perron (perron@mit.edu)
% Date: 19 March 2020



p.tr_run = 0; p.run_ic = 1; 

[p,g,speciesmat,tmat] = setupncm(p,g,niter); %set-up the initial condition
nsp = max(speciesmat);
nsplit = ceil(p.ss_gen/p.saveInt); %save frequency
for j = 1:nsplit
    [speciesmat,nsp,tmat] = run_ncm(p,g,niter,speciesmat,g,nsp,tmat,p.saveInt); %run the model 
    save([save_name,'_',num2str(j),'_IC.mat'],'speciesmat','tmat','p');  %save the initial condition
end

end