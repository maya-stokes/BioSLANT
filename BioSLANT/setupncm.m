function [p, g, varargout] = setupncm( p,g,niter)
%setupncm calls functions to set-up habitat capacity and dispersal kernels
%
% Syntax
%
%     setupncm(p,g,niter)
%
%
% Description
%
%     take the parameter structure/data structure and outputs 
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
%     p,g
%
%     if p.ic_run == 1:
%     speciesmat    matrix with species identity of every organism 
%     tmat          protracted speciation structure 
%
% Author: Maya Stokes (mstokes@mit.edu) and Taylor Perron (perron@mit.edu)
% Date: 19 March 2020

p.nhabitats = length(g.habitat);
[p,g] = habitatcapacity(p,g); %habitat capacity
p.mu = p.theta./p.nfishunits; %speciation frequency 
g = DispersalKernelTable(p,g); %create dispersal kernel

%set-up the initial condition (speciesmat and tmat for incipient
%speciation)
if p.run_ic == 1
    if p.ic_con == 1 % all the same species to start
        speciesmat = ones(p.nfishunits,niter);
             
    elseif p.ic_con == 0 %all different species to start
        speciesmat = repmat((1:p.nfishunits)',[1 niter]);
        
    end
    
    for j = 1:niter
        tmat(j).ind = []; tmat(j).val = [];
    end
    
    varargout{1} = speciesmat; 
    varargout{2} = tmat; 
end

end

