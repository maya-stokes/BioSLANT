function [speciesmat,p,g,tmat] = BioSLANTrun(p,g,niter,lem_dir,lem_name,speciesmat,tmat,save_name)
%BioSLANTrun runs BioSLANT
%
% Syntax
%
%     BioSLANTinitialize(p,g,niter,lem_dir, lem_name, speciesmat, tmat,
%     save_name); 
%
% Input arguments
%
%     p     parameter structure
%     g     landscape attributes created from function
%           "Bio_SLANT_Landscape_Initialize.m"
%     niter     number of random iterations
%     lem_dir   directory to find landscape attribute .mat files
%     lem_name  rootname for landscape attribute files
%     speciesmat    matrix with species identities (initial condition)
%     tmat          structure for incipient species
%     save_name basename of file where speciesmat's will be saved
%
% Author: Maya Stokes (mstokes@mit.edu) and Taylor Perron (perron@mit.edu)
% Date: 19 March 2020

p.tr_run = 1; p.run_ic = 0;  
nlandscapes = length(dir([lem_dir,'/*.mat'])); 
[p,g] = setupncm(p,g,niter); %get new habitat capacity, landscape attributes without 


for i = 2:nlandscapes
    g0 = g; 
    speciesmat0 = speciesmat; 
    nsp0 = max(speciesmat0); 
    
    load([lem_dir,'/',lem_name,'_',num2str(i),'_.mat'],'g'); % load the new landscape
    
    [p,g] = setupncm(p,g,niter); %get new habitat capacity, landscape attributes without 
    %setting up a new initial condition
    
    if p.ngen == p.saveInt %run the model 
        [speciesmat,~,tmat] = run_ncm(p,g,niter,speciesmat0,g0,nsp0,tmat,p.ngen);
        save([save_name,'_',num2str(i),'.mat'],'speciesmat','p');  %save the initial condition
   
    elseif p.ngen > p.saveInt %to save more frequently than "ngen"
      
        if rem(p.ngen,p.saveInt) == 1
            msg = 'Ngen must be divisible by the save interval'; 
            error(msg);
        end
        
        nsplit = p.ngen./p.saveInt;
        
        for j = 1:nsplit
            [speciesmat,~,tmat] = run_ncm(p,g,niter,speciesmat0,g0,nsp0,tmat,p.saveInt);
            save([save_name,'_lansdcape',num2str(i),'_saveint',num2str(j),'.mat'],'speciesmat','p');  %save the initial condition
        end        
    end
        
        
end

end