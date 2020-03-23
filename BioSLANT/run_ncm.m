function [speciesmat,nspecies,t_mat1,varargout] = run_ncm( p,g,niter,speciesmat0,g0,nsp0,t_mat1,nrun)
%run_ncm runs the model! 
%
% Syntax
%
%     run_ncm(p,g,niter,speciesmat0,g0,nsp0,t_mat1,nrun)
%
%
% Description
%
%     run_ncm takes information about the current landscape, the previous
%     landscape, and calls "newfishloop.c" where dispersal, speciation, and
%     random deaths. occur.
%
% Input arguments
%
%     p     parameter structure
%     g     landscape attributes of current landscape
%     niter     number of random iterations
%     speciesmat0   species identity from last landscape
%     g0            previous landscapes attributes; if just running on one
%                   landscape, set g0 = g.
%     nsp0          highest species number from previous landscape
%     t_mat1        structure used for protracted speciation 
%     nrun          number of generations to run "newfishloop" for
%
% Output arguments
%
%     speciesmat    matrix with species identity of every organism 
%     nspecies      highest species number in landscape
%     t_mat1        structure used for protracted speciation
%
%     optional outputs for tracking ancestry (this is still in development)
%
% Author: Maya Stokes (mstokes@mit.edu) and Taylor Perron (perron@mit.edu)
% Date: 19 March 2020

streamOffset = randi(10000); %random seed for the random number seed

nt =  round(p.nfishunits*nrun); %number of iterations to run between changing landscapes
t = round(p.tau*p.nfishunits); %number of iterations required for protracted speciation 


speciesmat = zeros(p.nfishunits,niter); %initialize the outputs
nspecies = zeros(niter,1); 

parfor jj = 1:niter
    
    s = RandStream.create('mrg32k3a','NumStreams',1,'Seed',streamOffset + jj);

    ppp = p; %for the parfor loop (to avoid getting an error because you load these variables
    ggg = g; %from a MAT file above)
    ic = speciesmat0(:,jj);
  
    if ppp.track_ancestry == 1 %if track ancestry initialize child and ancestor vectors
        child  = []; anc = [];
    end
    
    ppp.nspecies = nsp0(jj); %number of species
    t_inc0 = int32(zeros(length(t_mat1(jj).ind),2));
    t_inc0(:,1) = t_mat1(jj).ind; 
    t_inc0(:,2) = t_mat1(jj).val;
     
    if ppp.tr_run == 1 %function for transitioning betweeen landscapes
        [ppp,ggg,t_in, inc_in] = LandscapeTransition( ppp,ggg,g0,ic,t_inc0); %redistribute after landscape change, no speciation
    elseif ppp.tr_run == 0 %the inc_in and t_in are defined in the Redistribute funciton, if not calling assign here
         t_in = int32(t_mat1(jj).val);
         inc_in = int32(t_mat1(jj).ind);
         ggg.speciesvec = speciesmat0(:,jj); 
    end
    
    deathrand = int32(randi(s,ppp.nfishunits,[nt,1])); %fish to die, use independent stream for each worker
    randvec = rand(s,nt,1); %determines if speciation or dispersal occur
    
    if p.track_ancestry == 0
        [speciesmat(:,jj),nsp,inc_out,t_out] = newfishloop_anc_protracted(ggg.speciesvec,ggg.habitatcapacity,...
            ggg.SUMhc,ppp.nspecies,ggg.cdfP,ggg.cdfPind,ggg.fishhab,deathrand,randvec,t,p.track_ancestry,inc_in,t_in); %speciation and dispersal
    else
        [speciesmat(:,jj),nsp,inc_out,t_out,child,anc] = newfishloop_anc_protracted(ggg.speciesvec,ggg.habitatcapacity,...
            ggg.SUMhc,ppp.nspecies,ggg.cdfP,ggg.cdfPind,ggg.fishhab,deathrand,randvec,t,p.track_ancestry,inc_in,t_in); %speciation and dispersal
    end
    
    if any(inc_out)
        t_mat1(jj).ind = int32(inc_out); t_mat1(jj).val = int32(t_out) - nt; %save the incipient values in a structure
        t_mat1(jj).val(t_mat1(jj).val == 0) = -1; %cheating...but don't want any incipient species to be set to zeros (should make it NaN)
    else
        t_mat1(jj).ind = [];  t_mat1(jj).val = [];
    end
    
    nspecies(jj) = nsp; %update species identities

    if p.track_ancestry == 1 %are you tracking ancestry?
        if isfield(ggg,'anc_vic')
            ancestry(jj).child = [ ggg.child_vic; child];
            ancestry(jj).anc = [ggg.anc_vic; anc];
            vicariant_species(jj).list = ggg.child_vic; %if vicariance happens keep track of what species come from vicariance
        elseif ~isfield(ggg,'anc_vic')
            ancestry(jj).child = child;
            ancestry(jj).anc = anc;
            vicariant_species(jj).list = 0;
            
        end
    end
end

if nargout > 4 %if you want to save ancestry information
    varargout{1} = ancestry;
    varargout{2} = vicariant_species;
end

end

