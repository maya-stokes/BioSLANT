%Example: from Tadpole output to fully coupled run.

% --------------- Filenames ----------------------------------------------%

tadpole_output = 'example_fault.mat'; 
bioslant_savename = 'example-LEM';

% ---------------- BioSLANT parameters -----------------------------------%

dt = 2e5; %time interval between landscapes (yrs)
nbasins = 20; %number of basins to populate
dw = 1; %upstream v. downstream dispersal weighting 
Acf = 1e6; %minimum drainage area for habitat

% -------------- Landscape Attributes ------------------------------------%

% parallel option :
% npool = 2;
% BioSLANT_Initialize(tadpole_output,bioslant_savename,dt,nbasins,Acf,dw,'pool',npool);

% non parallel option:
% BioSLANT_Landscape_Initialize(tadpole_output,bioslant_savename,dt,nbasins,Acf,dw);

lem_dir = [pwd,'/',bioslant_savename]; 
load([lem_dir,'/',bioslant_savename,'_1_.mat'],'p','g')
% ----------- Neutral Model Parameters -----------------------------------%

niter = 100;                    % number of stochastic realizations of the model to run

% DISPERSAL

p.dispersalfunction = '2dt';    % 'uniform' = equal probability of dispersal to any habitat within the basin
                                % 'none' = no dispersal, every offspring
                                % stays where they were 'born'
                                % 'neighbor' = can only disperse to a
                                %  neighboring habitat node

p.u = 11;                       % scale of dispersal (2dt)                                    
p.p = 0.9;                      % shape of distribution (2dt)

% HABITAT DISTRIBUTION 

p.hab_dist_k = 4;               % minimum number of fish in a habitat node 
p.hab_dist = 0.5;               % exponent for relating drainage area to habitat capacity
p.Ch = 5;                       % average number of fish in a habitat node

% SPECIATION 

p.theta = 8;                    % fundamental biodiversity number
p.tau = 10;                     % protracted speciation (in generations)
p.vicariance = 1;               % vicariance on (1) off (0)
p.track_ancestry = 0;           % track ancestry (1) or don't track ancestry(0)
p.doSaveAncestry = 0;           % save ancestry yes (1); or not (0)


% COUPLING BIO to TADPOLE

p.ngen = 10;                    % number of generations to run the NCM after changing landscape    
p.run_ic = 1;                   % run initial condition (1) or don't (0)
p.ic_con = 0;                   % initial condition 
                                % (1) = 1 species to start 
                                % (0) = infinite species to start
                               
p.ss_gen = 5e3;                  % number of iterations to reach steady state for the initial condition
p.tr_run = 0;                   % (1) run with changing landscape
                                % (0) run only on one static landscape

% SAVING

p.saveInt = 10;                  % how frequently to save (generations); to save only once
                                 % set this equal to ss_gen

% ----------- RUN MODEL ----------------------------------- %

rng('shuffle');
p.ss_gen = 4; p.saveInt = 2;
savename = 'example_BioSLANT_IC' ;  %name for saving the initial condition
% [p,g,speciesmat,tmat] = BioSLANTinitialize(p,g,niter,savename); %initial condition

load('example_BioSLANT_IC_2_IC.mat'); 

p.tr_run = 1;
savename = 'example_bioSLANT_TR'; 
BioSLANTrun(p,g,niter,lem_dir,bioslant_savename,speciesmat,tmat,savename);


