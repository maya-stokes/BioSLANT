%Example: from Tadpole output to fully coupled run.

% --------------- Filenames ----------------------------------------------%

tadpole_output = 'example_fault1.mat'; 
bioslant_savename = 'example-LEM';
% ---------------- BioSLANT parameters -----------------------------------%

dt = 2e6; %time interval between landscapes (yrs); Stokes Perron (2020) use 2e5
nbasins = 20; %number of basins to populate
dw = 1; %upstream v. downstream dispersal weighting 
Acf = 1e6; %minimum drainage area for habitat

% -------------- Landscape Attributes ------------------------------------%

% parallel option :
npool = 2;
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
p.doSaveAncestry = 0;           % save ancestry yes (1); or not (0); 
                                % all model runs in Stokes and Perron
                                % (2020) did not save ancestry


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
%NOTE; parameters selected below (p.ss_gen, p.ngen, niter, and saveInt) are
%selected to ensure the model can run on a personal computer. To run for
%5000 generations/100 iterations at least 20 GB of RAM are recommended. 

% to experiment with the code on a dynamic landscape, use the example
% results in the StokesPerron2020_ResultExample folder: 
% load('IC_U_LEM1a_TAU11_P2.5_U2IC_r197.mat','speciesmat'); 
% ^ use the initial condition results in the StokesPerron2020_ResultExample\Example Steady State Output\'IC_U_LEM1a_TAU11_P2.5_U2IC_r197.mat
% for a ready-made NCM result at steady state. The landscape files
% associated with this example can be found in StokesPerron2020_ResultExample\Landscape Scenario 1 NCM\LEM1-a\LEM1-a_1_.mat'


% RUN THE CODE HERE 
rng('shuffle');
p.ss_gen = 4; p.saveInt = 2; niter = 1; %NOTE: short run-time specified here; long-run times can eat up quite a bit of RAM; we recommend trying these out on your computer first
savename = 'example_BioSLANT_IC' ;  %name for saving the initial condition
[p,g,speciesmat,tmat] = BioSLANTinitialize(p,g,niter,savename); %initial condition

load('example_bioSLANT_IC_2_IC.mat')

p.tr_run = 1;
p.saveInt = 1; p.ngen = 2;niter = 1; %specificed so that the code can run on a laptop; p.ngen = 10 and p.saveInt = 10 in Stokes and Perron (2020). 
savename = 'example_bioSLANT_TR'; 
BioSLANTrun(p,g,niter,lem_dir,bioslant_savename,speciesmat,tmat,savename);


