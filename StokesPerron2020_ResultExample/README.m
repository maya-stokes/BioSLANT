% README

% Directory for StokesPerron2020_ResultExample: 

% Example Steady State Output (unzip to use): results sampled every 10
% generations for a NCM at steady state (TAU = 11 and U = 2).

% Example Transient Output (unzip to use): NCM run for 10 generations on
% each landscape passed to BioSLANT. Includes vicariant speciation
% following river capture. 

%Landscape Scenario 1 NCM (unzip to use): habitat, basins and dispersal
%distances used by BioSLANT for LEM Scenario 1 in Stokes and Perron. 

%Landscape Scenario 1 Tadpole.mat: Tadpole output used to create the files
%in "Landscape Scenario 1 NCM" 

%Landscape Scenario 1 Tadpole Capture Stats.mat: Capture/drainage basin
%area exchange statistics output.

%Landscape Scenario 1 Tadpole Visualization: Movie highlighting fault
% trace migration and river captures.

% TimeSeriesFigure.m: script to calculate species richness, speciation, extinction and
% diversification (Figure 5 in Stokes and Perron). 

% SpeciesArea.m: script to calculate the species-area relationship (Figure S6-7 in Stokes
% and Perron). 

% NOTE ON PARAMETERS: 
% parameter space explored in Figure 6, Stokes and Perron is listed below. 
% The remainder of the parameters are identical to those listed in the "p" structure in the
% output files in this folder. 

% Tadpole: p.dip = [30 45 60]; 
% NCM: p.u = [1 1.2 1.4 1.6 1.8 2 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39]; 
%      p.tau = [1:10:200]; 

