function BioSLANT_Landscape_Initialize(tadpole_output,bioslant_savename,dt,nbasins,Acf,dw,varargin)
%BioSLANT_Initialize landscape files needed to run the neutral model
%
% Syntax
%
%     BioSLANT_Landscape_Initialize(tadpole_output,bioslant_savename,dt,nbasins,varargin)
%
%
% Description
%
%     BioSLANT_Initialize take the output files from Tadpole, saved in a
%     .mat file and calculates the habitat capacity and dispersal matrices
%     needed to run the neutral community model.
%
%     BioSLANT_Initialize(tadpole_output,bioslant_savename,dt,nbasins,'pool',np)
%     will run a parfor loop when calculating the dispersal matrices
%
% Input arguments
%
%     tadpole_output    string specifying the name of the .mat file where the
%                       Tadpole output is saved
%     bioslant_savename string specifying the name of a directory/root name 
%                               of the files where the dispersal matrices/habitat
%                               capacity information will be saved
%     dt                how frequently (in time) to sample the outputs from the landscape
%                       evolution model, in Stokes and Perron 2020, we used 2e5.
%     nbasins           specify number of basins to populate with
%                       organisms. To populate all available habitat with organisms
%                       set this equal to 0.
%     Acf               critical drainage area to define habitat (m^2)
%     dw                downstream dispersal weighting. 1 = no weighting. 
%
% Output arguments
%
%     creates a directory with ".mat" files containing "p" (the parameter
%     structure) and "g" (the data structure) 
%
% Author: Maya Stokes (mstokes@mit.edu) and Taylor Perron (perron@mit.edu)
% Date: 19 March 2020

% --------------- START PARPOOL ----------------------------------------- %
if nargin > 6 %delete and create a parallel pool
    delete(gcp('nocreate'))
    parpool(varargin{8});
end

% ------------ LOAD TADPOLE OUTPUT -------------------------------------- %
load(tadpole_output,'output','t','p'); 
int = round(dt./mean(diff(t))); %select Tadpole landscape timesteps at an equal interval in time
int1 = 1; %where to start the NCM (depending on the perturbation you put into Tadpole you may get a lot of captures right off the bat). 
ncm.Z = output(:,:,int1:int:end); 
nlandscapes = size(ncm.Z,3); 
t = t(int1:int:end); 

% ----------- HABITAT AND DISPERSAL PARAMETERS ----------------------------%
p.Acf = Acf; 
p.dw = dw; 

% ----------- SET-UP THE LANDSCAPES ------------------------------------- %
p.routing = 'D8'; %set some of the inputs necessary for drainage area calculation 
g.C = ones(p.Ny, p.Nx); g.C(1,:) =0; g.C(:,end) = 0; g.C(:,1) = 0; g.C(end,:) = 0; %for 4 fixed boundaries

for i = 1:nlandscapes 
    g.t = t(i); %time
    g.U = ncm.Z(:,:,i); %elevation
    g = DrainageArea(p,g); %drainage area
    ncm.A(:,:,i) = g.A;
end
ncm.C = g.C; %boundary conditions

% ------------ DEFINE BASINS TO POPULATE -------------------------------- %

alist = [ncm.A(1,:,1)'; ncm.A(p.Ny,:,1)';ncm.A(:,1,1); ...
    ncm.A(:,p.Nx,1)];
ind_edge = [ones(p.Nx,1) (1:p.Nx)'; p.Ny*ones(p.Nx,1) (1:p.Nx)'; (1:p.Ny)' ones(p.Ny,1) ; (1:p.Ny)' p.Nx*ones(p.Ny,1)];

if nbasins ~=0 %choose nbasins
    [~,i_sort] = sort(alist,'descend');
    basinij = ind_edge(i_sort(1:nbasins),:); %choose nbasin largest basins 
    p.bp = basinij; %basin pour points  
else
    p.bp = ind_edge(alist > Acf,:); %all edge nodes that have at least critical drainage area for habitat
end


setuplandscape(p,bioslant_savename,nlandscapes,ncm); %call function to set-up landscape
