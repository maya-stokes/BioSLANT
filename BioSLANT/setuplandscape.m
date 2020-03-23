function  setuplandscape( p,fname,nlandscapes,ncm)
%Set-up and save landscape attributes in mat-files

A = ncm.A; %drainage area
Z = ncm.Z; %elevation
C = ncm.C; %Tadpole boundary conditions

mkdir(fname); dr = pwd; %add lem file to path
fullpath = [dr '/',fname]; addpath(fullpath); 

pl = gcp('nocreate'); % If no pool, do not create new one.

% ------ REGULAR FOR LOOP ----------------------------------------------- %
if isempty(pl)
    for i = 1:nlandscapes
        g = struct();
        g.A = A(:,:,i); %area
        g.U = Z(:,:,i); %elevation
        g.C = C; %Tadpole boundary conditions
        
        basins = extractbasins(p,g); %define basins
        
        g.habitat = find(g.A > p.Acf & basins ~= 0); %linear index of habitat
        g.basinind = basins(g.A > p.Acf & basins ~= 0); %list of basin ID's for each habitat
        p.nhabitats = length(g.habitat); %number of habitats
        g = dispersaldistancetable(p,g); %dispersal distance matrix
        
        save([fullpath,'/',fname,'_',num2str(i),'_'],'p','g')
    end
    
else %PARALLEL LOOP; same loop but set-up for parfor loop
    parfor i = 1:nlandscapes
        g = struct();
        pp = p;
        
        g.A = A(:,:,i); %area
        g.U = Z(:,:,i); %elevation
        g.C = C;
        
        basins = extractbasins( pp,g); %define basins
        
        g.habitat = find(g.A > pp.Acf & basins ~= 0);
        g.basinind = basins(g.A > pp.Acf & basins ~= 0);
        pp.nhabitats = length(g.habitat);
        g = dispersaldistancetable(pp,g); %dispersal distance
        
        parsave([fullpath,'/',fname,'_',num2str(i),'_'],pp,g)
    end 
end
end

