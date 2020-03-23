function [basins] = extractbasins(p,g)
% Makes a matrix labelled by drainage basin 

%flatten sparse matrix g.W of D8 drainage directions
g = DrainageArea(p,g); 
D =g.W;
D(D==.5)=1;                       %change all direction originating from mirror boundaries to 1
drs = cumsum(ones(p.Ny,p.Nx,8),3);  %make a k by j matrix where each layer is 1 through 8
D=sum(D.*drs,3);                  %matrix of drainage directions

%find basins for each outlet
noutlet = length(p.bp);
basins = zeros(p.Ny,p.Nx);
for i = 1:noutlet
    [Basin, ~] = mexBasinModel(D,p.bp(i,1),p.bp(i,2),p.bvec); %calls Tadpole's BasinModel code
    basins(Basin==1) = i;
end
