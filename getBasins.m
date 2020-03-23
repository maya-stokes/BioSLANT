function Basins = getBasins(p,g)

% get D8 drainage directions
[D,~,~] = mexD8Dir(g.U,p.dy/p.dx,1*(~g.C),p.bvec,p.flood);

% find indices of outlet points based on the matrix g.C
[iout,jout] = find(1*(~g.C));

% get Basins that drain to the specified outlets
[Basins,~] = mexBasinModel(D,iout,jout,p.bvec);