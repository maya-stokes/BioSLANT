function [rec,length] = receiver(h,p,g)

% FIND RECEIVERS -- from FastScape workshop -- Nov. 4, 2015

nx = p.Nx;
ny = p.Ny;
dx = p.dx;
rec = g.habitat; %initialize
length = dx*ones(p.nhabitats,1); %length between node and it's receiver (currently doesn't allow different dx and dy)
[J,I] = meshgrid(1:nx,1:ny); % matrix subscripts

Iinterior = I(2:ny-1,2:nx-1); %don't consider boundaries...for now
Jinterior = J(2:ny-1,2:nx-1);
ainterior = h(2:ny-1,2:nx-1);

alogical = ainterior > p.Acf; %logical where a exceeds minimum drainage area a fish needs to live
ind2use = find(alogical); %find the index of habitable nodes
i2use = Iinterior(ind2use); %interior i (original matrix)
j2use = Jinterior(ind2use); %interior j (original matrix)
ind = sub2ind([ny,nx],i2use,j2use);
%n = length(ind2use); %number of nodes to consider (find the receiver's of)

% iterate through 3x3 matrices and find steepest slope between each node
% and its neighbors, this node is the "receiver"
dxdiag = dx*sqrt(2);

% this script only set up from fixed boundary conditions
for k = 1:size(ind2use,1)
    
    i = i2use(k);
    j = j2use(k);
    listind = find(g.habitat == ind(k));
    
    x = h(i,j);
    
    ii = i + (-1:1);
    jj = j + (-1:1);
    xn = h(ii,jj); % surrounding points
    
    s(1,1) = (x - xn(1))/dxdiag; %diagonal slope
    s(2,1) = (x - xn(2))/dx; %perpendicular slope (take diagonal length into account!)
    s(3,1) = (x- xn(3))/dxdiag;
    s(1,2) = (x - xn(4))/dx;
    s(3,2) = (x-xn(6))/dx;
    s(1,3) = (x-xn(7))/dxdiag;
    s(2,3) = (x- xn(8))/dx;
    s(3,3) = (x-xn(9))/dxdiag;
    
    [jin, iin] = meshgrid(jj, ii); %index surrounding points in original matrix
    [val, imin] = min(s(:)); %find minimum
    
    if val <= 0 %if everybody is uphill, node is own receiver
        rec(listind) = sub2ind([ny,nx],iin(imin),jin(imin)); %otherwise it is neighbor with max slope
        if mod(imin,2) == 1 %if receiver is on a diagonal multiply length by sqrt(2)
            length(listind) = dxdiag;
        end
    end
end

% Periodic boundaries, if locations on edge have a receiver on the other
% side
% 0 = fixed, 1 = mirror, 2 = periodic
% left, right, upper, lower

indall = h > p.Acf;
[J,I] = meshgrid(1:nx,1:ny); % matrix subscripts

if p.bvec(1) == 2
    j2use = J(indall(:,1),1);
    i2use = I(indall(:,1),1);
    indleft = sub2ind([ny,nx],i2use,j2use);
    for jk = 1:size(indleft,1)
        i = i2use(jk);
        j = j2use(jk);
        if i == 1 || i == ny
            continue
        end %corners are addressed below
        x = h(i,j);
       listind = find(g.habitat == indleft(jk));

        xn(:,1) = h((i + (-1:1)),nx);
        xn(:,2) = h((i + (-1:1)),1);
        xn(:,3) = h((i + (-1:1)),2);
        
        ii = i + (-1:1);
        jj = j + (-1:1);
        jj(:,1) = nx;
        
        s(1,1) = (x - xn(1))/dxdiag; %diagonal slope
        s(2,1) = (x - xn(2))/dx; %perpendicular slope (take diagonal length into account!)
        s(3,1) = (x- xn(3))/dxdiag;
        s(1,2) = (x - xn(4))/dx;
        s(3,2) = (x-xn(6))/dx;
        s(1,3) = (x-xn(7))/dxdiag;
        s(2,3) = (x- xn(8))/dx;
        s(3,3) = (x-xn(9))/dxdiag;
        [jin, iin] = meshgrid(jj, ii); %index surrounding points in original matrix
        [val, imin] = min(s(:)); %find minimum
        
        if val <= 0 %if everybody is uphill, node is own receiver
            rec(listind) = sub2ind([ny,nx],iin(imin),jin(imin)); %otherwise it is neighbor with max slope
            if mod(imin,2) == 1 %if receiver is on a diagonal multiply length by sqrt(2)
                length(listind) = dxdiag;
            end
        end
    end
    
    j2use = J(indall(:,nx),nx);
    i2use = I(indall(:,nx),nx);
    indright = sub2ind([ny,nx],i2use,j2use);
    for jk = 1:size(indright,1)
        listind = find(g.habitat == indright(jk));
        i = i2use(jk);
        j = j2use(jk);
        x = h(i,j);
        
        if i == 1 || i == ny
            continue
        end %corners are addressed below
        
        xn(:,1) = h((i + (-1:1)),nx-1);
        xn(:,2) = h((i + (-1:1)),nx);
        xn(:,3) = h((i + (-1:1)),1);
        
        ii = i + (-1:1);
        jj = j + (-1:1);
        jj(:,3) = 1;
        
        s(1,1) = (x - xn(1))/dxdiag; %diagonal slope
        s(2,1) = (x - xn(2))/dx; %perpendicular slope (take diagonal length into account!)
        s(3,1) = (x- xn(3))/dxdiag;
        s(1,2) = (x - xn(4))/dx;
        s(3,2) = (x-xn(6))/dx;
        s(1,3) = (x-xn(7))/dxdiag;
        s(2,3) = (x- xn(8))/dx;
        s(3,3) = (x-xn(9))/dxdiag;
        [jin, iin] = meshgrid(jj, ii); %index surrounding points in original matrix
        [val, imin] = min(s(:)); %find minimum
        
        if val <= 0 %if everybody is uphill, node is own receiver
            rec(listind) = sub2ind([ny,nx],iin(imin),jin(imin)); %otherwise it is neighbor with max slope
            if mod(imin,2) == 1 %if receiver is on a diagonal multiply length by sqrt(2)
                length(listind) = dxdiag;
            end
        end
    end
    
    %top left corner --- need to fix length's
    if h(1,1) > p.Acf
        tr = indall(1,1);
        trind = find(g.habitat == tr);
        nlist = [h(2,1) h(1,2) h(2,2)];
        indlist = [2 ny+1 ny+2];
        [~,indmax] = max(nlist);
        rec(trind) = indlist(indmax);
    end
    
    %bottom left
    if h(ny,1) > p.Acf
        br = indall(ny,1);
        brind = find(g.habitat == br);
        nlist = [h(ny-1,1) h(ny,2) h(ny-1,2)  h(ny,nx) h(ny-1,nx)];
        indlist = [(ny-1) (2*ny-1) (2*ny) ny*nx ny*nx -1];
        [~,indmax] = max(nlist);
        rec(brind) = indlist(indmax);
    end
    
    %top right
    if h(1,nx) > p.Acf
        tl = indall(1,nx);
        tlind = find(g.habitat == tl);
        nlist = [h(1,nx-1) h(2,nx-1) h(2,nx)];
        indlist = [(ny*(nx-2) +1) (ny*(nx-2) + 2) (ny*(nx-1) + 2)];
        [~,indmax] = max(nlist);
        rec(tlind) = indlist(indmax);
    end
    
    %bottom right
    if h(ny,nx) > p.Acf
        bl = indall(ny,nx);
        blind = find(g.habitat == bl);
        nlist = [h(ny-1,nx) h(ny-1,nx-1) h(ny,nx-1)];
        indlist = [(ny*nx-1) (ny*(nx-1)) (ny*(nx-1)-1)];
        [~,indmax] = max(nlist);
        rec(blind) = indlist(indmax);
    end
    
end



