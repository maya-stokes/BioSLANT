function [p,g] = FaultDist(p,g)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute distance of each point on land surface from a set of faults %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Give the coefficients defining the equation of each planar fault at t=0
% in terms of strike, dip, and elevations at given x, y, z.

% Example values where we see late stage drainage reorganization
% strike(1) = 60; dip(1) = 70; x0(1) = p.dx*0.5*(p.Nx-1); y0(1) = p.dy*0.5*(p.Ny-1); z0(1) = 0; % strike in degrees E of N, dip in degrees down from horizontal (right-hand rule), elevation at x=0,y=0 (meters)
% strike(1) = 60; dip(1) = 60; x0(1) = p.dx*0.3*(p.Nx-1); y0(1) = p.dy*0.7*(p.Ny-1); z0(1) = 0; % strike in degrees E of N, dip in degrees down from horizontal (right-hand rule), elevation at x=0,y=0 (meters)

% Update the constant term to account for cumulative uplift since t=0
z0 = p.z0 + p.t*p.E;

% The fault normal vector components are:
% x: sin(dip)*cos(strike)
% y: -sin(dip)*sin(strike)
% z: cos(dip)

strike = deg2rad(p.strike);
dip = deg2rad(p.dip);
nvec = [sin(dip(:)).*cos(strike(:)), -sin(dip(:)).*sin(strike(:)), cos(dip(:))]; % each row of nvec is the normal vector (nx,ny,nz) of a plane

% Points on the planes (x0,y0,z0)
P = [p.x0(:), p.y0(:), z0(:)]; % each row of P is a point on a plane

% The equation of each plane is
% a(x-x0) + b(y-y0) + c(z-z0) = 0
% or
% ax + by + cz + d = 0, with d = -(a*x0 + b*y0 + c*z0)
planes = [nvec, -sum(nvec.*P,2)]; % each row is a, b, c, d for a plane

% Establish x,y,z coordinates (in m) of each point on land surface
g.x = p.dx*(0:p.J-1);
g.y = p.dy*(0:p.K-1);
[X,Y] = meshgrid(g.x,g.y); 

nf = length(strike); % number of faults

sz = [1,1,nf];
a = reshape(planes(:,1),sz);
b = reshape(planes(:,2),sz);
c = reshape(planes(:,3),sz);
d = reshape(planes(:,4),sz);

g.FaultDist = (a.*repmat(X,sz) + b.*repmat(Y,sz) + c.*repmat(g.U,sz) + d)./sqrt(a.*a + b.*b + c.*c); % signed distances to faults (see below)
g.Dmin = min(abs(g.FaultDist),[],3); % matrix of minimum absolute distance to any fault

% calculate a multiplier (between 0 and 1) to give effective K at each point
g.FaultFactor = 1 + (p.maxFaultFactor-1)*exp(-g.Dmin/p.Dstar);
end