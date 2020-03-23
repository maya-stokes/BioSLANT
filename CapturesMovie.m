function CapturesMovie(output,p,t,frame_interval,basename,varargin)

% TadpoleMovie(output,p,t,frame_interval,basename)
%
% produces an image sequence from a saved Tadpole run using saved elevation
% grids in output, parameters struct p, time vector t, drawing an image every
% frame_interval saved grids with image file name root specified by basename.
% The image sequence can then be used to produce an animation with, e.g.,
% Quicktime's "Load image sequence" function. Using the default Quicktime
% codec at 15 fps produces a smooth animation suitable for Powerpoint or
% Keynote.

% Stats is a structure returned by TransientStats.m
% capdecay is the e-folding time for capture patch opacity as a fraction of
% the model run duration

if nargin > 5
    Stats = varargin{1}; 
    capdecay = varargin{2}; 
end
    
v = VideoWriter([basename,'.avi'],'Motion JPEG AVI'); v.FrameRate = 20; 
open(v)

[p.K,p.J,p.N] = size(output);

p.zmin = min(output(:));
p.zmax = max(output(:));


if ~isfield(p,'plottype')
    p.plottype = 'elevation';
end

if isempty(t)
    t = (0:p.N-1)*p.saveint;
end


% loop through the frames
framenum = 0;
% tprev = 0; % model time of previous frame

if p.plotCaptures == 1
    g = getCaptures(p,Stats); % g now contains g.tc, model times when captures occurred; g.xcb and g.ycb, cell arrays of x and y coordinates of capture-bounding polygons
    p.capdecay = capdecay*p.tf; % now capdecay is in years
end

for n=1:frame_interval:p.N

    framenum = framenum+1;
    
    g.U = output(:,:,n);
    p.t = t(n);
    
    if n==1
        [p,g] = SetUpPlot(p,g);
    else
        DrawPlot((n-1)*p.saveint,p,g); 
    end
   
   
   frame = getframe(gcf); 
   writeVideo(v,frame); 
end
    
% close the figure

if isfield(p,'fighandle') == 1
    close(p.fighandle);
end

close(v)

function g = getCaptures(p,Stats)

g.x = p.dx*(0:p.J-1);
g.y = p.dy*(0:p.K-1);

shrinkfactor = 1;

nf = length(Stats.t); % number of stats frames

nc = 0; % total number of captures

for f = 1:nf % loop through frames

    for c = 1:length(Stats.Acap{f}) % number of captures in this frame; loop through captures
        nc = nc + 1; % number of captures including this one
        g.tc(nc) = Stats.t(f); % model time for this model frame 
        ic = Stats.icap{f}{c};
        jc = Stats.jcap{f}{c};
        xc = g.x(jc)';
        yc = g.y(ic)';
        b = boundary(xc,yc,shrinkfactor);
        g.xcb{nc} = xc(b);
        g.ycb{nc} = yc(b);

    end
end

