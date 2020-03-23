function DrawPlot(n,p,g)

% to plot fault trace select 'elevation' as the choice for plotting, 
% and set p.plotFaults = 1

% to make a movie with captures, save results and run the script
% "Capture_Movie.m", currently set-up only to run this with elevation as
% the background...

switch p.plottype

    case 'drainage area' % map view with log(A) as color
        
        figure(p.fighandle)        
        imagesc(g.x,g.y,log10(eval(['mex' p.routing '(g.U,p.dy/p.dx,1*(~g.C),p.bvec,p.flood)'])))
        title(['Iteration=' num2str(n) ', t=' num2str(round(p.t)) ', dt=' num2str(round(p.dt))])
        drawnow                

    case 'curvature' % map view with normalized curvature as color
        
        figure(p.fighandle)        
        imagesc(g.x,g.y,curvn(g.U))
        title(['Iteration=' num2str(n) ', t=' num2str(round(p.t)) ', dt=' num2str(round(p.dt))])
        drawnow                

    case 'elevation' % map view with elevation as color
        
        figure(p.fighandle)
        cla
        
        
        imagesc(g.x,g.y,g.U)
        axis image
%         set(gca,'ydir','normal')
        title(['Iteration=' num2str(n) ', t=' num2str(round(p.t)) ', dt=' num2str(round(p.dt))])
        colorbar
        
         % faults
         if p.doFaults && p.plotFaults
             [~,g] = FaultDist(p,g);% Calculate fault plane locations and distances
             nf = size(g.FaultDist,3);
             for i = 1:nf
                 contour(g.x,g.y,g.FaultDist(:,:,i),[0 0],'w','LineWidth',3)
             end
         end
         
         if isfield(g,'tc') ==1 %if you want to draw captures
             % captures
             a0 = 1.0; % face alpha value at zero capture age (1 = opaque)
             cage = p.t - g.tc; % Capture age. The time now is p.t, the time of each capture is in g.tc
             xcb = g.xcb(cage>=0); % get captures that have already happened and their ages
             ycb = g.ycb(cage>=0);
             cage = cage(cage>=0);
             
             % draw a patch for each capture, with opacity proportional to age
             for c = 1:length(cage)
                 patch(xcb{c},ycb{c},'w','EdgeColor','none','FaceAlpha',a0*exp(-cage(c)/p.capdecay));
             end
         end
         
         %         set(gca,'visible','off')
         %         colorbar off
         
         drawnow
         
    case 'contour' % normalized contour map
        
        figure(p.fighandle)        
        ncont=10;
        contour(g.X,g.Y,(g.U-min(g.U(:)))/range(g.U(:)),ncont)
        title(['Iteration=' num2str(n) ', t=' num2str(round(p.t)) ', dt=' num2str(round(p.dt))])
        drawnow        

    case 'shade' % shaded relief
        
        figure(p.fighandle)
        Z.grid = flipud(g.U);
        Z.x = g.x;
        Z.y = g.x;
        Z.dx = p.dx;
        Z.dy = p.dy;
        az = 315;
        alt = 45;
        Shade(Z,az,alt,p.zfactor);
        title(['Iteration=' num2str(n) ', t=' num2str(round(p.t)) ', dt=' num2str(round(p.dt))])
        drawnow

    case 'color shade' % colored shaded relief
        
        figure(p.fighandle)
        Z.grid = flipud(g.U);
        Z.x = g.x;
        Z.y = g.x;
        Z.dx = p.dx;
        Z.dy = p.dy;
        az = 315;
        alt = 45;
%         ColorShade(Z,Z,az,alt,p.zfactor);
        [~,h] = ColorShade(Z,Z,az,alt,p.zfactor);
        A = zeros(size(p.F));
        A(p.F==1) = 0;
        A(p.F==0) = 1;
        set(h,'alphadata',flipud(A))
%         colorbar off
%         set(gca,'visible','off')
        title(['Iteration=' num2str(n) ', t=' num2str(round(p.t)) ', dt=' num2str(round(p.dt))])
        drawnow

    case 'mesh' % mesh plot

        figure(p.fighandle)
        surf(g.X,g.Y,flipud(g.U)); 
        vertexag = [1 max(g.Y(:))/max(g.X(:)) p.zfactor*range(g.U(:))/max(g.X(:))]; 
        set(gca, 'PlotBoxAspectRatio', vertexag);
        title(['Iteration=' num2str(n) ', t=' num2str(round(p.t)) ', dt=' num2str(round(p.dt))])
        drawnow
 

    otherwise % case 1 and default: perspective surface plot

        figure(p.fighandle)
        surf(g.X,g.Y,flipud(g.U)); 
        shading interp
        camlight left 
        lighting phong
        vertexag = [1 max(g.Y(:))/max(g.X(:)) p.zfactor*range(g.U(:))/max(g.X(:))]; 
        set(gca, 'PlotBoxAspectRatio', vertexag);
        title(['Iteration=' num2str(n) ', t=' num2str(round(p.t)) ', dt=' num2str(round(p.dt))])
        drawnow

end