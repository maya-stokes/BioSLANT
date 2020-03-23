function Stats = TransientStats(output,p,t,frame_interval,Amin,basename)

% TransientStats(output,p,t,frame_interval,basename)

% analyze changing drainage areas in fault runs
% Amin = number of pixels for drainage area to count as a big channel 

[~,~,N] = size(output);

if isempty(t)
    t = (0:N-1)*p.saveint;
end

g.U = output(:,:,1);
[~,g] = BoundaryMat(p,g);

numframes = length((frame_interval+1):frame_interval:N);
Stats.t = zeros(numframes,1);
Stats.dAdt = zeros(numframes,1); % d(drainage area in grid cells)/dt on boundaries from timestep n-1 to timestep n
Stats.dA = zeros(numframes,1); % total area exchanged (in grid cells) on boundaries from timestep n-1 to timestep n
Stats.Acaptot = zeros(numframes,1); % total area captured (in grid cells) from timestep n-1 to timestep n
Stats.Acap = cell(numframes,1); % each element is a vector containing areas of captures (in grid cells) from timestep n-1 to timestep n
Stats.icap = cell(numframes,1); % each element is itself a cell array of vectors containing row indices of captures from timestep n-1 to timestep n
Stats.jcap = cell(numframes,1); % each element is itself a cell array of vectors containing column indices of captures from timestep n-1 to timestep n
% Stats.vcap = cell(numframes,1); % each element is a vector containing unique values of captures from timestep n-1 to timestep n


% loop through the frames
framenum = 0;
tprev = 0;
% Zprev = output(:,:,1); % elevations at previous frame
g = DrainageArea(p,g);
Aprev = g.A/(p.dx*p.dy); % drainage areas at previous frame, in units of dx*dy

% watersheds for boundary points at previous frame
Bprev = getBasins(p,g);

for n=(frame_interval+1):frame_interval:N

    framenum = framenum+1;
    Stats.t(framenum) = t(n);
    g.U = output(:,:,n); % elevations at frame n
    g = DrainageArea(p,g);
    A = g.A/(p.dx*p.dy); % drainage areas at frame n
    dA = abs(A - Aprev);
    dA = [dA(1,:) dA(end,:) dA(2:end-1,1)' dA(2:end-1,end)']; % just boundary points
    Stats.dA(framenum) = 0.5*sum(dA); % total area exchanged between boundary points. We multiply by 1/2 because otherwise we could count each divide point displacement twice
    Stats.dAdt(framenum) = mean(dA)/(t(n) - tprev); % average dA/dt

    B = getBasins(p,g); % Basins that drain to boundaries, labeled with integers starting with 1 (zero for points that don't drain to boundaries)
        % Note: The max possible integer in B is 2*Nx + 2*(Ny-2) 
        
    % Assign unique integers to pairs of basin indices using the Cantor
    % pairing function 
    dB = 0.5*(Bprev + B).*(Bprev + B + 1) + B; % the pairing function
    dB(B == Bprev) = 0; % set locations where basin didn't change to zero
        
    % now record sizes of captures
    [icap,jcap,vcap] = find(dB); % indices and unique values of captures
    dBvals = unique(vcap); % dBvals is a list of unique capture values
                                    
    % calculate the area of each capture in grid cells [A1 A2 A3 etc]
    % should be able to do this with e.g. sum(vcap==dBvals(i)) for i=1:length(dBvals)
    cnum = 0; % counts captures larger than Amin for this framenum
    
    for i = 1:length(dBvals)
        thecap = (vcap==dBvals(i));
        
        theicap = icap(thecap);
        thejcap = jcap(thecap);
        capind = sub2ind(size(dB),theicap,thejcap); % linear indices of the capture in dB
        [~,maxAind] = max(A(capind));
        maxAi = theicap(maxAind); % the row and column indices of the maximum drainage area in the capture
        maxAj = thejcap(maxAind);
        
        CapBin = zeros(size(dB));
        CapBin(capind) = 1; % binary image that is 1 for the capture, zero elsewhere
        CapBin = bwselect(CapBin,maxAj,maxAi,8); % now CapBin is 1 only for points that are 8-connected to the maximum drainage area point
        
        [theicap,thejcap] = find(CapBin); % row and column indices of the contiguous capture
        
        acap = length(theicap); % number of grid points in the contiguous capture
        if acap >= Amin % don't record the ones that are too small
            cnum = cnum+1;
            Stats.Acap{framenum}(cnum) = acap;
            Stats.Acaptot(framenum) = Stats.Acaptot(framenum) + acap;
            Stats.icap{framenum}{cnum} = theicap;
            Stats.jcap{framenum}{cnum} = thejcap;
        end
    end
            
    tprev = t(n);
    Aprev = A;
    Bprev = B;
end

save([basename '_Stats.mat'],'Stats','Amin')

