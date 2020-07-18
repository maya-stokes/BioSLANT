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
    dA_all = abs(A - Aprev);
    dA = []; 
    for j = 1:4 %only use non-periodic boundaries
       if p.bvec(j) == 0 || p.bvec(j) == 1
           if j == 1
               vec = dA_all(2:end-1,1)'; 
           elseif j == 2
               vec = dA_all(2:end-1,end)'; 
           elseif j == 3
               vec = dA_all(1,:); 
           elseif j == 4
               vec = dA_all(end,:); 
           end
           dA = [dA vec]; 
       end
    end
%     dA = [dA(1,:) dA(end,:) dA(2:end-1,1)' dA(2:end-1,end)']; % just boundary points
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
    
    %find which ones are large enough; 
    capind = sub2ind(size(dB),icap,jcap); 
    aind = A(capind);  %indices that exceed drainage area and switched basins
    capind = capind(aind>Amin); aind = aind(aind>Amin); 
    cnum = 0; 
    
    if any(capind) %if captures
        [aind,ix] = sort(aind,'descend'); capind = capind(ix); %sort by drainage area
        [D,~,~] = mexD8Dir(g.U,p.dy/p.dx,1*(~g.C),p.bvec,p.flood);
        
        Basin_all = false(size(dB)); 
        for i = 1:length(aind) %go through each and if it has already been accounted for skip, save all pixels within the basin

            if Basin_all(capind(i)) %if in another basins don't count 
                continue
            else
                
                [iout,jout] = ind2sub(size(dB),capind(i)); %if not in another basin count, 
                [Basins,~] = mexBasinModel(D,iout,jout,p.bvec); %find basin 
                Basin_all(Basins == 1) = 1; %track what has already been labeled basin
                
                [theicap, thejcap] = find(Basins); 
                
                cnum = cnum+1;
                Stats.Acap{framenum}(cnum) = aind(i);
                Stats.Acaptot(framenum) = Stats.Acaptot(framenum) + aind(i);

                Stats.icap{framenum}{cnum} = theicap;
                Stats.jcap{framenum}{cnum} = thejcap;
 
            end
        end    
    end
    
    tprev = t(n);
    Aprev = A;
    Bprev = B;
end

save([basename '_Stats.mat'],'Stats','Amin')

