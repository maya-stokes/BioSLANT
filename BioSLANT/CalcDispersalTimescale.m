function td = CalcDispersalTimescale(rho,theta,p_p,p_u,niter,ngen,nx)
%calculate the generations it requires to traverse half of the domain on a
%straight 200x1 landscape 

 % Inputs: 
    %nfish = number of fish units in the model
    %theta = fundamental biodiversity number (mu*nfish)
    %tau = timescale of protracted speciation (generations)
    %p, u = parameters in the dispersal kernel 
    %niter = number of iterations to run the dispersal simulation for 
    %ngen = number of generations to let the model run for
    %nx = length of the domain 
    
 % Outputs: 
    %ts = speciation timescale (generations)
    
rng('shuffle'); 
nfish = rho*nx;
mu = theta/nfish; %calculate probability of speciation

p.p = p_p; 
p.u = p_u; 

% build a non-spatially explicit simulation where at each timestep we
% choose one fish to die. Then, we sample the dispersal kernel and either
% allow a fish to disperse or speciate to populate the empty spot. If the
% fish speciate, then we set the total distance to zero. If they disperse,
% we add the total distance the disperse has already moved during the
% simulation and the distance they are moving during this timestep. It's
% not exactly that an individual is moving every timestep, more like they
% are moving to a spot, laying an egg, then going back home and we
% accumulate all of the movements of that lineage together. 

meand = zeros(niter,1); 
distance_matrix = zeros(nx+1,nx); 
kernel_matrix = zeros(nx+1,nx); 
index_matrix = zeros(nx+1,nx); 
distance = 1:nx; 

for i = 1:nx
dist = abs(i - distance); %right now I am treating a dispersal event from the location of the death = distance of 0, this could be edited.
K = p.p./(pi*p.u^2*(1+(dist./p.u).^2).^(p.p+1)); %calculate dispersal function
P = [(1-mu).*(K.*nfish/sum(K.*nfish)) mu]; %probability
dist_r = [dist NaN]; %add in the probability of speciation 
[P,is] = sort(P); %sort and save the index
kernel_matrix(:,i) = P; %save matrices to reference in loop below
distance_matrix(:,i) = dist_r(is); 
index_matrix(:,i) = is; 
end

[sp_ind,~] = find(isnan(distance_matrix)); %this is the index where speciation would occur

for jj = 1:niter
    fish_distance = zeros(nfish,1); %initalize the distance traveled by a fish and its ancestors vector
    
    %assign a location to each of the fish, distributed uniformly
    %throughout the number of nodes; divide evenly at first and then randomly distribute the remainder  
    fl = floor(nfish./nx); rem = nfish - nx*fl;     
    hab = repmat(1:nx,fl,1);     
    fish_hab = [hab(:); randi(nx,rem,1)]; 

    
   for j = 1:nfish*ngen
        fish_ind = randi(nfish,1,1); %first choose a random fish
        fish_ind_hab = fish_hab(fish_ind); 
        fish_distance(fish_ind) = 0; %set her distance to zero      
      
        P = kernel_matrix(:,fish_ind_hab); 
        val = min(P) + (max(P) - min(P)).*rand(1,1); %select random number in the range of the probability of dispersal
        [~,~,bin] = histcounts(val,[0 ; P]); %choose a random number to figure out which habitat bin to choose from
        
        if bin ~= sp_ind(fish_ind_hab) %if speciation didn't happen....do dispersal
            mig_ix = index_matrix(bin,fish_ind_hab); %find the position of the habitat you're looking at (that is the correct distance away from the fish death)
            fish_bin = find(fish_hab == mig_ix); %find the fish that live in that position
            fish_ix = randi(length(fish_bin),1,1); %randomly choose one to be the replacement             
            fish_distance(fish_ind) = fish_distance(fish_bin(fish_ix)) + distance_matrix(bin,fish_ind_hab); %distance already recorded by the fish and its ancestors + this distance it travels this timestep
        end
        
   end
    
    meand(jj) = mean(fish_distance); %mean distance travled by a fish and its descendants
end

td = ngen*(nx/mean(meand))*(1/2); %number of generations required to go half the domain
end
