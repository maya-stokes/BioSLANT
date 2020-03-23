function [ p1,g1,t_inc1,inc_ind1] = LandscapeTransition( p1,g1,g0,ic,t_inc0_array)
% LandscapeTransition transitions to the next landscape timestep. 
% Description
%
%       1) Does vicariance in the case of river capture
%       2) Randomly selects fish to die in habitat nodes with drainage area
%           loss
%       3) Allows fish to disperse throughout the landscape to occupy newly
%           available habitat
%
% Input arguments
%
%     p1    parameter structure for current landscape timestep
%     g1    data structure for current landscape timestep
%     g0    data structure for previous landscape timestep 
%     ic    list of species identity from previous timestep
%     t_inc0_arry   vector used to carry-over information about protracted
%                   speciation
%


listind = zeros(p1.nhabitats,1); %vector that tells you where you're allowed to disperse from

[~,g1ind, g0ind] = intersect(g1.habitat,g0.habitat); %habitats in initial landscape that are also in new landscape
nfish0 = numel(ic); %number of fish in initial landscape
nhab = length(g1ind); %number of habitats shared
oldhc = g0.habitatcapacity(g0ind); newhc = g1.habitatcapacity(g1ind); %old and new habitat capacities
capdif = newhc - oldhc; %difference between the two
disp_ind = find(capdif > 0); %where you'll need to do dispersal
[~,disp_fish_ind] = ismember(g1ind(disp_ind),g1.fishhab);
ndis = length(disp_ind);

n_incipient = length(t_inc0_array); % number of indices you need to transfer over to t_inc1
row = zeros(n_incipient,1);
val = zeros(n_incipient,1); count = 1;

row_in = t_inc0_array(:,1); 
time_in = t_inc0_array(:,2);

% FIRST do vicariance
if p1.vicariance == 1
    capind = g1.basinind(g1ind) - g0.basinind(g0ind)~=0; %not in same basin
    if any(capind) 
        fish_ind = find(ismember(g0.fishhab,g0ind(capind))); %fish in captured reach
        inc_ind = ismember(row_in,fish_ind); % reset incipient speciation
        if any(inc_ind)
            time_in(inc_ind) = 0; 
        end
        
        fish_cap = ic(fish_ind); %species in captured reach
        specieslist = unique(fish_cap); %unique list of species in captured reach
        nsp = length(specieslist); %number of unique species in captured reach
        g1.anc_vic = specieslist; g1.child_vic = zeros(nsp,1); %initialize ancestry vectores
        
        not_cap = true(nfish0,1); not_cap(fish_ind) = 0; sp_notincap = ic(not_cap); %list of species not in the capture
        for kk = 1:nsp %for each species captured
            %if fish occurs outside of the captured reach (otherwise all
            %are transfered and there is no "splitting of a population)"
            if any(sp_notincap == specieslist(kk))
                ic(fish_ind(fish_cap == specieslist(kk))) = p1.nspecies+1; %speciate
                p1.nspecies = p1.nspecies+1; %track number of species that have existed
                g1.child_vic(kk) = p1.nspecies; %children
            end
        end
    end
end


for kk = 1:nhab
    i1 = find(g1.fishhab == g1ind(kk)); i0 = find(g0.fishhab == g0ind(kk)); %index of fish living in shared habitat
    if capdif(kk) == 0 %if no difference in habitat capacity, stays the same
        g1.speciesvec(i1) = ic(i0); %assign fish
        i = i1(ismember(i0,row_in));
        v = time_in(ismember(row_in,i0));
        if any(i) %if incipient
            n_add_ind = length(i);
            for jik = 1:n_add_ind
                row(count) = i(jik); val(count) = v(jik);
                count = count +1;
            end
        end
        
    elseif capdif(kk) < 0 % if habitat loss
        fdind = randperm(oldhc(kk),newhc(kk)); %randomly choose newhc integers from 1:oldhc
        ind_old2keep = i0(fdind); %the lucky ones
        g1.speciesvec(i1) = ic(ind_old2keep); %fill new speciesvec
        i = i1(ismember(ind_old2keep,row_in));
        v = time_in(ismember(row_in,ind_old2keep));
        if any(i)
            n_add_ind = length(i);
            for jik = 1:n_add_ind
                row(count) = i(jik);val(count) = v(jik);
                count = count +1;
            end
        end
        
    elseif capdif(kk) > 0 %if habitat capacity has increased
        ind1 = i1(1); %first fish that lives in the habitat
        g1.speciesvec(ind1: ind1 + oldhc(kk) -1) = ic(i0); %replace everything you can
        i = i1(ismember(i0,row_in));
        v = time_in(ismember(row_in,i0));
        if any(i)
            n_add_ind = length(i);
            for jik = 1:n_add_ind
                row(count) = i(jik);
                val(count) = v(jik);
                count = count +1;
            end
        end
    end
    
end

%--- DISPERSAL --------------------------------------------------

for kk = 1:p1.nhabitats
    if any(g1.speciesvec(g1.fishhab ==kk)) % you can call this habitat for dispersal
        listind(kk) = 1;
    end
end

%then call dispersal (still will be empty places maybe, but that is
%now minimized) We don't call this in the first loop because we want to
%fill every available spot w/out dispersal first.

for  kk = 1:ndis %first disperse into habitats that exist in both landscapes
    habind = g1ind(disp_ind(kk)); %habitat
    dif = capdif(disp_ind(kk));  %difference in habitat capacity
    ind1 = disp_fish_ind(kk); %fish index of first fish in the habitat
    hc0 = oldhc(disp_ind(kk)); %old habitat capacity
   
    for iii = 1:dif %for each fish
        fish = ind1 + hc0 + iii - 1; %empty fish unit that needs to be filled
        [ p1,g1,fish_d ] = onlydispersal_protracted( p1,g1,fish,habind,listind); %disperse, and track which fish dispersed there for incipient speciation matrix       
       
        inc_ind = find(row == fish_d); %incipient speciation
        if any(inc_ind)
            n_inc = length(inc_ind);
            for jj = 1:n_inc %for however many incipient species columns there are
                row(count) = fish; 
                val(count) = val(inc_ind(jj));
                count = count +1;
            end
        end
        
    end
    
    if listind(habind) == 0 %allow dispersal from this habitat if it has now been filled
        listind(habind) = 1;
    end
end

%now for totally new habitats
newhabind = find(~ismember(g1.habitat,g0.habitat)); %totally new habitat
for kk = 1:length(newhabind)
    habind = newhabind(kk);
    hc1 = g1.habitatcapacity(habind); %new habitat capacity
    ind1 = find(g1.fishhab == habind,1); %find first instance
    for iii = 1:hc1
        fish = ind1 + iii - 1; %index of fish that needs to be replaced
        [ p1,g1,fish_d ] = onlydispersal_protracted( p1,g1,fish,habind,listind); %disperse, and track which fish dispersed there for incipient speciation matrix
  
        inc_ind = find(row == fish_d); %protracted speciation
        if any(inc_ind)
            n_inc = length(inc_ind);
            for jj = 1:n_inc %for however many incipient species columns there are
                row(count) = fish;
                val(count) = val(inc_ind(jj));
                count = count +1;
            end            
        end
    end
    
    listind(habind) = 1; %you can now choose from this habitat for dispersal
end

t_inc1 = int32(val(val~=0)); 
inc_ind1 = int32(row(val~=0)); 
end

