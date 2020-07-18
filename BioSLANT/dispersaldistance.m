function [distance,updist,downdist] = dispersaldistance( p,g,ind )

%function that calculates number of upstream and downstream steps between
%two habitat nodes in a landscape; draws on FastScape algorithm.

Ad8 = g.A;
[rec,~] = receiver(Ad8,p,g); %receivers
[ndon,don] = donor(rec,p,g); %donors
updist = zeros(p.nhabitats,1); %initialize upstream distance
downdist = zeros(p.nhabitats,1); %initialize downstream distance
ind1 = ind; 
upstreamdistance(ind,don,ndon); %call upstream distance function, this until you reach somebody with no donors (headwaters)
downstreamdistance(ind); %call downstream distance function, this will you go until your each somebody who is her own receiver (outlet)

    function upstreamdistance(ind,don,ndon)
        downdist(ind) = downdist(g.habitat == rec(ind));  
        for k = 1:ndon(ind) %for number of donors
            ijk = don(ind,k); %donor
            if ijk ~= g.habitat(ind)            
            %downdist(ind) = downdist(g.habitat == rec(ind));            
            hab = find(g.habitat == ijk);
            % add length of donor to receiver to distance already calculated for receiver to previous index (death index)
            updist(hab) = 1 + updist(ind); %length is distance between a pixel and it's receiver
            upstreamdistance(hab,don,ndon)
            end
        end
    end

    function downstreamdistance(ind)
        
        ijk = rec(ind); %index is the receiver
        hab = find(g.habitat == ijk); %habitat of the receiver
        downdist(hab) = 1 + downdist(ind);  % length(ind) + downdist(ind); %length, index to its receiver + downstream distance already calculated from first pixel to receiver
        
        if ijk == g.habitat(ind) %if you reach an indnd2ex that it is it's own receiver (outlet), end
           return
        else
            downstreamdistance(hab) %otherwise call again!
        end
        
    end

trunk = downdist > 0; %where downdist is distance downstream from dead fish node
trunkhab = g.habitat(trunk);
donhab = don(trunk,:); %donor of each trunk node
trunkhabnew = donhab(~ismember(donhab,trunkhab)); % linear index of donors 
indlist = find(ismember(g.habitat,trunkhabnew)); % index of donors in g.habitat list 
indlist = indlist(indlist~=ind1); 

for jj = 1:size(indlist,1) %for fish in tribs that need to go down then upstream 
    ind = indlist(jj);
    updist(ind) = 1;
    downdist(ind) = downdist(g.habitat == rec(ind));
    upstreamdistance(ind,don,ndon); 
end

% updist(ind1) = 0; downdist(ind1) = 0; 
distance = updist + p.dw*downdist; %weighted by p.dw 
end

