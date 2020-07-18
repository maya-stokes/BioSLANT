function [ p,g,ind_new ] = onlydispersal_protracted( p,g,fish,fishind,listind)
%does dispersal with no speciation 

ind = 1+p.nhabitats*(fishind -1);
cdfind = g.cdfPind_no_speciation(ind:ind+p.nhabitats-1); %index vector for habitat you will disperse to

ind = 1+(p.nhabitats+1)*(fishind-1); %index for cdf
cdf = g.cdfP_no_speciation(ind:ind+p.nhabitats);


list2keep = listind; %habitats that you can pull from (there are fish there)
hab2use = find(list2keep); %habitat index

list2keep = ismember(cdfind,hab2use); %cdf index
cdf2use = cumsum(cdf(list2keep))./sum(cdf(list2keep)); %only consider places where there are fish and re-scale the probability
cdfhab = cdfind(list2keep);

num = rand; %choose random number
[~,i] = histc(num,cdf2use);
ind = find(g.fishhab == cdfhab(i)); %fish that live in index
list = g.speciesvec(ind);
list2use = list(list > 0);
ind2use = ind(list>0);
nfishunits = length(list2use);
nfish = randi(nfishunits); %randomly choose which fish unit within node will disperse
nnew = list2use(nfish); % fish species to use for replacement
ind_new = ind2use(nfish);
g.speciesvec(fish) = nnew; %replace with the species of the immigrant
end

