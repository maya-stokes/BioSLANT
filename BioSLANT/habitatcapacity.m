function [p,g] = habitatcapacity(p,g)
% Calculate habitat capacity for each node; 
% Create vectors for indexing in "newfishloop.c"
% Save the location of each organism 

% if habitat field doesn't already exit
if ~isfield(p,'nhabitats') || ~isfield(g,'habitat') 
    channels = g.A > p.Acf; %predesignate critical drainage area that a fish could live in
    g.habitat = find(channels); %linear index of the channels
    p.nhabitats = length(g.habitat); %number of habitats
end

areavec = g.A(g.habitat); %area of each habitat
a = areavec.^p.hab_dist; %relationship between habitat and area
habitatcapacity = ceil(p.Ch*(a./max(a))); %habitat capacity (amnt of resources) in each habitat

%make sure minimum number of fish is in each pixel 
add_hc = p.hab_dist_k - min(habitatcapacity); %make sure minimum number of fish is in each pixel
if add_hc > 0
    habitatcapacity = habitatcapacity + add_hc;
end

%Get closer to average fish density
add_hc = p.Ch - mean(habitatcapacity);
g.habitatcapacity = ceil(habitatcapacity + add_hc);

p.nfishunits = sum(g.habitatcapacity);

%assign habitat locatio to each organism
fishunits = zeros(p.nhabitats,max(habitatcapacity)); %

for i = 1:p.nhabitats
    fishunits(i,1:habitatcapacity(i)) = 1; %populate each free habitat resource with fish 
end

g.fishhab = zeros(p.nfishunits,1);

counter = 0;

for j = 1:p.nhabitats %each fish has individually assigned habitat in vector
    g.fishhab(counter + 1: counter + g.habitatcapacity(j)) = j; %assign each fish unit with a river pixel
    counter = counter + g.habitatcapacity(j);
end


g.speciesvec = zeros(p.nfishunits,1); %species identity of each organism
g.SUMhc = cumsum(g.habitatcapacity);  %for indexing in "newfishloop.c"
end




