function [ndon, don] = donor(rec,p,g)
%Make a donor list for each node (don) and record how many donors each node
%has (ndon) - From Jean Braun's Fastscape Workshop November, 2015

ndon = zeros(p.nhabitats,1); %initialize
don = zeros(p.nhabitats,8);

for k = 1:p.nhabitats
    ind = g.habitat(k);
    dij = find(rec == ind); %find everyone for whom I am a receiver
    ndon(k) = length(dij); %number of donors I have
    don(k,1:length(dij)) = g.habitat(dij); %add to donor list  
end
end


