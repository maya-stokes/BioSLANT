% time-series of species richness, speciation, extinction, diversification
% and river capture (e.g. Figure 5)

%To use, unzip and add "Example Transient Output" and "Landscape Scenario
%1" to path

nspecies = zeros(196,1);
speciation = zeros(195,1);
extinction = zeros(195,1);

nsp = zeros(100,1);
sp = zeros(100,1);
ex = zeros(100,1);

%find average number of species, speciations, and extinctions 
for i = 1:196
    load(['IC_U_LEM1a_TAU11_P2.5_U2_TR_',num2str(i+1),'.mat'],'speciesmat');
    
    for j = 1:100
        specieslist = unique(speciesmat(:,j)); %number of species
        
        if i > 1
            specieslist0 =unique(speciesmat0(:,j));
            sp(j) = sum(~ismember(specieslist,specieslist0)); %speciation
            ex(j) = sum(~ismember(specieslist0,specieslist)); %extinction
        end
        nsp(j) = length(specieslist);
    end
    
    if i > 1 %averages
        speciation(i) = mean(sp);
        extinction(i) = mean(ex);
    end
    nspecies(i) = mean(nsp);
    
    speciesmat0 = speciesmat;
end

%%
nlandscapes = 197; %find river captures
load(['LEM-1-a_1_.mat'],'g'); g0  = g; ncap = zeros(nlandscapes-1,1);
for i = 2:nlandscapes
    load(['LEM-1-a_',num2str(i),'_.mat'],'g');
    [~,g1ind,g0ind] = intersect(g.habitat,g0.habitat); %any habitat node that switches basins between timesteps is defined as a river capture
    capind = g.basinind(g1ind) - g0.basinind(g0ind)~=0;
    ncap(i-1,1) = nnz(capind);
    g0 = g;
    
end

%%
figure; %plot results
subplot(5,1,1);
plot(nspecies(2:end));
ylabel('Number of Species'); 
subplot(5,1,2);
plot(speciation(2:end));
ylabel('Speciations');
subplot(5,1,3);
plot(extinction(2:end));
ylabel('Extinctions');
subplot(5,1,4);
plot(speciation(2:end) - extinction(2:end));
ylabel('Diversification');
subplot(5,1,5);
bar(ncap(2:end))
ylabel('Captures'); 
xlabel('10 Generations');