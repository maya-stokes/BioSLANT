%Species-area relationship. 
%To use, unzip and add "Example Steady State Output" to path as well as
%"Landscape Scenario 1 NCM" 

load('IC_U_LEM1a_TAU11_P2.5_U2IC_r197.mat'); 
load('LEM-1-a_1_.mat','g'); 
p.nhabitats = length(g.habitat);
[p,g] = habitatcapacity(p,g); %habitat capacity
nsp = zeros(20,100); alist = zeros(20,1); %initialize area and species lists

for i = 1:20 %number of basins
    hablist = find(g.basinind == i); %habitats in the basin
    fishlist = ismember(g.fishhab,hablist); 
    
    for j = 1:100 %number of stochastic realizations of the NCM
        nsp(i,j) = length(unique(speciesmat(fishlist,j))); 
    end
    
    alist(i) = max(g.A(g.habitat(g.basinind == i))); %drainage area
    
end

x = log10(alist); y = log10(mean(nsp,2));

b = regress(y,[ones(length(x),1) x]); %regress log basin richness and log species richness

figure; %plot data and regression
scatter(x,y); hold on; 

xx = min(x):max(x);
yy = b(1)+b(2).*xx; 

plot(xx,yy);
text(7.8,0.7,['log(S) = ', num2str(round(b(1),2)),' + ', num2str(round(b(2),2)),'log(A)']); 
xlabel('log(Drainage Area (km^2))'); ylabel('log(Species)'); 