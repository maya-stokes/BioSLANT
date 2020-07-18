function [ g ] = DispersalKernelTable( p,g )
%Calculates the probability of dispersal; with and without speciation; 
%
% Dispersal Kernel options: 
%   2dt : see Clarke et al., 1999
%   uniform: all dispersal within drainage basin has same probability
%   none: no dispersal 
%   neighbor: can only disperse to closest neighbor 
%

g.cdfP = zeros(p.nhabitats + 2, p.nhabitats);
g.cdfPind = zeros(p.nhabitats +1,p.nhabitats); %index of sorted cumulative distribution function (
g.cdfP_no_speciation = zeros(p.nhabitats+1,p.nhabitats);
g.cdfPind_no_speciation = zeros(p.nhabitats,p.nhabitats);

for i = 1:p.nhabitats
    distance = g.ddtable(:,i); %find dispersal distance in look-up table (just in number of ssteps, not units)    
    
    switch p.dispersalfunction
        
        case '2dt'
            K = p.p./(pi*p.u^2*(1+(distance./p.u).^2).^(p.p+1));
            % Lynch et al. (2011) version: K = p.p/(pi*p.u*(1+(distance.^2/p.u)).^(p.p+1); 
        case 'uniform'
            K = ones(p.nhabitats,1)./p.nhabitats; %uniform dispersal within basin
            K(distance == 0) = 0;
        case 'none'
            K = zeros(p.nhabitats,1); %can't disperse
            K(i) = 1;
        case 'neighbor' %can only disperse to neighbor
            K = zeros(p.nhabitats,1);
            ne = distance == 1; ne(i) = 1; nne = nnz(ne);
            K(ne) = 1./nne;
    end
    
    indr = distance == 0;  indr(i) = 0;
    K(indr) = 0; %take care of 1's or NaN's where distance == 0
    
    P_no_speciation = K.*g.habitatcapacity/sum(K.*g.habitatcapacity); % no chance of speciation (this is for only dispersal, when landscape changes)
    P = [(1-p.mu).*(K.*g.habitatcapacity/sum(K.*g.habitatcapacity)) ; p.mu]; %probability, sum == 1
    
    
    [P_nosp,ind_nosp] = sort(P_no_speciation,'ascend');
    [Psort,ind] = sort(P,'ascend'); %sort so that distance = 0 are in the front (and cdf stays 0).
   
    g.cdfP(:,i) =  [0; cumsum(Psort)]; %cumulative probability distribution
    g.cdfPind(:,i) = ind; %for proper indexing 
    g.cdfP_no_speciation(:,i) = [0; cumsum(P_nosp)];
    g.cdfPind_no_speciation(:,i) = ind_nosp;   
end


g.cdfP = g.cdfP(:); %vectorize for effiency in C code
g.cdfPind = g.cdfPind(:);
g.cdfP_no_speciation = g.cdfP_no_speciation(:);
g.cdfPind_no_speciation = g.cdfPind_no_speciation(:);

