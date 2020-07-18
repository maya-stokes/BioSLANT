function [ g ] = dispersaldistancetable( p,g )
%First calculate drainage area 

g.ddtable = zeros(p.nhabitats);

for i = 1:p.nhabitats
   g.ddtable(:,i) = dispersaldistance(p,g,i); %make a lookup table of dispersal distance
end

end

