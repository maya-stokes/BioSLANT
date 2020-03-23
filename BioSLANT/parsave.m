function parsave( fname,p,g)
%Save outputs from parfor runs of the NCM 
save(strcat(fname,'.mat'),'p','g');
end

