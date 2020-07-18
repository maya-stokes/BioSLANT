%example script to call the functions used to calculate statistics on
%divide migration and capture size/frequency, and make a movie highlighting
%the captures with white polygons. 


runname = 'example_fault1.mat';
load(runname);
p.plotFaults = 1; Amin = 20; capdecay = 0.15; cap_interval = 1;frame_interval = 10;
p.plotCaptures = 1;
Stats = TransientStats(output,p,t,cap_interval,Amin,'example_fault1'); 
CapturesMovie(output,p,t,frame_interval,runname,Stats,capdecay)
