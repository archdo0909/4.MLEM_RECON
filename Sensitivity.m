close all;
clear all;
clc


Det.Num_scatter = 64;
Det.Num_absorber = 64;
Det.WIDTH = 8; %cm
Det.RESO = 1; %cm
Det.Distance_scatter_absorber = -65; %cm\
Det.nGrid = 64*64;
result.detect = [];
result.coordinate = [];

for i=1:Det.nGrid
   index=observation_mat(i,Det);
   result.detect = [result.detect; index GetXYZFromScatterIndex(index(1),Det) GetXYZFromAbsorberIndex(index(2),Det)];
end
csvwrite('detector_index.csv',result.detect);

function index=observation_mat(compton_ind, Det)
    index = compton_ind - 1;
    ind_scat = fix(index/Det.Num_scatter);
    ind_abso = rem(index,Det.Num_scatter);
    
    index = [ind_scat ind_abso];
end
function xyz=GetXYZFromScatterIndex(ig,Det)
    index = ig;
    indx = fix(index/Det.WIDTH);
    indy = rem(index,Det.WIDTH);
    
    x = GetXYPosition(indx,Det.WIDTH,Det.RESO);
    y = GetXYPosition(indy,Det.WIDTH,Det.RESO);
    
    xyz = [x*10 y*10 0];
end
function xyz=GetXYZFromAbsorberIndex(ig,Det)
    indx = fix(ig/Det.WIDTH);
    indy = rem(ig,Det.WIDTH);
    
    x = GetXYPosition(indx,Det.WIDTH,Det.RESO);
    y = GetXYPosition(indy,Det.WIDTH,Det.RESO);
    
    xyz = [x*10 y*10 Det.Distance_scatter_absorber];
end
function position=GetXYPosition(index, width, resolution)
    position=(index*resolution-width/2)+(resolution)/2;
end