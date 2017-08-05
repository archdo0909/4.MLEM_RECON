
clc
result.data = [];
load('E:\DoYeon\Document\5. Program\Reconstruction\4.MLEM_RECON\LM_MLEM_3D\result\20170716\miniron_data3.mat');

max = -1000;
for j=1:td.nGrid
    if max < td.data(1,j)
        max = td.data(1,j);
        disp(j);
        
    end
end
for j = 1:td.nGrid
    td.data(1,j)=td.data(1,j)/max;
end
for j=1:td.nGrid
    J = rem(j,(td.WIDTH*td.HEIGHT/(td.RESO^2)));
    xy=GetXYFromDataIndex(J,td);
    z = fix(j/(td.WIDTH*td.HEIGHT/(td.RESO^2))) + 3;
    pixel = [xy(1) xy(2) z];
    plot3(pixel(1),pixel(2),pixel(3),'--ks',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerFaceColor',[td.data(1,j),0.1,1-td.data(1,j)]);
%     scatter3(pixel(1),pixel(2),pixel(3),'d','MarkerFaceColor',[td.data(1,j),0.1,1-td.data(1,j)]...
%     ,'MarkerFaceAlpha',0.2);
    grid on;
    hold on;
    %result.coordinates = [result.coordinates; j pixel];
end
for j=1:td.nGrid
    J = rem(j,(td.WIDTH*td.HEIGHT/(td.RESO^2)));
    xy=GetXYFromDataIndex(J,td);
    z = fix(j/(td.WIDTH*td.HEIGHT/(td.RESO^2))) + 3;
    pixel = [xy(1) xy(2) z];
    
    result.data = [result.data; pixel td.data(1,j)];

end
csvwrite('E:\DoYeon\Document\5. Program\Reconstruction\4.MLEM_RECON\LM_MLEM_3D\result\20170716\miniron_data3.csv',result.data);
function xy=GetXYFromDataIndex(index,rm)
    ig = index -1;
    %Get coordinates of x and y from the data index
    %Get x,y index
    indy = rem(ig,(rm.WIDTH/rm.RESO));
    indx = fix(ig/(rm.WIDTH/rm.RESO));

    x = GetXYPosition(indx,rm.WIDTH,rm.RESO);
    y = GetXYPosition(indy,rm.WIDTH,rm.RESO);
    xy=[x y];
end
function position=GetXYPosition(index, width, resolution)
    position=(index*resolution-width/2)+(resolution)/2;
end