close all;
clear all;
clc

x1 = csvread('E:\DoYeon\Document\5. Program\Reconstruction\4.LM_MLEM\data\20170626\Sullivan2.csv');
x2 = csvread('E:\DoYeon\Document\5. Program\Reconstruction\4.LM_MLEM\data\20170626\(-9,-3,0)_10^6.csv');
%x1 = csvread('sample2.csv');
% x2 = csvread('(0,0,3)_test.csv');
%Radiation Energy
E_o = 0.66; %MeV
%mass of electron
m=9.10938356*10^-31;  %kg 
%speed of light
c=299792458; %m/s
%eV->J

eV=1.60218*10^-19; %J
%foscal 
f=10; %cm
%coordinate reconstruction
rm.CENTER = [0 0];
rm.RESO = 0.01;
rm.WIDTH = 2; 
rm.HEIGHT = 2;
rm.nGrid = round(rm.WIDTH*rm.HEIGHT/(rm.RESO^2));
rm.data = zeros(1,rm.nGrid);
rm.Index = zeros(1,rm.nGrid);
rm.Sph_elevation = 90;
rm.Sph_azimuth= 90;
rm.Sph_RESO = 1;
Det.Num_scatter = 64;
Det.Num_absorber = 64;
Det.WIDTH = 8; %cm
Det.RESO = 1; %cm
Det.Distance_scatter_absorber = -65; %mm\
Det.nGrid = 64*64;
compton.scatter=[];
compton.absorber=[];
compton.angle=[];
compton.cnt = 1;
result.Index = [];
robot_po1=[0 0 -5];
scope = 10^-1;
camera_direction1 = [0 0 1];
Index_compton = zeros(2,Det.nGrid*2);
for i=1:1:(size(x1,1)-1)  
    if(x1(i,1)<64 && x1(i+1,1) > 63 && x1(i,5) < 0.66)   %select first layer coincidence
        if((x1(i,5)+x1(i+1,5)) <= 0.66 && (x1(i,5)+x1(i+1,5)) >= 0.6555)  %compared deposited energy
            compton.scatter(compton.cnt, 1) = x1(i,1);         
            compton.scatter(compton.cnt, 2) = x1(i,2)*scope+robot_po1(1);
            compton.scatter(compton.cnt, 3) = x1(i,3)*scope+robot_po1(2);
            compton.scatter(compton.cnt, 4) = x1(i,4)*scope+robot_po1(3);
            compton.scatter(compton.cnt, 5) = x1(i,5);
            compton.scatter(compton.cnt, 6) = camera_direction1(1);
            compton.scatter(compton.cnt, 7) = camera_direction1(2);
            compton.scatter(compton.cnt, 8) = camera_direction1(3);
            
            compton.absorber(compton.cnt, 1) = x1(i+1,1);         
            compton.absorber(compton.cnt, 2) = x1(i+1,2)*scope+robot_po1(1);
            compton.absorber(compton.cnt, 3) = x1(i+1,3)*scope+robot_po1(2);
            compton.absorber(compton.cnt, 4) = x1(i+1,4)*scope+robot_po1(3);
            compton.absorber(compton.cnt, 5) = x1(i+1,5);
            compton.cnt = compton.cnt + 1;
        end
    end
end
j = compton.cnt;
for i=1:compton.cnt-1
    index = Det.Num_scatter * compton.scatter(i,1) + compton.absorber(i,1) - (Det.Num_absorber - 1);
    Index_compton(1,index) = Index_compton(1,index) + 1;
    Index_compton(2,index) = i;
end
robot_po2 = [10 0 5];
camera_direction2 = [-1 0 0];
for i=1:1:(size(x2,1)-1)  
    if(x2(i,1)<64 && x2(i+1,1) > 63 && x2(i,5) < 0.66)   %select first layer coincidence
        if((x2(i,5)+x2(i+1,5)) <= 0.66 && (x2(i,5)+x2(i+1,5)) >= 0.6555)  %compared deposited energy
            compton.scatter(compton.cnt, 1) = x2(i,1);         
            compton.scatter(compton.cnt, 2) = -x2(i,4)*scope+robot_po2(1);
            compton.scatter(compton.cnt, 3) = x2(i,3)*scope+robot_po2(2);
            compton.scatter(compton.cnt, 4) = x2(i,2)*scope+robot_po2(3);
            compton.scatter(compton.cnt, 5) = x2(i,5);
            compton.scatter(compton.cnt, 6) = camera_direction2(1);
            compton.scatter(compton.cnt, 7) = camera_direction2(2);
            compton.scatter(compton.cnt, 8) = camera_direction2(3);
            
            compton.absorber(compton.cnt, 1) = x2(i+1,1);         
            compton.absorber(compton.cnt, 2) = -x2(i+1,4)*scope+robot_po2(1);
            compton.absorber(compton.cnt, 3) = x2(i+1,3)*scope+robot_po2(2);
            compton.absorber(compton.cnt, 4) = x2(i+1,2)*scope+robot_po2(3);
            compton.absorber(compton.cnt, 5) = x2(i+1,5);
            compton.cnt = compton.cnt + 1;
        end
    end
end
for i=j:compton.cnt-1
    index = Det.Num_scatter * compton.scatter(i,1) + compton.absorber(i,1) - (Det.Num_absorber - 1);
    index2 = index + 64*64;
    Index_compton(1,index2) = Index_compton(1,index2) + 1;
    Index_compton(2,index2) = i;
end
%Calculating cos(Theta)
for i=1:1:compton.cnt-1 
    %angle=1-((0.511875*compton.scatter(i,5))/(0.66*compton.absorber(i,5)));  %cos
    angle = 1 - 0.511*(1/compton.absorber(i,5) - 1/(E_o));
     compton.angle(i,1)=angle;
end
Recon_R=1;
E_ang = 1;

% Index_compton = zeros(2,Det.nGrid);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:compton.cnt-1
%     index = Det.Num_scatter * compton.scatter(i,1) + compton.absorber(i,1) - (Det.Num_absorber - 1);
%     Index_compton(1,index) = Index_compton(1,index) + 1;
%     Index_compton(2,index) = i;
% end

for i=1:Det.nGrid*2
   result.Index =[result.Index; i Index_compton(1,i) Index_compton(2,i)];
end
csvwrite('compton_index.csv', result.Index);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Shperical reconstruction
% for i=1:1:compton.cnt-1
%    for x = -rm.Sph_azimuth:rm.Sph_RESO:rm.Sph_azimuth
%       for y = -rm.Sph_elevation:rm.Sph_RESO:rm.Sph_elevation
%         [a,b,c]=sph2cart(toRadian(x), toRadian(90-y),Recon_R);
%         u1 = [a,b,c];
%         d = compton.scatter(i, 2:4)-compton.absorber(i,2:4);
%         if(Inner_product(u1, d) <= compton.angle(i,1)*cos(toRadian(E_ang)) + sin(acos(compton.angle(i,1)))*sin(toRadian(E_ang))...
%            && Inner_product(u1, d) >= compton.angle(i,1)*cos(toRadian(E_ang)) - sin(acos(compton.angle(i,1)))*sin(toRadian(E_ang)))
%             r = sqrt( (a/Recon_R)^2 + (b/Recon_R)^2 );
%             r_prime = sqrt(1 - sqrt(1 - r^2));
%             %if(r~=0)
%             new_x = (a/Recon_R)*r_prime/r;
%             new_y = (b/Recon_R)*r_prime/r;
%             R_prime = sqrt(new_x^2 + new_y^2);
%             Index = GetIndexFromXY(new_x, new_y, rm.WIDTH,rm.RESO);
%             rm.Index(1,Index)=rm.Index(1,Index)+1;
%             if(rm.Index(1,Index) == 1)
%                 rm.data(1,Index)=rm.data(1,Index)+1;
%             end
%         end
%       end
%    end
%    rm.Index = zeros(1,rm.nGrid);
% end
% %%%%%%%%%%%%%%%%LM-MLEM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for it=1:2
%     compton.revise_factor = zeros(1,compton.cnt-1);
%     for i=1:1:compton.cnt-1
%        for k=1:1:rm.nGrid
%            compton.revise_factor(1,i)=compton.revise_factor(1,i) + System_mat(i,k,compton,rm,E_o)*rm.data(1,k);
%        end
%     end 
%     A=zeros(1,rm.nGrid);
%     for j=1:1:rm.nGrid    
%        for i=1:1:compton.cnt-1
%            if compton.revise_factor(1,i) ~= 0
%             A(1,j)=A(1,j)+ System_mat(i,j,compton,rm,E_o)/compton.revise_factor(1,i);
%            end
%        end
%             rm.data(1,j)=rm.data(1,j)*A(1,j);
%     end
% end
% for j=1:rm.nGrid
%     rm.data(1,j)=round(rm.data(1,j));
% end

%3d-reconstruction
td.CENTER = [0 0];
td.RESO = 1;
td.WIDTH = 20; 
td.HEIGHT = 20;
td.nGrid = round(td.WIDTH*td.HEIGHT/(td.RESO^2));
td.data = zeros(1,td.nGrid);
td.sen = zeros(1,td.nGrid);
td.Index = zeros(1,td.nGrid);
result.test = [];
result.test_2 = [];
%robot_po = [0 0 -5];
%scope = 10^-1;
%camera_direction = [ 0 0 1 ];

%%%%% y = 0;
for i=1:compton.cnt-1
    if (compton.angle(i,1) <= 1 && compton.angle(i,1) >= -1)
        for j=1:td.nGrid
            xy=GetXYFromDataIndex(j,td);
            pixel = [xy(1) xy(2) 5];
            %d = compton.scatter(i, 2:4)-compton.absorber(i,2:4);
            scatter = compton.scatter(i, 2:4);
            absorber = compton.absorber(i,2:4);
            d = (scatter - absorber);
            u1 = (pixel - scatter);
            %d2 = [d(1) d(3)]; 
            if Cal_Angle(u1, d) <= compton.angle(i,1)*cos(toRadian(E_ang))+ sin(acos(compton.angle(i,1)))*sin(toRadian(E_ang))...
               && Cal_Angle(u1, d) >= compton.angle(i,1)*cos(toRadian(E_ang)) - sin(acos(compton.angle(i,1)))*sin(toRadian(E_ang)) 
            %if Cal_Angle(u1,d) <= acos(compton.angle(i,1)) + toRadian(E_ang) && Cal_Angle(u1,d) >= acos(compton.angle(i,1)) - toRadian(E_ang)
                td.data(1,j) = td.data(1,j) + 1;
            end 
        end
    end
end
figure(100);
ShowImage(td);
shading interp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%LM
%%%%%%%%%%%%%%%%Sensitivity%%%%%%%%%%%%%%%%%%%%%
%S = zeros(1,td.nGrid);
% 
% %%%%%%%%%%%%%%%%%%
% A_scat = 1*1; %cm
% A_absor = 1*1; %cm
% % %%%%%%%%%%%%%%%%%%
% Z_det = 1; %cm
% X_air = 1432.8 * 16 / ( 8*(8+1)*(11.319 - log(8)));
% %X_CsI = 1432.8 * (132*0.5 + 125.9*0.5) / ( (55*0.5 + 53*0.5)*(55*0.5 + 53*0.5+1)(11.319 - log(55*0.5 + 53*0.5)))
% X_CsI = 171.5; % g cm^-2  
% norm_compton = [0 0 1];
% sen_abso=0;
% sen_scat=0;
% for j = 1:td.nGrid
%    xy=GetXYFromDataIndex(j,td);
%    pixel_y_zero = [xy(1) -3  xy(2)];
%    
%    for N=0:Det.Num_scatter-1
%        scatter = GetXYZFromScatterIndex(N,Det)*scope + robot_po;
%        cos_01 = Cal_Angle((pixel_y_zero - scatter), norm_compton);
%        R_01 = Cal_distance(pixel_y_zero,scatter);
%        for K=0:Det.Num_absorber-1
%           absorber = GetXYZFromAbsorberIndex(K,Det)*scope + robot_po; 
%           cos_12 = Cal_Angle((scatter-absorber), norm_compton);
%           cos_compton = Cal_Angle( (scatter-absorber) , (pixel_y_zero-scatter) );
%           R_12 = Cal_distance(scatter,absorber);
%           sen_abso = sen_abso + exp(-R_12/X_air)*(1-exp(-A_absor/X_CsI))*Klein_Nishina(cos_compton,E_o)*A_absor*cos_12/R_12^2;
%        end
%        sen_scat = sen_scat + exp(-R_01/X_air)*(1-exp( -A_scat/X_CsI))*A_scat*cos_01/(4*pi*R_01^2)*sen_abso;
%    end
%     
%    S(1,j)=sen_scat;
%    result.Pro = [result.Pro; j S(1,j)];
%    sen_abso=0;
%    sen_scat=0;
%     
% end

 sen_max = -1000;
for j=1:td.nGrid
    xy=GetXYFromDataIndex(j,td);
    pixel = [xy(1) xy(2) 5];
       for N=0:Det.Num_scatter-1
           scatter = GetXYZFromScatterIndex(N,Det)*scope + robot_po1;
           for K=0:Det.Num_absorber-1
               absorber = GetXYZFromAbsorberIndex(K,Det)*scope + robot_po1; 
               td.sen(1,j)=td.sen(1,j)+System_mat_Sen(scatter,absorber,j,td,E_o,camera_direction1);
           end
           if sen_max < td.sen(1,j)
              sen_max = td.sen(1,j);
           end
       end
end
for j = 1:td.nGrid
    xy=GetXYFromDataIndex(j,td);
    pixel = [xy(1) xy(2) 5];
       for N=0:Det.Num_scatter-1
           %scatter = GetXYZFromScatterIndex(N,Det)*scope + robot_po2;
           scat = GetXYZFromScatterIndex(N,Det);
           scatter_turn = [-scat(3) scat(2) scat(1)];
           scatter = scatter_turn*scope + robot_po2;
           for K=0:Det.Num_absorber-1
               abso = GetXYZFromAbsorberIndex(K,Det)*scope + robot_po2;
               abso_turn = [-abso(3) abso(2) abso(1)];
               absorber = abso_turn*scope + robot_po2;
               td.sen(1,j)=td.sen(1,j)+System_mat_Sen(scatter,absorber,j,td,E_o,camera_direction2);
           end
           if sen_max < td.sen(1,j)
              sen_max = td.sen(1,j);
           end
       end
    
end
for j = 1:td.nGrid 
   td.sen(1,j)=td.sen(1,j)/sen_max; 
end
% figure(50)
% ShowGrahp(td);
% shading interp;


%%%%%%%%%%%%%%%%%%%%%%%constant z%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for j = 1:td.nGrid
%     xy=GetXYFromDataIndex(j,td);
%     pixel = [xy(1) xy(2) 5 ];
%     % (0 0 0) : 1.00
%     for k =0:Det.Num_scatter-1
%         xyz=GetXYZFromScatterIndex(k,Det);
%         scatter = [xyz(1) xyz(2) xyz(3)]*scope + robot_po;
%         S(1,j) = S(1,j) + ((1-exp(-4.51*1))/Cal_distance(pixel, scatter)^3);
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%constant y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for j = 1:td.nGrid
%     xy=GetXYFromDataIndex(j,td);
%     pixel = [xy(1) 0 xy(2) ];
%     % (0 0 0) : 1.00
%     S(1,j)=(1.5/Cal_distance(pixel,robot_po))^2*Cal_Angle(pixel-robot_po, camera_direction);
%     result.Pro=[result.Pro; j xy(1) 0 xy(2) S(1,j)];
% end

%csvwrite('result_Pro.csv',result.Pro);
% for j=1:td.nGrid
%     for i=1:Det.nGrid
%         S(1,j)=S(1,j)+System_mat_Sen(i,j,Det,td,E_o,robot_po,scope);
%     end
% end
%csvwrite('result_Pro.csv',result.Pro);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%Iteration%%%%%%%%%%%%%%%%%%%%%%%
compton_cnt=0;
for i=1:Det.nGrid*2
    if Index_compton(1,i) > 0
        compton_cnt = compton_cnt + 1;
    end
end
max_data = -1000;
% for j = 1:td.nGrid
%     if max_data < td.data(1,j)
%         max_data=td.data(1,j);
%     end
% end
% for it = 1:1:1
%     %max = -100;
%     for j = 1:td.nGrid 
%         B=0;
%             for i=1:Det.nGrid*2
%                 A = 0;
%                 if (Index_compton(1,i) > 0)
%                     if(compton.angle(Index_compton(2,i),1) <= 1 && compton.angle(Index_compton(2,i),1) >= -1)
%                     %if (compton.angle(i,1) <= 1 && compton.angle(i,1) >= -1)
%                         for k = 1:td.nGrid 
%                             A = A + System_mat_2d(Index_compton(2,i),k,compton,td,E_o)*td.data(1,k);
%                         end
%                         result.test = [result.test; i A];
%                         B = B + Index_compton(1,i)*System_mat_2d(Index_compton(2,i),j,compton,td,E_o)/A;
%                     end
%                 end
%             td.data(1,j) = td.data(1,j)*B;
%             end
%     end
%     figure(it);
%     ShowImage(td);
%     shading interp;
% end
% csvwrite('testtest.csv',result.test);
for it=1:10
    compton.revise_factor = zeros(1,compton.cnt-1);
    max = -100;
    for i=1:Det.nGrid*2
        if (Index_compton(1,i) > 0)
            if (compton.angle(Index_compton(2,i),1) <= 1 && compton.angle(Index_compton(2,i),1) >= -1)
                for k=1:td.nGrid
                    compton.revise_factor(1,Index_compton(2,i))=compton.revise_factor(1,Index_compton(2,i))...
                        + (System_mat_2d(Index_compton(2,i),k,compton,td,E_o)*td.data(1,k));
                end
            end
            result.test = [result.test; i compton.revise_factor(1,Index_compton(2,i))];
        end
    end
    
    A=zeros(1,td.nGrid);
    for j=1:1:td.nGrid   
        for i=1:Det.nGrid*2
           if (Index_compton(1,i) > 0)
              if (compton.angle(Index_compton(2,i),1) <= 1 && compton.angle(Index_compton(2,i),1) >= -1)
                  if compton.revise_factor(1,Index_compton(2,i)) ~= 0
                    A(1,j)=A(1,j)+ (Index_compton(1,i)*System_mat_2d(Index_compton(2,i),j,compton,td,E_o)/compton.revise_factor(1,Index_compton(2,i)));
                  end
              end
           end
        end
        result.test_2 = [result.test_2; j, A(1,j)];
         td.data(1,j)=td.data(1,j)*A(1,j)/td.sen(1,j);
         
%          if max < td.data(1,j)
%             max = td.data(1,j); 
%          end
    end
%     thread = max - 20*sqrt(max);
%     disp(thread);
%     disp(max);
%     for j=1:td.nGrid
%         if td.data(1,j) < thread
%             td.data(1,j) = 0;
%         end
%     end
    
    figure(it);
    ShowImage(td);
    shading interp;
    csvwrite('test_1.csv',result.test);
    csvwrite('test_2.csv',result.test_2);
end



% for it=1:2
%     compton.revise_factor = zeros(1,compton.cnt-1);
%     for i=1:1:compton.cnt-1
%         if (compton.angle(i,1) <= 1 && compton.angle(i,1) >= -1)
%            for k=1:1:td.nGrid
%                compton.revise_factor(1,i)=compton.revise_factor(1,i) + (System_mat_2d(i,k,compton,td,E_o,robot_po,scope)*td.data(1,k));
%            end
%           % result.data_it = [result.data_it; i compton.revise_factor(1,i)];
%         end
%     end 
%     
%     A=zeros(1,td.nGrid);
%     for j=1:1:td.nGrid    
%        for i=1:1:compton.cnt-1
%            if (compton.angle(i,1) <= 1 && compton.angle(i,1) >= -1)
%                if compton.revise_factor(1,i) ~= 0
%                 A(1,j)=A(1,j)+ (System_mat_2d(i,j,compton,td,E_o,robot_po,scope)/compton.revise_factor(1,i));
%                end
%            end
%        end
%             td.data(1,j)=td.data(1,j)*A(1,j)*S(1,j);
%     end
%     figure(it);
%     ShowImage(td);
%     shading interp;
% end
% for it=1:2
%     compton.revise_factor = zeros(1,compton.cnt-1);
%     for i=1:1:compton.cnt-1
%         if (compton.angle(i,1) <= 1 && compton.angle(i,1) >= -1)
%            for k=1:1:td.nGrid
%                compton.revise_factor(1,i)=compton.revise_factor(1,i) + (System_mat_2d(i,k,compton,td,E_o,robot_po,scope)*td.data(1,k));
%            end
%           % result.data_it = [result.data_it; i compton.revise_factor(1,i)];
%         end
%     end 
%     
%     A=zeros(1,td.nGrid);
%     for j=1:1:td.nGrid    
%        for i=1:1:compton.cnt-1
%            if (compton.angle(i,1) <= 1 && compton.angle(i,1) >= -1)
%                if compton.revise_factor(1,i) ~= 0
%                 A(1,j)=A(1,j)+ (System_mat_2d(i,j,compton,td,E_o,robot_po,scope)/compton.revise_factor(1,i));
%                end
%            end
%        end
%             td.data(1,j)=td.data(1,j)*A(1,j)*S(1,j);
%     end
%     figure(it);
%     ShowImage(td);
%     shading interp;
% end
%csvwrite('Compton_revise_factor.csv',result.data_it);
% for j = 1:td.nGrid 
%    result.data = [result.data; j td.data(1,j)]; 
% end
%  csvwrite('result_data_it.csv',result.data_it);
% 
% figure(13);
% ShowImage(td);
% shading interp;
function ShowImage(rm)
    xmin=rm.CENTER(1)-rm.WIDTH/2;
    xmax=rm.CENTER(1)+rm.WIDTH/2-rm.RESO;
    ymin=rm.CENTER(2)-rm.HEIGHT/2;
    ymax=rm.CENTER(2)+rm.HEIGHT/2-rm.RESO;
    
    [X,Y]=meshgrid(xmin:rm.RESO:xmax, ymin:rm.RESO:ymax);
    sizeX = size(X);
    data=reshape(rm.data,sizeX(1),sizeX(2));
    surf(X,Y,data);
    
end
function [system_mat]=System_mat_2d(indexforcompton,indexforpixel,compton,td,E_o)
    var = 3;
    P = 2;  
    xy=GetXYFromDataIndex(indexforpixel,td);
    pixel_y_zero = [xy(1) xy(2) 5];
    scatter = compton.scatter(indexforcompton, 2:4);
    absorber = compton.absorber(indexforcompton,2:4);
    beta_frompixel = Cal_Angle(pixel_y_zero-scatter,scatter-absorber);
    beta = acos(compton.angle(indexforcompton,1));
    norm_compton = compton.scatter(indexforcompton, 6:8);
    %norm_compton = [0 0 1];
    %center = [0 0 0];
    cos_theta_m = Cal_Angle((pixel_y_zero - scatter), norm_compton);
    cos_theta_g = Cal_Angle((scatter-absorber), norm_compton);
    %theta_g = Cal_Angle((Transform_cart(pixel)-center),norm_compton);
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if abs(beta_frompixel - beta) < P*var
%         system_mat = Klein_Nishina(compton.angle(indexforcompton,1), E_o)*abs(cos_theta_m)*(1/(sqrt(2*pi)*var))*exp( (-1/2)*( (acos(beta_frompixel) - beta)/var )^2  )...
%                      /((Cal_distance(pixel_y_zero,scatter))^2); 
%     else 
%         system_mat = 0;
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m_c_2 = 0.511;
    E_prime = E_o*m_c_2/(E_o*(1-compton.angle(indexforcompton,1)) + m_c_2);
    %171.5
    system_mat = Klein_Nishina(compton.angle(indexforcompton,1),E_o)*(1 - exp(-1/(1.243*(compton.absorber(indexforcompton,5))*171.5)))...
        /((Cal_distance(scatter,absorber))^2 * cos_theta_g);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function [system_mat]=System_mat_Sen(scatter,absorber,indexforpixel,td,E_o,camera_direction)
    var = 3.0;
    P = 2;
    %%%%%%%%%%%%%%%%%%
    A_scat = 1*1; %cm
    A_absor = 1*1; %cm
    %%%%%%%%%%%%%%%%%%
    xy=GetXYFromDataIndex(indexforpixel,td);
    pixel_y_zero = [xy(1) xy(2) 5];
    norm_compton = [camera_direction(1) camera_direction(2) camera_direction(3)];
    cos_theta_g = Cal_Angle((scatter-absorber), norm_compton);
    cos_theta_m = Cal_Angle((pixel_y_zero - scatter), norm_compton);
    beta_frompixel = Cal_Angle(pixel_y_zero-scatter,scatter-absorber);
    beta = acos(beta_frompixel);
    %index = observation_mat(indexforcompton,Det);
    %scatter = GetXYZFromScatterIndex(index(1),Det)*scope + robot_po;
    %absorber = GetXYZFromAbsorberIndex(index(2),Det)*scope + robot_po;
    %theta_g = Cal_Angle((Transform_cart(pixel)-center),norm_compton);
%     if abs(beta_frompixel - beta) < P*var
%         system_mat = Klein_Nishina(cos_compton, E_o)*abs(cos_theta_m)*abs(cos_theta_m)*(1/(sqrt(2*pi)*var))*exp( (-1/2)*( (acos(beta_frompixel) - beta)/var )^2  )...
%             /((Cal_distance(pixel_y_zero,scatter))^2); 
%           %  /((Cal_distance(pixel_y_zero,scatter))^2*(Cal_distance(scatter,absorber))^2);
%     else 
%         system_mat = 0;
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %m_c_2 = 0.511;
    %E_prime = E_o*m_c_2/(E_o*(1-cos_compton) + m_c_2);
    %system_mat = Klein_Nishina(cos_compton,E_o)/((Cal_distance(pixel_y_zero,scatter))^2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     m_c_2 = 0.511;
    %E_prime = E_o*m_c_2/(E_o*(1-compton.angle(indexforcompton,1)) + m_c_2);
    %171.5
    E_2=1/((-beta_frompixel+1)/m_c_2 + 1/E_o); 
    system_mat = Klein_Nishina(beta_frompixel,E_o)*(1 - exp(-1/(1.243*(E_2)*171.5)))...
        /((Cal_distance(scatter,absorber))^2* cos_theta_g);
%     if abs(beta_frompixel - beta) < P*var
%         system_mat = Klein_Nishina(beta_frompixel, E_o)*abs(cos_theta_m)*(1/(sqrt(2*pi)*var))*exp( (-1/2)*( (acos(beta_frompixel) - beta)/var )^2  )...
%                      /((Cal_distance(pixel_y_zero,scatter))^2); 
%     else 
%         system_mat = 0;
%     end


end
function [system_mat]=System_mat(indexforcompton,indexforpixel,compton,rm,E_o)
    var = 5.0;
    P = 3;
    pixel = GetThetaPhiFromIndex(indexforpixel,rm);
    beta_frompixel = Cal_Scatter_ang(pixel,compton.scatter(indexforcompton,2:4),compton.absorber(indexforcompton,2:4));
    beta = acos(compton.angle(indexforcompton,1));
    norm_compton = [0 0 1];
    %center = [0 0 0];
    theta_g = Cal_Angle(compton.scatter(indexforcompton,2:4)-compton.absorber(indexforcompton,2:4), norm_compton);
    %theta_g = Cal_Angle((Transform_cart(pixel)-center),norm_compton);
    if abs(beta_frompixel - beta) < P*var
        system_mat = Klein_Nishina(compton.angle(indexforcompton,1), E_o)*10^4*abs(cos(theta_g))*(1/sqrt(2*pi*var))*exp( (-1/2)*( (beta_frompixel - beta)/var )^2  )...
            /(Cal_distance(Transform_cart(pixel),compton.scatter(indexforcompton,2:4)))^2;
    else 
        system_mat = 0;
    end
end
function [distance]=Cal_distance(x,y)
    distance = sqrt( (x(1) - y(1))^2 + (x(2) - y(2))^2 + (x(3) - y(3))^2 );
end
function xyz = Transform_cart(pixel)
    R=1;
    [x,y,z] = sph2cart(toRadian(pixel(1)),toRadian(90 - pixel(2)),R);
    xyz = [x,y,z];
end
function Angle=Cal_Scatter_ang(pixel,scatter,absorber)
    [x, y, z] = sph2cart(toRadian(pixel(1)),toRadian(90-pixel(2)),1);
    cart = [x, y, z];
    %center = [0 0 0];
    X = cart-scatter;
    Y = scatter - absorber;
    Angle = Cal_Angle(X,Y);
end
function Angle=Cal_Angle(x,y)
    cos_radian=(x(1)*y(1) + x(2)*y(2) + x(3)*y(3))/(sqrt(x(1)^2 + x(2)^2 + x(3)^2)*sqrt(y(1)^2+y(2)^2+y(3)^2));
    %Angle = acos(cos_radian);
    Angle = cos_radian;
end
function Angle=Cal_Angle_2d(x,y)
    cos_radian=(x(1)*y(1) + x(2)*y(2))/(sqrt(x(1)^2 + x(2)^2)*sqrt(y(1)^2+y(2)^2));
    %Angle = acos(cos_radian);
    Angle = cos_radian;
end
function [Pro_scat]=Klein_Nishina(theta, E_o)
    %r_e = 2.817940 * 10^-13;  %m   %theta => cosin
    %r_e =1;
    %m_c_2 = 0.511;
    %E_prime = E_o*m_c_2/(E_o*(1-theta) + m_c_2);
    %Pro_scat = (r_e^2/2)*(E_prime/E_o)^2*( (E_prime/E_o) + (E_o/E_prime) -sin(acos(theta))^2  ); 
    a = E_o/0.511;
    Pro_scat = (1/(1+a*(1-theta)))^2*((1 + theta^2)/2)*( 1+ (a^2*(1-theta)^2/((1+theta^2)*(1+a*(1-theta))))   );
end
function Angle=GetThetaPhiFromIndex(ig,rm)
    indx = fix((ig-1)/rm.WIDTH);
    indy = rem((ig-1),rm.HEIGHT);
    
    x = (rm.CENTER(1) - rm.WIDTH/2 + rm.RESO/2) + indx*rm.RESO ;
    y = (rm.CENTER(2) - rm.HEIGHT/2 + rm.RESO/2) + indy*rm.RESO; 
    Angle = [x*90 y*90];
end
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
function ig=GetIndexFromXY(x,y,width,resolution)
%     indx = (x - (resolution/2) + width/2)/resolution;
%     indy = (y - (resolution/2) + width/2)/resolution;
%     
%     ig = fix(indx * width + (indy + 1)); 
 
    indx = fix((x-resolution/2+width/2)/resolution);
    indy = fix((y-resolution/2+width/2)/resolution);
    ig=fix(indx*width/resolution+indy);
end
function A=Inner_product(x, y)
    
    inner = x(1)*y(1) + x(2)*y(2) + x(3)*y(3);
    length = sqrt( x(1)^2 + x(2)^2 + x(3)^2 )*sqrt( y(1)^2 + y(2)^2 + y(3)^2);
    %%%% cos
    A = inner / length;
end
function radian = toRadian(degree)
    radian = degree/180*pi;
end
function Angle = toDegree(Radian)
    Angle = Radian/pi*180;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    index = ig;
    indx = fix(index/Det.WIDTH);
    indy = rem(index,Det.WIDTH);
    
    x = GetXYPosition(indx,Det.WIDTH,Det.RESO);
    y = GetXYPosition(indy,Det.WIDTH,Det.RESO);
    
    xyz = [x*10 y*10 Det.Distance_scatter_absorber];
end






