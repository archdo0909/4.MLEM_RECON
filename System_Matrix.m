function System_matrix()

x1 = csvread('compton_scatter_0_1.csv');
x2 = csvread('compton_absorber_0_1.csv');
x3 = csvread('compton_angle_0_1.csv');
x4 = csvread('compton_cnt_0_1.csv');
x5 = csvread('rm_data_0_1.csv');

compton.scatter = [];
compton.absorber = [];
compton.angle = [];
compton.cnt = 0;
for i=1:size(x1,1)
    compton.scatter = [compton.scatter; x1(i,2:4)];
end
for i=1:size(x2,1)
    compton.absorber = [compton.absorber; x2(i,2:4)];
end
for i=1:size(x3,1)
    compton.angle = [compton.angle; x3(i,1)];
end
compton.cnt = x4;

%compton = matfile(compton.mat);
rm.CENTER = [0 0];
rm.WIDTH = 2;
rm.HEIGHT = 2;
rm.RESO = 0.1;
rm.nGrid = round(rm.WIDTH*rm.HEIGHT/(rm.RESO^2));
rm.data = zeros(1,rm.nGrid);
result.pixel = [];
result.revise = [];
for i=1:rm.nGrid
   rm.data(1,i)=x5(i);
end

%data
Det.Num_scatter = 64;
Det.Num_absorber = 64;
Det.WIDTH = 1; %cm
Det.RESO = 1; %cm
Det.Distance_scatter_absorber = 80; %cm\
%Num_angle = 90;
%Num_data = Det.Num_scatter*Det.Num_absorber*Num_angle;
E_o = 0.66;

% for i = 1:rm.nGrid
%     xyz=GetThetaPhiFromIndex(i,rm);
%     result.pixel = [result.pixel; i xyz];
% end
% csvwrite('result_pixel.csv',result.pixel);

for it=1:20
    compton.revise_factor = zeros(1,compton.cnt-1);
    for i=1:1:compton.cnt-1
       for k=1:1:rm.nGrid
           compton.revise_factor(1,i)=compton.revise_factor(1,i) + System_mat(i,k)*rm.data(1,k);
       end
       %result.pixel = [result.pixel; i compton.revise_factor(1,i)];
    end
% csvwrite('result_pixel_A.csv',result.pixel);
% revise = 0;
%     for j=1:1:rm.nGrid    
%         for i=1:1:compton.cnt-1
%             if compton.revise_factor(1,i) ~= 0
%                 revise=revise + System_mat(i,j)/compton.revise_factor(1,i);
%             end
%         end
%           result.revise = [result.revise; j revise];
%           revise = 0;
%     end
%     csvwrite('result_revise.csv',result.revise);
    
    A=zeros(1,rm.nGrid);
    for j=1:1:rm.nGrid    
       for i=1:1:compton.cnt-1
           if compton.revise_factor(1,i) ~= 0
            A(1,j)=A(1,j)+ System_mat(i,j)/compton.revise_factor(1,i);
           end
       end
            rm.data(1,j)=rm.data(1,j)*A(1,j);
    end
end

for j=1:1:rm.nGrid
   rm.data(1,j)=round(rm.data(1,j)); 
end

ShowImage(rm);
% %scatter.coordinate 
% for i=0:Det.Num_scatter-1
%     compton.scatter = [compton.scatter; GetXYZFromScatterIndex(i,Det)];
% end
% %absorber.coordinate
% for i=0:Det.Num_absorber-1
%     compton.absorber = [compton.absorber; GetXYZFromAbsorberIndex(i,Det)];
% end

% for j=1:rm.nGrid
%    for i=0:Det.Num_scatter-1
%        System_mat(i,j)
%    end
% end

function [system_mat]=System_mat(indexforcompton,indexforpixel)
var = 1.5;
P = 1;
pixel = GetThetaPhiFromIndex(indexforpixel,rm);
%compton.scatter = GetXYZFromScatterIndex(indexforcompton,Det);
%compton.absorber = GetXYZFromAbsorberIndex(indexforcompton,Det);
beta_frompixel = Cal_Scatter_ang(pixel,compton.scatter(indexforcompton,1:3),compton.absorber(indexforcompton,1:3));
beta = acos(compton.angle(indexforcompton,1));
norm_compton = [0 0 1];
%center = [0 0 0];
theta_g = Cal_Angle(compton.scatter(indexforcompton,1:3)-compton.absorber(indexforcompton,1:3), norm_compton);
%theta_g = Cal_Angle((Transform_cart(pixel)-center),norm_compton);
if abs(beta_frompixel - beta) < P*var
    system_mat = Klein_Nishina(compton.angle(indexforcompton,1), E_o)*abs(cos(theta_g))*(1/sqrt(2*pi*var))*exp( (-1/2)*( (beta_frompixel - beta)/var )^2  )...
        /(Cal_distance(Transform_cart(pixel),compton.scatter(indexforcompton,1:3)))^2;
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
    [x y z] = sph2cart(toRadian(pixel(1)),toRadian(90-pixel(2)),1);
    cart = [x y z];
    %center = [0 0 0];
    X = cart-scatter;
    Y = scatter - absorber;
    Angle = Cal_Angle(X,Y);
end
function Angle=Cal_Angle(x,y)
    cos_radian=(x(1)*y(1) + x(2)*y(2) + x(3)*y(3))/(sqrt(x(1)^2 + x(2)^2 + x(3)^2)*sqrt(y(1)^2+y(2)^2+y(3)^2));
    Angle = acos(cos_radian);
    %Angle = cos_radian;
end
function [Pro_scat]=Klein_Nishina(theta, E_o)
    r_e = 2.817940 * 10^-15;  %cm
    m_c_2 = 0.511;
    %theta =-> cos value
    %a = E_o/m_c_2;
%     Pro_scat = (r_e^2/2) * (1/ (1 + a(1 - cos(theta))))^2 *...
%         (1 + cos(theta)*cos(theta) + a^2*(1-cos(theta))^2/(1+a*(1-cos(theta)))); 

%%%% the value is radian
%     E_prime = E_o*m_c_2/(E_o*(1-cos(theta)) + m_c_2);
%     Pro_scat = (r_e^2/2)*(E_prime/E_o)^2*( (E_prime/E_o) + (E_o/E_prime) -sin(theta)^2  );
    
%%% the value is cos 
    E_prime = E_o*m_c_2/(E_o*(1-theta) + m_c_2);
    Pro_scat = (r_e^2/2)*(E_prime/E_o)^2*( (E_prime/E_o) + (E_o/E_prime) -sin(acos(theta))^2  );
    
end
function Angle=GetThetaPhiFromIndex(ig,rm)
    indx = fix((ig-1)/rm.WIDTH);
    indy = rem((ig-1),rm.HEIGHT);
    
    x = (rm.CENTER(1) - rm.WIDTH/2 + rm.RESO/2) + indx*rm.RESO ;
    y = (rm.CENTER(2) - rm.HEIGHT/2 + rm.RESO/2) + indy*rm.RESO; 
    Angle = [x*90 y*90];
end
%scatter coordinate
function xyz=GetXYZFromScatterIndex(ig,Det)
    indx = fix(ig/Det.WIDTH);
    indy = rem(ig,Det.WIDTH);
    
    x = GetXYPosition(indx,Det.WIDTH,Det.RESO);
    y = GetXYPosition(indy,Det.WIDTH,Det.RESO);
    
    xyz = [x y 0];
end
function xyz=GetXYZFromAbsorberIndex(ig,Det)
    indx = fix(ig/Det.WIDTH);
    indy = rem(ig,Det.WIDTH);
    
    x = GetXYPosition(indx,Det.WIDTH,Det.RESO);
    y = GetXYPosition(indy,Det.WIDTH,Det.RESO);
    
    xyz = [x y Det.Distance_scatter_absorber];
end
function position=GetXYPosition(index,width,resolution)
    position = (index*resolution - width/2)+resolution/2;
end
function radian = toRadian(degree)
    radian = degree/180*pi;
end
function Angle = toDegree(Radian)
    Angle = Radian/pi*180;
end
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
end


