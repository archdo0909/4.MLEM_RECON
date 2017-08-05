clear all;
close all;
clc


function [X_sp,Y_sp,Z_sp]=CamtoWorld(X_w,Y_w,Z_w,X_s,Y_s,Z_s,alpha)
    
    X_sp = X_s - X_w;
    Y_sp = Y_s - Y_w;
    Z_sp = Z_s - Z_w;
    
    R=[cos(alpha) -sin(alpha);
     sin(alpha) cos(alpha)];
    
    [X_cp, Y_cp] = R*[X_w,Y_w]; 
   
end
function radian = ToAngle(x1,y1,z1,x2,y2,z2)
    cos=(x1*x2 + y1*y2 + z1*z2)/(sqrt(x1^2+y1^2+z1^2)*sqrt(x2^2+y2^2+z2^2);
    radian = arc(cos);
end
function radian = ToRadian(degree)
    radian = degree/180*pi;
end