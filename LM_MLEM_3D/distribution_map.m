close all; 
clear all;
clc

x1 = csvread('E:\DoYeon\Document\5. Program\Reconstruction\4.MLEM_RECON\LM_MLEM_3D\Meeting_result\result_recon_ob7.csv');

figure(1); 
plot(x1(:,1),x1(:,4),'-b','linewidth', 4);
title('Intensity distribution in X-axis', 'fontsize', 16, 'fontname', 'times');
xlabel('X (m)', 'fontsize', 16, 'fontname', 'times');
ylabel('Intensity', 'fontsize', 16, 'fontname', 'times');
figure(2);
plot(x1(:,2),x1(:,4),'-r','linewidth', 4);
title('Intensity distribution in Z-axis', 'fontsize', 16, 'fontname', 'times');
xlabel('Z (m)', 'fontsize', 16, 'fontname', 'times');
ylabel('Intensity', 'fontsize', 16, 'fontname', 'times');
figure(3);
plot(x1(:,3),x1(:,4),'-g','linewidth', 4);
title('Intensity distribution in Y-axis', 'fontsize', 16, 'fontname', 'times');
xlabel('Y (m)', 'fontsize', 16, 'fontname', 'times');
ylabel('Intensity', 'fontsize', 16, 'fontname', 'times');
