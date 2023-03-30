% test the performance and give an example of the using of MC

clear all;
close all;
clc;

% initialize the input data structure
input= InitializeInput('C:\Users\纳豆够巷\Desktop\生物医学光子学实验\MC_dll\inputdatafile_1.m');
% call MC
[output] = MC(input);

%plot reflectance
dr = input.dr;
nr = input.nr;
rd_r = output.rd_r;
log_rd_r = log(rd_r);
rindex = (dr:dr:nr*dr)-dr/2;

figure('name','reflectance')

subplot(2,1,1)
plot(rindex,rd_r)
grid on
title('reflectance')
xlabel('distance[cm]')
ylabel('reflectance')

subplot(2,1,2)
plot(rindex,log_rd_r)
grid on
title('LOG(reflectance)')
xlabel('distance[cm]')
ylabel('LOG(reflectance)')
