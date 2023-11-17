clc;
clear;
close all;
warning off;
addpath(genpath(pwd));
rng('default')


NUM_signal = 1;                             % 信号个数
NUM_zy     = 3;                             % 阵元个数
samples    = 1000;                          % 采集信号的数量
theta      = [[45,60]];                     % 多目标信号的来波方向,夹角方向为以x轴夹角为第一夹角，以水平面的夹角为第二夹角
theta_pi   = theta.*pi/180;                 % 将角度转换为弧度制   
fil        = 1;                             % 与x轴夹角    
thta       = 2;                             % 与水平面夹角 
A          = zeros(NUM_zy,NUM_signal);      % 目标方向的导向矢量，用于模拟真实信号
cos_thetam     = zeros(1,NUM_signal);       % 角度转换器
SNR        = 10;                            % 信号信噪比 

for m=1:NUM_signal
    cos_thetam(m) = sin(theta_pi(m,fil))*cos(theta_pi(m,thta));                                % 求解方向矢量与x轴方向的cos值 
    A(:,m) = exp(j*pi*cos_thetam(m)*(0:NUM_zy-1));             % 其中一个有效信号的导向矢量
end 

Sn = randn(NUM_signal,samples)+j*randn(NUM_signal,samples);     % 模拟真实信号中的有效信号
Vn = randn(NUM_zy,samples)+j*randn(NUM_zy,samples);             % 模拟真是信号中的噪声信号
xt = A*Sn+Vn/1000;                                              % xt就是实际接收到的信号数据
Rjn = (xt*xt')/samples;                                         % 计算自身信号的自相关系数

fil_angle   =  -90:1:90;
theta_angle =  0:1:90;

A_Direction = zeros(NUM_zy,1);
% beam_power = zeros(length(fil_angle),length(theta_angle));
for ii = 1:length(fil_angle)
   for jj = 1:length(theta_angle)
       Cos_Theta_Direction = sind(fil_angle(ii))*cosd(theta_angle(jj));
       for kk = 1:NUM_zy
       A_Direction(kk) = exp(j*pi*Cos_Theta_Direction*(kk-1));
       end
%        beam_power(ii,jj) = A_Direction'*;
%        beam_power(ii,jj) = A_Direction'*Rjn*A_Direction;
       beam_power(ii,jj) = 1/norm(A_Direction'*inv(Rjn)*A_Direction);
%        beam_power(ii,jj) = norm(A_Direction'*xt)/norm(Sn);
   end
end

beam_power=abs(beam_power);
max_value = max(beam_power(:));
beam_power = beam_power/max_value;
% beam_power = 20*log10(beam_power);
[r, c] = find(beam_power == max(beam_power(:)));
% 
[X, Y] = meshgrid(theta_angle,fil_angle);
% 绘制三维图形
waterfall(X, Y,beam_power);          % 画出立体网状图
ylabel('fil');
xlabel('thta');
zlabel('功率谱');
% plot3(fil_angle,theta_angle,beam_power,'o');

