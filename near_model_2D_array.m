clc;
clear;
close all;
warning off;
addpath(genpath(pwd));
rng('default')

mic_position = [30,1];                       % 声源位置,极坐标形式
NUM_signal = 1;                             % 信号个数
NUM_zy     = 9;                             % 阵元个数，阵列以（0,0）为原点，沿着x轴正向移动
samples    = 1000;                          % 采集信号的数量
% theta      = [30];                          % 多目标信号的来波方向
A          = zeros(NUM_zy,NUM_signal);      % 目标方向的导向矢量，用于模拟真实信号
thetam     = zeros(1,NUM_signal);           % 角度转换器
c          = 340;
f          = 1000;
lamba      = c/f;  
r          = mic_position(2);
theta      = mic_position(1);
for m=1:NUM_zy
%     thetam(m) = theta(m)*pi/180;                                % 将角度转换为弧度制
    r_m_square = r.^2+((lamba/2)*(m-1)).^2-lamba*(m-1)*r*cosd(theta); 
    r_m        = sqrt(r_m_square);
    A(m)       = exp(j*2*pi*f*(r-r_m)/c);             % 其中一个有效信号的导向矢量
end 

Sn = randn(NUM_signal,samples)+j*randn(NUM_signal,samples);     % 模拟真实信号中的有效信号
Vn = randn(NUM_zy,samples)+j*randn(NUM_zy,samples);             % 模拟真是信号中的噪声信号
xt = A*Sn+Vn;                                                   % xt就是实际接收到的信号数据
Rjn = (xt*xt')/samples;                                         % 计算自身信号的自相关系数