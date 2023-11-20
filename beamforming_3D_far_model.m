clc;
clear;
close all;
warning off;
addpath(genpath(pwd));
rng('default')

NUM_signal = 1;                             % 信号个数
NUM_zy     = 9;                             % 阵元个数，均匀线阵
NUM_zy_2   = 9;                             % 阵元个数，均匀线阵
samples    = 1000;                          % 采集信号的数量
theta      = [[20,10]];                     % 多目标信号的来波方向,夹角方向为以x轴夹角为第一夹角，以水平面的夹角为第二夹角
theta_pi   = theta.*pi/180;                 % 将角度转换为弧度制   
fil        = 1;                             % 与x轴夹角，方位角    
thta       = 2;                             % 与水平面夹角,俯仰角
A          = zeros(NUM_zy,NUM_signal*NUM_zy_2);      % 目标方向的导向矢量，用于模拟真实信号
SNR        = 10;                            % 信号信噪比 
cos_thetam     = zeros(1,NUM_signal);       % 角度转换器
sin_thetam     = zeros(1,NUM_signal);       % 角度转换器


for m=1:NUM_signal*NUM_zy_2
    cos_thetam = sin(theta_pi(1,thta))*cos(theta_pi(1,fil));                   % 求解方向矢量与x轴方向的cos值 
    sin_thetam = sin(theta_pi(1,thta))*sin(theta_pi(1,fil));                   % 求解方向矢量与x轴方向的cos值 
    A(:,m) = exp(j*pi*(cos_thetam*(0:NUM_zy-1)+sin_thetam*(m-1)));             % 其中一个有效信号的导向矢量
end 

Sn = randn(NUM_signal,samples)+j*randn(NUM_signal,samples);     % 模拟真实信号中的有效信号
Vn = randn(NUM_signal,samples)+j*randn(NUM_signal,samples);             % 模拟真是信号中的噪声信号
% Sn = repmat(Sn,NUM_zy_2,1);                                     % 将Sn信号进行垂直复制9次
% xt = A*Sn;                                     % xt就是实际接收到的信号数据
% 矩阵中创造矩阵
xt = zeros(NUM_zy*NUM_zy_2,samples);
for AA = 1:NUM_zy
    for BB = 1: NUM_zy_2
        num = (AA-1)*9+BB;
        xt(num,:) = A(AA,BB)*Sn+Vn/(10.^(SNR/10));
    end
end

Rjn = (xt*xt')/samples;                                         % 计算自身信号的自相关系数

fil_angle   =  -90:1:90;
theta_angle =  -90:1:90;

A_Direction = zeros(81,1);
for ii = 1:length(fil_angle)
   for jj = 1:length(theta_angle)
       Cos_Theta_Direction = sind(theta_angle(jj))*cosd(fil_angle(ii));
       Sin_Theta_Direction = sind(theta_angle(jj))*sind(fil_angle(ii));
       for kk = 1:NUM_zy
           for mm = 1:NUM_zy_2   
            num_2 = (kk-1)*9+mm;
            A_Direction(num_2) = exp(j*pi*(Cos_Theta_Direction*(kk-1)+Sin_Theta_Direction*(mm-1)));
           end
       end
%        beam_power(ii,jj) = abs(sum(sum(conj(A_Direction).*A)));       % (1).用于最简单的波束形成，导向矢量求其相关性
       beam_power(ii,jj) = abs(A_Direction'*Rjn*A_Direction);         % (2).用于传统的波束形成，求其空间谱
%        beam_power(ii,jj) = 1/abs(A_Direction'*inv(Rjn)*A_Direction);  % (3).用于MVDR的波束形成
%        beam_power(ii,jj) = sum(sum((A_Direction'*xt)));               % (4).用于传统的波束形成，求其空间谱与(2)类似

   end
end

beam_power=abs(beam_power);
max_value = max(beam_power(:));
% beam_power = beam_power/max_value;3D-
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
