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
theta      = [[30,20,1.2]];                   % 多目标信号的来波方向,夹角方向为以x轴夹角为第一夹角，以z轴面的夹角为第二夹角,之后为距离
theta_pi   = theta.*pi/180;                 % 将角度转换为弧度制   
fil        = 1;                             % 与x轴夹角，方位角    
thta       = 2;                             % 与水平面夹角,俯仰角
r          = theta(3);
A          = zeros(NUM_zy,NUM_signal*NUM_zy_2);      % 目标方向的导向矢量，用于模拟真实信号
SNR        = 10;                            % 信号信噪比 
cos_thetam     = zeros(1,NUM_signal);       % 角度转换器
sin_thetam     = zeros(1,NUM_signal);       % 角度转换器
c          = 340;
f          = 1000;
lamba      = c/f;

for m=1:NUM_zy
    for m_x = 1:NUM_zy_2
    cos_thetam = sin(theta_pi(1,thta))*cos(theta_pi(1,fil))*(m-1);                   % 求解方向矢量与x轴方向的cos值 
    sin_thetam = sin(theta_pi(1,thta))*sin(theta_pi(1,fil))*(m_x-1);                   % 求解方向矢量与x轴方向的cos值 
    fen_zi     = sqrt(((m-1).^2+(m_x-1).^2));
    cos_lambda = (cos_thetam+sin_thetam)/fen_zi; 
    
    r_di_side  = sqrt(((lamba/2)*(m-1)).^2+((lamba/2)*(m_x-1)).^2);
    r_m_square = r.^2+(r_di_side).^2-2*r*r_di_side*cos_lambda; 
    r_m        = sqrt(r_m_square);
    
    if m == 1 && m_x == 1   %由于存在原点位置的计算问题，因此需要认为补偿
        r_m = r;
    end
    A(m,m_x)       = exp(1i*2*pi*f*(r-r_m)/c);             % 其中一个有效信号的导向矢量
%     A(:,m) = exp(j*pi*(cos_thetam*(0:NUM_zy-1)+sin_thetam*(m-1)));             % 其中一个有效信号的导向矢量
    end
end 

Sn = randn(NUM_signal,samples)+1i*randn(NUM_signal,samples);     % 模拟真实信号中的有效信号
Vn = randn(NUM_signal,samples)+1i*randn(NUM_signal,samples);             % 模拟真是信号中的噪声信号
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

fil_angle   =  0:1:90;
theta_angle =  0:1:90;
R_all       =  0:0.01:1.5;

A_Direction = zeros(81,1);
for ii = 1:length(fil_angle)
   for jj = 1:length(theta_angle)
       for rr = 1:length(R_all)
           for kk = 1:NUM_zy
               for mm = 1:NUM_zy_2 
                   Cos_Theta_Direction = sind(theta_angle(jj))*cosd(fil_angle(ii))*(kk-1);
                   Sin_Theta_Direction = sind(theta_angle(jj))*sind(fil_angle(ii))*(mm-1);
                   fen_zi_Direction    = sqrt(((kk-1).^2+(mm-1).^2));
                   cos_lambda_Direction   = (Cos_Theta_Direction+Sin_Theta_Direction)/fen_zi_Direction;
                   
                   r_di_side_Direction  = sqrt(((lamba/2)*(kk-1)).^2+((lamba/2)*(mm-1)).^2);
                   r_m_square_Direction = R_all(rr).^2+(r_di_side_Direction).^2-2*R_all(rr)*r_di_side_Direction*cos_lambda_Direction; 
                   r_m_Direction        = sqrt(r_m_square_Direction);
                   
                   if kk == 1 && mm == 1                                    %由于存在原点位置的计算问题，因此需要认为补偿
                        r_m_Direction = R_all(rr);
                   end
                   num_2 = (kk-1)*9+mm;
                   A_Direction(num_2)   = exp(1i*2*pi*f*(R_all(rr)-r_m_Direction)/c);             % 其中一个有效信号的导向矢量    
               end
           end
%            beam_power(ii,jj,rr) = abs(sum(sum(conj(A_Direction).*A)));       % (1).用于最简单的波束形成，导向矢量求其相关性
           beam_power(ii,jj,rr) = abs(A_Direction'*Rjn*A_Direction);         % (2).用于传统的波束形成，求其空间谱
%          beam_power(ii,jj) = 1/abs(A_Direction'*inv(Rjn)*A_Direction);  % (3).用于MVDR的波束形成
%          beam_power(ii,jj) = sum(sum((A_Direction'*xt)));               % (4).用于传统的波束形成，求其空间谱与(2)类似

       end
   end
end

beam_power=abs(beam_power);
max_value = max(beam_power(:));
beam_power = beam_power/max_value;
beam_power = 20*log10(beam_power);
[maxValue, maxIndex] = max(beam_power(:)); % 获取最大值及其在矩阵中的线性索引
[X, Y, Z] = ind2sub(size(beam_power), maxIndex); % 将线性索引转换为三维坐标


B = reshape(beam_power, [], size(A, 100)); % 将A转换为一个二维矩阵B，其中每一行对应A中的一个二维平面
disp(['最大值：', num2str(maxValue)]);
disp(['最大值坐标：(', num2str(X), ',', num2str(Y), ',', num2str(Z), ')']);


