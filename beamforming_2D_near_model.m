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

%% mvdr算法
ths=0:1:80;
R_all = 0.5:0.01:2;
a=zeros(NUM_zy,1);
for k=1:length(ths)
    for jj =1:length(R_all)
        for i=1:NUM_zy
            r_m_square_2 = R_all(jj).^2+((lamba/2)*(i-1)).^2-lamba*(i-1)*R_all(jj)*cosd(ths(k)); 
            r_m_2        = sqrt(r_m_square_2);
            a(i)         = exp(j*2*pi*f*(R_all(jj)-r_m_2)/c);             % 其中一个有效信号的导向矢量
        end                      
          beam1(k,jj) = 1/(a'*inv(Rjn)*a);       % 根据mvdr得到的结论，进行功率谱的计算，得到最大值便是想要的来波角度    
    end
%     w = (inv(Rjn)*a)/((a')*inv(Rjn)*a);
%     baeam1(k) = norm(conj(w')*a);

end
beam1=abs(beam1);
% beam1=beam1/max(beam1);
% [V,I] = max(beam1);
% err=abs(ths(I)-theta)/theta;                               % 误差计算
% beam1 = 10*log10(beam1);

[X, Y] = meshgrid(R_all,ths);
% 绘制三维图形
waterfall(X, Y,beam1);          % 画出立体网状图
ylabel('ths');
xlabel('R_all');
zlabel('功率谱');
hold on;
