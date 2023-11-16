clc;
clear;
close all;
warning off;
addpath(genpath(pwd));
rng('default')


NUM_signal = 2;                             % 信号个数
NUM_zy     = 9;                             % 阵元个数
samples    = 1000;                          % 采集信号的数量
theta      = [30,60];                       % 多目标信号的来波方向
A          = zeros(NUM_zy,NUM_signal);      % 目标方向的导向矢量，用于模拟真实信号
thetam     = zeros(1,NUM_signal);           % 角度转换器

for m=1:NUM_signal
    thetam(m) = theta(m)*pi/180;                                % 将角度转换为弧度制    
    A(:,m) = exp(j*pi*sin(thetam(m))*(0:NUM_zy-1));             % 其中一个有效信号的导向矢量
end 

Sn = randn(NUM_signal,samples)+j*randn(NUM_signal,samples);     % 模拟真实信号中的有效信号
Vn = randn(NUM_zy,samples)+j*randn(NUM_zy,samples);             % 模拟真是信号中的噪声信号
xt = A*Sn+Vn;                                                   % xt就是实际接收到的信号数据
Rjn = (xt*xt')/samples;                                         % 计算自身信号的自相关系数

%% mvdr算法
ths=-90:1:90;
a=zeros(NUM_zy,1);
for k=1:length(ths)
    for i=1:NUM_zy
        a(i)=exp(j*pi*sin(ths(k)*pi/180)*(i-1))           % 遍历每个角度
    end
    beam1(k) = 1/(a'*inv(Rjn)*a);                         % 根据mvdr得到的结论，进行功率谱的计算，得到最大值便是想要的来波角度
end
beam1=abs(beam1);
beam1=beam1/max(beam1);
[V,I] = max(beam1);
err=abs(ths(I)-theta)/theta                               % 误差计算
% beam1 = 10*log10(beam1);
figure (1);
plot(ths,beam1,'b-')
hold on; 

%% 最简单的DOA估计
a_test = zeros(9,1);
y_result = zeros(1,length(ths));
ths=-90:1:90;
Rjn1 = (xt*xt')/samples;                                                % 计算信号的自相关系数

for jj = 1:length(ths)
    for jjj = 1:NUM_zy
        a_test(jjj) = exp(j*pi*sin(ths(jj)*pi/180)*(jjj-1));            % 遍历每个角度      
    end
                                                                        % 上述两种波束形成均可以，分别是最简单的波束形成和MVDR波束形成。
%     y_result(jj) = a_test'*Rjn1*a_test;
    y_result(jj) = norm(a_test'*xt)/norm(Sn);
end

y_result = abs(y_result);
y_result = y_result/(max(y_result));
[V1,I1] = max(y_result);
err1=abs(ths(I)-theta)/theta

figure (2);
plot(ths,y_result,'b-');

save R.mat ths beam1 err