clear all;
% 信号参数设置
symbol_rate = 625e6;                %原始码元速率(625MHz)
sps = 32;                           %每符号采样点数(插值因子)
Fs = sps * symbol_rate;             %原始采样率
N = 1e5;                            %符号数量
Ts = 1/Fs;                          %原始采样周期
target_rate = 131.072 * 16 * 1e6;   %目标采样率
% 生成初始随机 BPSK 信号
rng(0);
msg = randi([0 1], 1, N);    % 生成随机数据
bmsg = 2 * msg - 1;          % BPSK 调制，将0映射为-1，1映射为1

% 滚降滤波器
rolloff_factor = 0.5;        % 滚降系数                  
span = 16;                   % 滚降点数                      
rcos_fir = rcosdesign(rolloff_factor, span, sps, 'sqrt');  % 平方根滚降滤波器
upmsg = upsample(bmsg, sps); % 上采样
rcosmsg = conv(upmsg, rcos_fir, 'same');  % 滤波成型

% 初始化变换后的信号数组
M = N;                      % 保持符号数量不变
y = zeros(1, N);            % 初始化输出信号数组
len = length(rcosmsg);
nco = [0.889,zeros(1,N-1)];

current_doppler = 31.25e3;
w = target_rate / (symbol_rate+current_doppler) / sps;
temp_D2 = 0;
temp_D1 = 0;
temp_D0 = 0; 
k = 1;
for i = 1 : len - 1
    nco(i + 1) = nco(i) - w;
    if(nco(i + 1) <= 0)
        u = nco(i) / w;
        nco(i + 1) = mod(nco(i + 1), 1);
        v3 = 1/6 * temp_D2 - 1/2*temp_D1 + 1/2 * temp_D0 + -1/6 * rcosmsg(i);
        v2 = 0 + 1/2 * temp_D1 - temp_D0 + 1/2 * rcosmsg(i);
        v1 = -1/6 * temp_D2 + temp_D1 - 1/2 * temp_D0 - 1/3 * rcosmsg(i);
        v0 = temp_D0;
        y(k) = ((v3 * u + v2) * u + v1) * u + v0;
        k=k+1;
    end
    temp_D2 = temp_D1;
    temp_D1 = temp_D0;
    temp_D0 = rcosmsg(i);
end
x = 1:1:len;
figure;
plot(x(1:1000),rcosmsg(1:1000),'-+',x(1:1000),y(2:1001),'-^');

