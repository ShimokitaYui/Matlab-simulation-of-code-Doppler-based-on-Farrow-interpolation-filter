# Matlab-simulation-of-code-Doppler-based-on-Farrow-interpolation-filter
# 拉格朗日Farrow滤波器与码多普勒模拟
## 重采样过程
### 目标
- 原始信号: $x(m)$ 采样周期: $T_s$
- 目标信号: $y(k)$ 采样周期: $T_i$
### 步骤
1. $x(t)$ 采样至$x(mT_s)$
2. $x(mTs)$数模转换变为$y(t)$ : $y(t)=\sum_m{x(mT_s)h(t-mT_s)}$
3. $y(t)$采样至: $y(kT_i)$: $y(kT_i)=\sum_m{x(mT_s)h(kT_i-mT_s)}$
## Farrow算法
### 变量定义
$kT_i-mT_s=Ts(k\frac{T_i}{T_s}-m)=Ts(i+\mu_k)$  
- $m_k=\lfloor k\frac{T_i}{T_s} \rfloor$
- $\mu_k=k\frac{Ti}{Ts}-m_k$
- $i=m_k-m$
### 算法表达式
理想滤波器:    
$$y(kT_i)=\sum_i{x\left[(m_k-i)T_s\right]h\left[(i + \mu_k)T_s\right]}$$  
由于插值用的点数不可能是无穷的所以可以得到实际有限点滤波器表达如下所示:  
$$y(kT_i)=\sum_{i=I_1}^{I_2}{x\left[(m_k-i)T_s\right]h\left[(i + \mu_k)T_s\right]}$$  
## 拉格朗日多项式插值
拉格朗日插值法表达式:  
$$L_n(x)=\sum_{j=0}^{n-1}{y_ip_j(x)}$$  
$$p_j(x)=\prod_{i \in B_k(i\ne k)}{\frac{x-x_i}{x_k-x_i}}$$  
采样拉格朗日插值法来代替Farrow中的滤波器系数。当$t\in \left[m_kT_s, (m_k+1)T_s\right]$时,利用$m_k$周围的点进行拉格朗日插值法得到一下表达式:  
$$y(t)=\sum_{i=I_1}^{I_2}C_ix\left[m_k-i\right]=\sum_{i=I_1}^{I_2}C_ix((m_k-i)T_s)$$  
根据拉格朗日插值法可以得到$C_i$:  
$$C_i=h(i+\mu_k)=\prod_{j=I_1,j\ne i}^{I_2}{\frac{t-t_j}{t_i-t_j}}$$
将Farrow的变量带入可以得到:    
$$C_i(\mu_k)=\prod_{j=I_1,j\ne i}^{I_2}{\frac{j+\mu_k}{j-i}}$$  
因此当I固定时可以得到固定的系数,当点数为4时可以得到系数为:  
$$C_{-2}=\frac{1}{6}\mu_k^3-\frac{1}{6}\mu_k$$  
$$C_{-1}=-\frac{1}{2}\mu_k^3+\frac{1}{2}\mu_k^2+\mu_k$$  
$$C_{0}=\frac{1}{2}\mu_k^3-\mu_k^2-\frac{1}{2}\mu_k+1$$  
$$C_{1}=-\frac{1}{6}\mu_k^3+\frac{1}{2}\mu^2-\frac{1}{3}\mu_k$$  
## 码多普勒模拟
使用拉格朗日Farrow滤波器进行码多普勒仿真,首先需要对原始信号进行高采样点数的上采样,对原始信号的采样点进行插值抽取获得目标码率的信号。下面是matlab仿真程序:
```matlab
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
    w = target_rate / (symbol_rate+current_doppler) / sps;
    temp_D2 = temp_D1;
    temp_D1 = temp_D0;
    temp_D0 = rcosmsg(i);
end
```
- 程序中目标输出的采样频率是定值,也符合实际系统中DA转换器件的特性,那么码速率是如何发生变化的,从本质上来说其实改变的是输入信号的采样频率,也就是时钟。`w = target_rate / (symbol_rate+current_doppler) / sps;`起到了关键作用。无论码速率是多少,一个码元过采样32个点并成形滤波后得到的数据是一样的,那么码率体现在什么地方？这就要考虑到`rcosmsg`的下标也就是`i`代表的含义,`i`代表的是输入的采样频率,`for i = 1 : len - 1`可以看成时钟,原始码率和发生码多普勒的后的码率所对应的`rcosmsg`是一样的,但是`i`所对应的时钟周期是不一样,也就是说该代码的码多普勒本质上是通过改变时钟来实现的。  
- 但是在FPGA的实现中这是很难实现的,那么假设`i`代表的时钟周期不变,这段代码还有意义吗？  
    可以这么看,输出信号的`y(k)`代表的是输出的第k个周期的信号大小,那么`i`的物理意义的变化与`y(k)`的大小值并无关系。那么是不是可以说即使是我时钟周期不变,只要我的`y(k)`的值大小一一对应了就没有问题？确实,无论你`i`对应的周期是否发生变化,你的`y(k)`都不会发生变化,但是在具体实现中生成`y(k)`的速度会发生变化,在`i`的周期不变的情况下可能会出现`y(k)`的值堆积或不够用的情况。那么只要解决了这个问题就可以实现完全的码多普勒。  
