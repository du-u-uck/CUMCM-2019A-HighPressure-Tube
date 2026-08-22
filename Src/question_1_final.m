function question_1_final
clc; close all;
format longG;

%% ====================== 1. 全局参数 ======================
L = 500; D = 10;
V = pi * (D/2)^2 * L;          % 油管容积 mm³
dA = 1.4;
AA = pi*(dA/2)^2;               % 入口小孔面积 mm²
C = 0.85;                       % 流量系数
P_pump = 160;                   % 油泵恒压 MPa
rho100 = 0.850;                 % 100MPa基准密度 mg/mm³
T_close = 10;                   % 单向阀固定关闭时长 ms

dt_opt = 0.05;                  % 寻优仿真步长 ms
t_sim_opt = 20000;              % 寻优总仿真时长 ms
t_steady_start = 15000;         % 从15秒后开始算稳态

%% ====================== 2. 预计算等间隔查表======================
% 读取原始E-P数据
data_E = readmatrix("附件3-弹性模量与压力.xlsx");
P_raw = data_E(:,1);
E_raw = data_E(:,2);

P_step = 0.01;
P_table = (0:P_step:200)';
N_table = length(P_table);

% 插值得到对应E数组
E_table = interp1(P_raw, E_raw, P_table, 'linear', 'extrap');

% 数值积分计算rho_table（dP = E/rho * drho）
rho_table = zeros(N_table, 1);
idx_100 = round(100 / P_step) + 1;
rho_table(idx_100) = rho100;

% 向高压积分
for i = idx_100+1 : N_table
    dP = P_step;
    E_mid = (E_table(i-1) + E_table(i)) / 2;
    rho_mid = rho_table(i-1);
    drho = rho_mid / E_mid * dP;
    rho_table(i) = rho_table(i-1) + drho;
end
% 向低压积分
for i = idx_100-1 : -1 : 1
    dP = P_step;
    E_mid = (E_table(i) + E_table(i+1)) / 2;
    rho_mid = rho_table(i+1);
    drho = rho_mid / E_mid * dP;
    rho_table(i) = rho_table(i+1) - drho;
end

% 封装极速查表函数
get_rho = @(P) lookup_linear(P_table, rho_table, P_step, P);
get_E = @(P) lookup_linear(P_table, E_table, P_step, P);

rho160 = get_rho(160);          

%% ====================== 3. 寻优：稳态100MPa最优Ton ======================
fprintf('===== 寻优单向阀开启时长 =====\n');
T_low = 0.20; T_high = 0.40;
iter_max = 25;
loss_history = zeros(iter_max,1);

for iter = 1:iter_max
    T_mid = (T_low + T_high) / 2;
    [P_mean, loss] = sim_steady(T_mid, 100, t_sim_opt, dt_opt, t_steady_start, ...
        V, C, AA, P_pump, rho160, get_rho, get_E);
    loss_history(iter) = loss;
    fprintf('迭代%2d | Ton=%.5f ms | 稳态均值=%.4f MPa\n', iter, T_mid, P_mean);
    
    if P_mean > 100
        T_high = T_mid;
    else
        T_low = T_mid;
    end
end
Ton_opt_100 = (T_low + T_high) / 2;
fprintf('\n【最优开启时长】Ton = %.4f ms\n\n', Ton_opt_100);

%% ====================== 4. 生成20000ms输出文件 ======================
fprintf('正在生成20000ms压力曲线并保存至 result.xlsx...\n');
dt_out = 0.1;
t_total_out = 20000;
[t_out, P_out] = sim_full(Ton_opt_100, 100, t_total_out, dt_out, ...
    V, C, AA, P_pump, rho160, get_rho, get_E);

% 每10ms采样
t_sample = (0:10:t_total_out)';
P_sample = interp1(t_out, P_out, t_sample, 'linear');
T_out = table(t_sample, P_sample, 'VariableNames', {'时间_ms','压力_MPa'});
writetable(T_out, 'result.xlsx');
disp('已输出0~20000ms每10ms压力至 result.xlsx');

%% ====================== 5. 求解150MPa稳态Ton ======================
fprintf('\n===== 求解150MPa稳态开启时长 =====\n');
T15_low = 0.6; T15_high = 0.9;
for iter=1:20
    Tm = (T15_low + T15_high)/2;
    [P_mean, ~] = sim_steady(Tm, 150, t_sim_opt, dt_opt, t_steady_start, ...
        V, C, AA, P_pump, rho160, get_rho, get_E);
    if P_mean > 150
        T15_high = Tm;
    else
        T15_low = Tm;
    end
end
Ton_150_stable = (T15_low + T15_high)/2;
fprintf('150MPa稳态Ton = %.4f ms\n', Ton_150_stable);

%% ====================== 6. 分阶段升压策略（2s, 5s, 10s） ======================
fprintf('\n===== 分阶段升压策略 =====\n');
adjust_time = [2000, 5000, 10000];
for ta = adjust_time
    Tl = Ton_150_stable;
    Th = Ton_150_stable * 1.3;
    for iter=1:15
        Tr = (Tl + Th)/2;
        [~, ~, P_end] = sim_full(Tr, 100, ta, dt_opt, ...
            V, C, AA, P_pump, rho160, get_rho, get_E);
        if P_end > 150
            Th = Tr;
        else
            Tl = Tr;
        end
    end
    Tr_opt = (Tl + Th)/2;
    fprintf('调整时间 %d ms：升压阶段 Ton=%.4f ms，稳定后 Ton=%.4f ms\n',...
        ta, Tr_opt, Ton_150_stable);
end

%% ====================== 7. 绘图 ======================
figure('Position',[100,100,900,400]);
plot(t_sample/1000, P_sample, 'b-', 'LineWidth',1.1);
hold on;
yline(100, 'r--', 'LineWidth',1.2, 'DisplayName', '目标100MPa');
xlabel('时间 (s)'); ylabel('压力 (MPa)');
title(sprintf('高压油管压力稳定曲线 (Ton=%.4f ms)', Ton_opt_100));
legend('show', 'Location', 'best');
grid on; ylim([95, 105]);

saveas(gcf, '第一问压力曲线.png');

disp('图1给出了压力波动的细节，result.xlsx已生成完毕。');
disp('命令行输出的"升压阶段Ton"即为2s/5s/10s快充时应设置的开启时长。');

%% ====================== 工具函数 ======================
    function val = lookup_linear(x_arr, y_arr, dx, x)
        if x < x_arr(1)
            val = y_arr(1);
        elseif x > x_arr(end)
            val = y_arr(end);
        else
            idx = floor((x - x_arr(1)) / dx) + 1;
            frac = (x - x_arr(idx)) / dx;
            val = y_arr(idx) + frac * (y_arr(idx+1) - y_arr(idx));
        end
    end

    function Q = q_inject(t)
        tm = mod(t, 100);
        if tm < 0.2
            Q = 100 * tm;
        elseif tm < 2.2
            Q = 20;
        elseif tm < 2.4
            Q = 240 - 100 * tm;
        else
            Q = 0;
        end
    end

    function flag = valve_state(t, Ton)
        tm = mod(t, Ton + T_close);
        flag = tm < Ton;
    end

    function dP = ode_p(t, P, Ton, V, C, AA, Pp, rho_p, f_rho, f_E)
        rho = f_rho(P);
        E = f_E(P);
        
        dm_in = 0;
        if valve_state(t, Ton)
            dP_delta = Pp - P;
            if dP_delta > 0
                q_vol = C * AA * sqrt(2 * dP_delta / rho_p);
                dm_in = rho_p * q_vol;
            end
        end
        
        q_out_vol = q_inject(t);
        dm_out = rho * q_out_vol;
        
        dP = (dm_in - dm_out) * E / (V * rho);
    end

    function [P_mean, loss_val] = sim_steady(Ton, P0, t_total, dt, t_start, V, C, AA, Pp, rho_p, f_rho, f_E)
        N = floor(t_total / dt) + 1;
        P = zeros(N, 1);
        P(1) = P0;
        
        for i = 1:N-1
            t = (i-1)*dt;
            Pi = P(i);
            
            k1 = ode_p(t, Pi, Ton, V, C, AA, Pp, rho_p, f_rho, f_E);
            k2 = ode_p(t+dt/2, Pi + dt/2*k1, Ton, V, C, AA, Pp, rho_p, f_rho, f_E);
            k3 = ode_p(t+dt/2, Pi + dt/2*k2, Ton, V, C, AA, Pp, rho_p, f_rho, f_E);
            k4 = ode_p(t+dt, Pi + dt*k3, Ton, V, C, AA, Pp, rho_p, f_rho, f_E);
            
            P(i+1) = Pi + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
        end
        
        idx_start = floor(t_start / dt) + 1;
        P_steady = P(idx_start:end);
        P_mean = mean(P_steady);
        loss_val = sum((P_steady - 100).^2) * dt;
    end

    function [t_arr, P_arr, P_end] = sim_full(Ton, P0, t_total, dt, V, C, AA, Pp, rho_p, f_rho, f_E)
        N = floor(t_total / dt) + 1;
        t_arr = (0:N-1)' * dt;
        P_arr = zeros(N, 1);
        P_arr(1) = P0;
        
        for i = 1:N-1
            t = t_arr(i);
            Pi = P_arr(i);
            
            k1 = ode_p(t, Pi, Ton, V, C, AA, Pp, rho_p, f_rho, f_E);
            k2 = ode_p(t+dt/2, Pi + dt/2*k1, Ton, V, C, AA, Pp, rho_p, f_rho, f_E);
            k3 = ode_p(t+dt/2, Pi + dt/2*k2, Ton, V, C, AA, Pp, rho_p, f_rho, f_E);
            k4 = ode_p(t+dt, Pi + dt*k3, Ton, V, C, AA, Pp, rho_p, f_rho, f_E);
            
            P_arr(i+1) = Pi + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
        end
        P_end = P_arr(end);
    end

end

%运行结果：
%  ===== 寻优单向阀开启时长 =====
% 迭代 1 | Ton=0.30000 ms | 稳态均值=108.5064 MPa
% 迭代 2 | Ton=0.25000 ms | 稳态均值=83.4413 MPa
% 迭代 3 | Ton=0.27500 ms | 稳态均值=99.9191 MPa
% 迭代 4 | Ton=0.28750 ms | 稳态均值=100.3629 MPa
% 迭代 5 | Ton=0.28125 ms | 稳态均值=97.5239 MPa
% 迭代 6 | Ton=0.28438 ms | 稳态均值=99.5983 MPa
% 迭代 7 | Ton=0.28594 ms | 稳态均值=99.3590 MPa
% 迭代 8 | Ton=0.28672 ms | 稳态均值=99.7830 MPa
% 迭代 9 | Ton=0.28711 ms | 稳态均值=99.7877 MPa
% 迭代10 | Ton=0.28730 ms | 稳态均值=99.9185 MPa
% 迭代11 | Ton=0.28740 ms | 稳态均值=99.8994 MPa
% 迭代12 | Ton=0.28745 ms | 稳态均值=99.9314 MPa
% 迭代13 | Ton=0.28748 ms | 稳态均值=99.9246 MPa
% 迭代14 | Ton=0.28749 ms | 稳态均值=99.9314 MPa
% 迭代15 | Ton=0.28749 ms | 稳态均值=99.9340 MPa
% 迭代16 | Ton=0.28750 ms | 稳态均值=99.9340 MPa
% 迭代17 | Ton=0.28750 ms | 稳态均值=99.9340 MPa
% 迭代18 | Ton=0.28750 ms | 稳态均值=99.9340 MPa
% 迭代19 | Ton=0.28750 ms | 稳态均值=99.9340 MPa
% 迭代20 | Ton=0.28750 ms | 稳态均值=99.9340 MPa
% 迭代21 | Ton=0.28750 ms | 稳态均值=99.9340 MPa
% 迭代22 | Ton=0.28750 ms | 稳态均值=99.9340 MPa
% 迭代23 | Ton=0.28750 ms | 稳态均值=99.9340 MPa
% 迭代24 | Ton=0.28750 ms | 稳态均值=99.9340 MPa
% 迭代25 | Ton=0.28750 ms | 稳态均值=99.9340 MPa
% 
% 【最优开启时长】Ton = 0.2875 ms
% 
% 正在生成20000ms压力曲线并保存至 result.xlsx...
% 已输出0~20000ms每10ms压力至 result.xlsx
% 
% ===== 求解150MPa稳态开启时长 =====
% 150MPa稳态Ton = 0.7525 ms
% 
% ===== 分阶段升压策略 =====
% 调整时间 2000 ms：升压阶段 Ton=0.8770 ms，稳定后 Ton=0.7525 ms
% 调整时间 5000 ms：升压阶段 Ton=0.7525 ms，稳定后 Ton=0.7525 ms
% 调整时间 10000 ms：升压阶段 Ton=0.7525 ms，稳定后 Ton=0.7525 ms
% 图1给出了压力波动的细节，result.xlsx已生成完毕。
% 命令行输出的"升压阶段Ton"即为2s/5s/10s快充时应设置的开启时长。