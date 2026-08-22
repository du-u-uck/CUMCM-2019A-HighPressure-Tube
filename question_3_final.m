function question_3_final
clear; clc; close all;
format longG;

%% ====================== 1. 全局参数定义 ======================
% 高压油管参数
D_pipe = 10;    L_pipe = 500;
V_pipe = pi * (D_pipe/2)^2 * L_pipe;

% 柱塞参数
d_plunger = 5;
A_plunger = pi * (d_plunger/2)^2;
V0 = 20;        % 残余容积 mm³
P_low = 0.5;    % 低压供油压力 MPa

% 供油小孔参数
d_hole = 1.4;
A_hole = pi * (d_hole/2)^2;
C = 0.85;

% 针阀参数
d_needle = 2.5;
alpha = 9 * pi / 180;   % 半角转弧度
d_jet = 1.4;
A_jet = pi * (d_jet/2)^2;
injector_num = 2;       % 喷油嘴数量：2个

% 减压阀参数
d_relief = 1.4;
A_relief = pi * (d_relief/2)^2;
P_relief_set = 101;     % 减压阀开启压力 MPa
P_relief_low = 0.5;     % 减压阀出口低压 MPa

% 仿真参数
dt = 0.01;               % 仿真步长 ms
t_total = 10000;         % 总仿真时长 ms
t_steady_start = 5000;   % 稳态统计起始时刻 ms

%% ====================== 2. 预处理所有附件数据 ======================
% --- 附件3：E-P、rho-P 快速查表 ---
data_E = readmatrix("附件3-弹性模量与压力.xlsx");
P_raw_E = data_E(:,1);
E_raw = data_E(:,2);
P_step = 0.01;
P_table = (0:P_step:200)';
N_table = length(P_table);
E_table = interp1(P_raw_E, E_raw, P_table, 'linear', 'extrap');

rho_table = zeros(N_table, 1);
idx_100 = round(100 / P_step) + 1;
rho_table(idx_100) = 0.850;
for i = idx_100+1 : N_table
    E_mid = (E_table(i-1) + E_table(i)) / 2;
    drho = rho_table(i-1) / E_mid * P_step;
    rho_table(i) = rho_table(i-1) + drho;
end
for i = idx_100-1 : -1 : 1
    E_mid = (E_table(i) + E_table(i+1)) / 2;
    drho = rho_table(i+1) / E_mid * P_step;
    rho_table(i) = rho_table(i+1) - drho;
end

get_E = @(P) lookup_linear(P_table, E_table, P_step, P);
get_rho = @(P) lookup_linear(P_table, rho_table, P_step, P);

% --- 附件1：凸轮极径曲线 ---
data_cam = readmatrix("附件1-凸轮边缘曲线.xlsx");
theta_raw = data_cam(:,1);
r_raw = data_cam(:,2);
r_min = min(r_raw);
r_max = max(r_raw);
x_max = r_max - r_min;

theta_step = 0.001;
theta_table = (0:theta_step:2*pi)';
r_table = interp1(theta_raw, r_raw, theta_table, 'linear', 'extrap');
dr_dtheta_table = gradient(r_table, theta_step);

get_r = @(theta) lookup_linear(theta_table, r_table, theta_step, mod(theta, 2*pi));
get_drdtheta = @(theta) lookup_linear(theta_table, dr_dtheta_table, theta_step, mod(theta, 2*pi));

% --- 附件2：针阀升程曲线 ---
raw_needle = readmatrix("附件2-针阀运动曲线.xlsx");
rise_col = raw_needle(:,1:2);
rise_col(any(isnan(rise_col),2), :) = [];
t_rise = rise_col(:,1);
h_rise = rise_col(:,2);

fall_col = raw_needle(:,4:5);
fall_col(any(isnan(fall_col),2), :) = [];
t_fall = fall_col(:,1);
h_fall = fall_col(:,2);

t_needle_raw = [t_rise; 2.0; t_fall; 100.0];
h_needle_raw = [h_rise; 2.0; h_fall; 0.0];

t_needle_step = 0.01;
t_needle_table = (0:t_needle_step:100)';
h_table = interp1(t_needle_raw, h_needle_raw, t_needle_table, 'linear', 0);

get_h = @(t) lookup_linear(t_needle_table, h_table, t_needle_step, mod(t, 100));

%% ====================== 3. 工况1：双喷油 无减压阀 寻优 ======================
fprintf('===== 【工况1】双喷油嘴 无减压阀 寻优 =====\n');
omega_low = 0.050;
omega_high = 0.065;
iter_max = 20;

for iter = 1:iter_max
    omega_mid = (omega_low + omega_high) / 2;
    P_mean = sim_steady(omega_mid, t_total, dt, t_steady_start, injector_num, 0);
    fprintf('迭代%2d | ω=%.6f rad/ms | 稳态均值=%.4f MPa\n', iter, omega_mid, P_mean);
    
    if P_mean > 100
        omega_high = omega_mid;
    else
        omega_low = omega_mid;
    end
end
omega_opt_no_relief = (omega_low + omega_high) / 2;
fprintf('无减压阀最优角速度：ω = %.6f rad/ms = %.1f rpm\n\n', ...
    omega_opt_no_relief, omega_opt_no_relief*1000*60/(2*pi));

%% ====================== 4. 工况2：双喷油+减压阀 补偿寻优 ======================
fprintf('===== 【工况2】双喷油嘴+减压阀 补偿寻优 =====\n');
omega_low_r = omega_opt_no_relief * 1.0;  % 从无阀最优值开始搜
omega_high_r = omega_opt_no_relief * 1.05;

for iter = 1:iter_max
    omega_mid = (omega_low_r + omega_high_r) / 2;
    P_mean = sim_steady(omega_mid, t_total, dt, t_steady_start, injector_num, 1);
    fprintf('迭代%2d | ω=%.6f rad/ms | 稳态均值=%.4f MPa\n', iter, omega_mid, P_mean);
    
    if P_mean > 100
        omega_high_r = omega_mid;
    else
        omega_low_r = omega_mid;
    end
end
omega_opt_with_relief = (omega_low_r + omega_high_r) / 2;
fprintf('带减压阀最优角速度：ω = %.6f rad/ms = %.1f rpm\n\n', ...
    omega_opt_with_relief, omega_opt_with_relief*1000*60/(2*pi));

%% ====================== 5. 对比仿真与绘图 ======================
fprintf('正在生成对比曲线...\n');
[t1, P_pipe1, ~] = sim_full(omega_opt_no_relief, t_total, dt, injector_num, 0);
[t2, P_pipe2, ~] = sim_full(omega_opt_with_relief, t_total, dt, injector_num, 1);

figure('Position',[100,100,900,600]);
subplot(2,1,1);
plot(t1/1000, P_pipe1, 'b-', 'LineWidth',1);
yline(100, 'r--', '目标100MPa');
xlabel('时间 / s'); ylabel('油管压力 / MPa');
title('双喷油嘴 无减压阀 压力曲线');
grid on; ylim([95, 105]);

subplot(2,1,2);
plot(t2/1000, P_pipe2, 'g-', 'LineWidth',1);
yline(100, 'r--', '目标100MPa');
yline(P_relief_set, 'k--', '减压阀开启压力');
xlabel('时间 / s'); ylabel('油管 / MPa');
title('双喷油嘴+减压阀 压力曲线（角速度补偿后）');
grid on; ylim([95, 105]);

saveas(gcf, '第三问压力曲线.png');

%% ====================== 6. 波动指标统计 ======================
idx_steady = find(t1 >= t_steady_start, 1);
P1_steady = P_pipe1(idx_steady:end);
P2_steady = P_pipe2(idx_steady:end);

fprintf('\n===== 稳态压力指标对比 =====\n');
fprintf('无减压阀：均值=%.4f MPa，峰峰值=%.2f MPa，标准差=%.3f MPa\n', ...
    mean(P1_steady), max(P1_steady)-min(P1_steady), std(P1_steady));
fprintf('有减压阀：均值=%.4f MPa，峰峰值=%.2f MPa，标准差=%.3f MPa\n', ...
    mean(P2_steady), max(P2_steady)-min(P2_steady), std(P2_steady));
fprintf('减压阀使波动峰峰值降低：%.1f%%\n', ...
    100*( (max(P1_steady)-min(P1_steady)) - (max(P2_steady)-min(P2_steady)) ) / (max(P1_steady)-min(P1_steady)) );

%% ====================== 7. 核心子函数 ======================
    function val = lookup_linear(x_arr, y_arr, dx, x)
        if x < x_arr(1)
            val = y_arr(1);
        elseif x > x_arr(end)
            val = y_arr(end);
        else
            idx = floor((x - x_arr(1)) / dx) + 1;
            idx = min(idx, length(y_arr)-1);
            frac = (x - x_arr(idx)) / dx;
            val = y_arr(idx) + frac * (y_arr(idx+1) - y_arr(idx));
        end
    end

    function dY = ode_system(t, Y, omega, n_inj, relief_flag)
        Pp = Y(1);
        Ppipe = Y(2);
        
        theta = omega * t;
        dr_dtheta = get_drdtheta(theta);
        dV_dt = -A_plunger * dr_dtheta * omega;
        Vp = V0 + A_plunger * (x_max - (get_r(theta) - r_min));
        
        rho_p = get_rho(Pp);
        E_p = get_E(Pp);
        rho_pipe = get_rho(Ppipe);
        
        % 柱塞腔供油
        dm_supply = 0;
        if dV_dt < 0
            if Pp > Ppipe
                dP_delta = Pp - Ppipe;
                dm_supply = C * A_hole * sqrt(2 * rho_p * dP_delta);
            end
            dPp_dt = E_p / Vp * ( -dm_supply / rho_p - dV_dt );
        else
            if Pp <= P_low
                dPp_dt = 0;
            else
                dPp_dt = E_p / Vp * ( -dV_dt );
            end
        end
        
        % 双喷油总流量
        h = get_h(t);
        if h <= 0
            dm_inj_single = 0;
        else
            An = pi * d_needle * h * sin(alpha);
            A_inj = min(An, A_jet);
            dm_inj_single = C * A_inj * sqrt(2 * rho_pipe * Ppipe);
        end
        dm_inj_total = n_inj * dm_inj_single;
        
        % 减压阀流量
        dm_relief = 0;
        if relief_flag && Ppipe > P_relief_set
            dP_rel = Ppipe - P_relief_low;
            dm_relief = C * A_relief * sqrt(2 * rho_pipe * dP_rel);
        end
        
        % 油管压力微分
        E_pipe = get_E(Ppipe);
        dPpipe_dt = E_pipe / (V_pipe * rho_pipe) * (dm_supply - dm_inj_total - dm_relief);
        
        dY = [dPp_dt, dPpipe_dt];
    end

    function P_mean = sim_steady(omega, t_total, dt, t_start, n_inj, relief_flag)
        N = floor(t_total / dt) + 1;
        Y = zeros(N, 2);
        Y(1,:) = [0.5, 100];
        
        for i = 1:N-1
            t = (i-1)*dt;
            k1 = ode_system(t, Y(i,:), omega, n_inj, relief_flag);
            k2 = ode_system(t+dt/2, Y(i,:)+dt/2*k1, omega, n_inj, relief_flag);
            k3 = ode_system(t+dt/2, Y(i,:)+dt/2*k2, omega, n_inj, relief_flag);
            k4 = ode_system(t+dt, Y(i,:)+dt*k3, omega, n_inj, relief_flag);
            Y(i+1,:) = Y(i,:) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
            
            Y(i+1,1) = max(Y(i+1,1), 0.1);
            Y(i+1,2) = max(Y(i+1,2), 0.1);
        end
        
        idx_start = floor(t_start / dt) + 1;
        P_mean = mean(Y(idx_start:end, 2));
    end

    function [t_arr, P_pipe, P_plunger] = sim_full(omega, t_total, dt, n_inj, relief_flag)
        N = floor(t_total / dt) + 1;
        t_arr = (0:N-1)' * dt;
        Y = zeros(N, 2);
        Y(1,:) = [0.5, 100];
        
        for i = 1:N-1
            t = t_arr(i);
            k1 = ode_system(t, Y(i,:), omega, n_inj, relief_flag);
            k2 = ode_system(t+dt/2, Y(i,:)+dt/2*k1, omega, n_inj, relief_flag);
            k3 = ode_system(t+dt/2, Y(i,:)+dt/2*k2, omega, n_inj, relief_flag);
            k4 = ode_system(t+dt, Y(i,:)+dt*k3, omega, n_inj, relief_flag);
            Y(i+1,:) = Y(i,:) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
            
            Y(i+1,1) = max(Y(i+1,1), 0.1);
            Y(i+1,2) = max(Y(i+1,2), 0.1);
        end
        
        P_plunger = Y(:,1);
        P_pipe = Y(:,2);
    end
end