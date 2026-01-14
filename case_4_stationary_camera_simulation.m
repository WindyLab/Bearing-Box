clc; clear; close all

% Trajectory Generation Parameters
delta_t = 0.05;
t = linspace(0, 20, 20 / delta_t + 1);
N = length(t);

% Generate reference trajectory (3D spiral)
r_o = 5;
theta = 0:pi/100:4*pi;
x_box = r_o*cos(theta);
y = r_o*sin(theta);
z = zeros(1, N);
p_o = [x_box; y; z];

% Compute velocity and acceleration numerically
v_o = zeros(3, N);
a_o = zeros(3, N);
for i = 1:N
    if i == N
        v_o(:, i) = (p_o(:, 2) - p_o(:, i)) / delta_t;
        continue
    end
    v_o(:, i) = (p_o(:, i+1) - p_o(:, i)) / delta_t;
end

for i = 1:N
    if i == N
        a_o(:, i) = (v_o(:, 2) - v_o(:, i)) / delta_t;
        continue
    end
    a_o(:, i) = (v_o(:, i+1) - v_o(:, i)) / delta_t;
end

% Static camera setup
v_c = [0; 0; 0];
p_c_init = [0; 0; -5];
p_c = p_c_init + v_c * t;

% System parameters
alpha = 0.7;
alphas = ones(1, N) * alpha;
e_3g = [0; 0; 9.8]; % Gravity vector

% Compute measurement vectors
n = zeros(3, N);
for i = 1:N
    temp_n = a_o(:, i) - e_3g;
    temp_n = temp_n / norm(temp_n);
    n(:, i) = temp_n;
end

T_bar = (p_o - p_c) / alpha;

g = zeros(3, N);
for i = 1:N
    temp_g = p_o(:, i) - p_c(:, i);
    temp_g = temp_g / norm(temp_g);
    g(:, i) = temp_g;
end

% Monte Carlo simulation parameters
times = 100;
nees_bo_100 = zeros(times, N);
nees_bb_100 = zeros(times, N);
pred_x_bo_100 = zeros(times, 9, N);
pred_x_bb_100 = zeros(times, 10, N);

% Main simulation loop
for jj = 1:times
    % Add measurement noise
    n_std = 0.01;
    T_std = 0.1;
    g_std = 0.01;
    var_n = n_std^2;
    var_T = T_std^2;
    var_g = g_std^2;
    
    n_noise = zeros(3, N);
    T_bar_noise = zeros(3, N);
    g_noise = zeros(3, N);
    
    for i = 1:N
        n_noise_temp = n(:, i) + [n_std * randn(1, 1); n_std * randn(1, 1); n_std * randn(1, 1)];
        T_noise_temp = T_bar(:, i) + [T_std * randn(1, 1); T_std * randn(1, 1); T_std * randn(1, 1)];
        g_noise_temp = g(:, i) + [g_std * randn(1, 1); g_std * randn(1, 1); g_std * randn(1, 1)];
        
        n_noise(:, i) = n_noise_temp / norm(n_noise_temp);
        T_bar_noise(:, i) = T_noise_temp;
        g_noise(:, i) = g_noise_temp / norm(g_noise_temp);
    end
    
    % Kalman filter initialization
    initial_scale = 1;
    initial_pos = p_c(:, 1) + (p_o(:, 1) - p_c(:, 1)) * 1.5;
    x_init = [initial_pos(1); initial_pos(2); initial_pos(3); 1; 0; 0; 1; 0; 0; 1.5*alpha];
    P_init = eye(10) * 0.1;
    
    Q = eye(10) * 0.05;
    Q(1, 1) = 0.0001^2; % Position noise
    Q(2, 2) = 0.0001^2;
    Q(3, 3) = 0.0001^2;
    
    Q(4, 4) = 0.001^2; % Velocity noise
    Q(5, 5) = 0.001^2;
    Q(6, 6) = 0.001^2;
    
    Q(10, 10) = 0.0001^2; % Alpha noise
    
    A = [eye(3), delta_t*eye(3), 0.5*delta_t^2*eye(3), zeros(3, 1);
         zeros(3), eye(3), delta_t*eye(3), zeros(3, 1);
         zeros(3), zeros(3), eye(3), zeros(3, 1);
         zeros(1, 9), 1];
    
    x_box = x_init;
    P_box = P_init;
    
    % Bearing-Box estimation
    pred_x_bb = zeros(10, N);
    e_nees_bb = zeros(1, N);
    
    for i = 1:N
        % Measurement at time i
        temp_n = n_noise(:, i);
        temp_T_bar = T_bar_noise(:, i);
        P_n = eye(3) - temp_n * temp_n';
        
        % Prediction step
        x_box = A * x_box;
        P_box = A * P_box * A' + Q;
        
        % Measurement noise covariance
        temp_f = x_box(7:9) - e_3g;
        temp_f_value = norm(temp_f);
        V_k = [-x_box(10) * eye(3), zeros(3);
               zeros(3), temp_f_value * P_n];
        Sigma_k = [var_T * eye(3), zeros(3);
                   zeros(3), var_n * eye(3)];
        R_k = V_k * Sigma_k * V_k';
        
        % Measurement model
        H_k = [eye(3), zeros(3), zeros(3), -temp_T_bar;
               zeros(3), zeros(3), P_n, zeros(3, 1)];
        z_k = [p_c(:, i); P_n * e_3g];
        
        % True state for validation
        true_x = [p_o(:, i); v_o(:, i); a_o(:, i); alpha];
        
        % Kalman update
        K_k = P_box * H_k' * pinv(H_k * P_box * H_k' + R_k);
        x_box = x_box + K_k * (z_k - H_k * x_box);
        P_box = (eye(10) - K_k * H_k) * P_box;
        
        % NEES calculation
        e_nees = (true_x - x_box)' * pinv(P_box) * (true_x - x_box) / 10;
        
        % Store results
        pred_x_bb(:, i) = x_box;
        e_nees_bb(i) = e_nees;
    end
    
    % Bearing-Only estimation
    x_only = x_init(1:9);
    P_only = P_init(1:9, 1:9);
    Q_only = Q(1:9, 1:9);
    A_only = A(1:9, 1:9);
    
    pred_x_bo = zeros(9, N);
    e_nees_bo = zeros(1, N);
    
    for i = 1:N
        % Measurement at time i
        temp_g = g_noise(:, i);
        P_g = eye(3) - temp_g * temp_g';
        
        % Prediction step
        x_only = A_only * x_only;
        P_only = A_only * P_only * A_only' + Q_only;
        
        % Measurement noise covariance
        temp_r = norm(x_only(1:3) - p_c(:, i));
        V_k = temp_r * P_g;
        Sigma_k = var_g * eye(3);
        R_k = V_k * Sigma_k * V_k';
        
        % Measurement model
        H_k = [P_g, zeros(3), zeros(3)];
        z_k = P_g * p_c(:, i);
        
        % True state for validation
        true_x = [p_o(:, i); v_o(:, i); a_o(:, i)];
        
        % Kalman update
        K_k = P_only * H_k' * pinv(H_k * P_only * H_k' + R_k);
        x_only = x_only + K_k * (z_k - H_k * x_only);
        P_only = (eye(9) - K_k * H_k) * P_only;
        
        % Store results
        pred_x_bo(:, i) = x_only;
        e_nees = (true_x - x_only)' * pinv(P_only) * (true_x - x_only) / 9;
        e_nees_bo(i) = e_nees;
    end
    
    % Store Monte Carlo results
    nees_bb_100(jj, :) = e_nees_bb;
    nees_bo_100(jj, :) = e_nees_bo;
    pred_x_bo_100(jj, :, :) = pred_x_bo;
    pred_x_bb_100(jj, :, :) = pred_x_bb;
end

% Compute average results
avg_nees_bb = mean(nees_bb_100);
avg_nees_bo = mean(nees_bo_100);
avg_pred_bo = mean(pred_x_bo_100);
avg_pred_bb = mean(pred_x_bb_100);
avg_pred_bo = reshape(avg_pred_bo, 9, N);
avg_pred_bb = reshape(avg_pred_bb, 10, N);


% Calculate average prediction across MC runs
avg_pred_bo_test = zeros(9, N);
for i = 1:times
    temp_Matrix = reshape(pred_x_bo_100(i, :, :), 9, N);
    avg_pred_bo_test = avg_pred_bo_test + temp_Matrix;
end
avg_pred_bo_test = avg_pred_bo_test / times;


% Error analysis
dis_error_bb = zeros(1, N);
dis_error_bo = zeros(1, N);
avg_dis_error_bb = zeros(1, N);
avg_dis_error_bo = zeros(1, N);
relative_dis_error_bb = zeros(1, N);
relative_dis_error_bo = zeros(1, N);

for i = 1:N
    dis_error_bb(i) = norm(pred_x_bb(1:3, i) - p_o(:, i));
    dis_error_bo(i) = norm(pred_x_bo(1:3, i) - p_o(:, i));
    
    avg_dis_error_bb(i) = norm(avg_pred_bb(1:3, i) - p_o(:, i));
    avg_dis_error_bo(i) = norm(avg_pred_bo(1:3, i) - p_o(:, i));
    
    true_dis = norm(p_o(1:3, i) - p_c(:, i));
    relative_dis_error_bb(i) = dis_error_bb(i) / true_dis;
    relative_dis_error_bo(i) = dis_error_bo(i) / true_dis;
end

% Display average relative errors
fprintf('Average relative error (BO): %.4f\n', sum(relative_dis_error_bo) / length(relative_dis_error_bo));
fprintf('Average relative error (BB): %.4f\n', sum(relative_dis_error_bb) / length(relative_dis_error_bb));

% Color scheme for plots
x_color = [31, 119, 180]/255;
y_color = [44, 160, 44]/255;
z_color = [255, 127, 14]/255;
bo_color = [158, 24, 157]/255;
fontscale = 15;

% Main visualization figure
figure(1)
subplot(2, 3, [1 4])
h1 = plot3(p_c(1, :), p_c(2, :), p_c(3, :), 'Color', x_color, 'MarkerEdgeColor', 'auto', 'LineWidth', 2); hold on
h2 = plot3(p_o(1, :), p_o(2, :), p_o(3, :), 'Color', y_color, 'LineWidth', 2); hold on
h3 = plot3(pred_x_bb(1, :), pred_x_bb(2, :), pred_x_bb(3, :), 'Color', z_color, 'LineStyle', '-', 'LineWidth', 2); hold on
h4 = plot3(pred_x_bo(1, :), pred_x_bo(2, :), pred_x_bo(3, :), 'Color', bo_color, 'LineStyle', '-', 'LineWidth', 2); hold on

% Mark starting points
plot3(p_o(1, 1), p_o(2, 1), p_o(3, 1), 'Color', y_color, 'LineWidth', 2, 'Marker', 'diamond'); hold on
plot3(pred_x_bb(1, 1), pred_x_bb(2, 1), pred_x_bb(3, 1), 'Color', z_color, 'LineStyle', '-.', 'LineWidth', 2, 'Marker', 'pentagram'); hold on
plot3(pred_x_bo(1, 1), pred_x_bo(2, 1), pred_x_bo(3, 1), 'Color', bo_color, 'LineStyle', '-.', 'LineWidth', 2, 'Marker', 'pentagram'); hold on
plot3(p_c(1, 1), p_c(2, 1), p_c(3, 1), 'Color', x_color, 'MarkerEdgeColor', 'auto', 'LineWidth', 2, 'Marker', 'diamond'); hold on

l1 = legend([h1, h2, h3, h4], '$p_c$', '$p_o$', '$\hat{p}_o$ (Ours)', '$\hat{p}_o$ (BO)', 'Interpreter', 'latex');
set(l1, 'Fontname', 'Times New Roman', 'FontAngle', 'italic', 'FontSize', fontscale)
xlabel('$x(m)$', 'Fontname', 'Times New Roman', 'FontSize', fontscale, 'Interpreter', 'latex')
ylabel('$y(m)$', 'Fontname', 'Times New Roman', 'FontSize', fontscale, 'Interpreter', 'latex')
zlabel('$z(m)$', 'Fontname', 'Times New Roman', 'FontSize', fontscale, 'Interpreter', 'latex')
grid on
zlim([-5 5])
set(gca, 'Fontname', 'Times New Roman', 'FontSize', fontscale)

% Position error subplot
subplot(2, 3, 2)
h1 = plot(t, dis_error_bo, 'Color', bo_color, 'LineStyle', '-', 'LineWidth', 2); hold on
h2 = plot(t, dis_error_bb, 'Color', z_color, 'LineStyle', '-', 'LineWidth', 2);
xlabel('$t(s)$', 'Fontname', 'Times New Roman', 'FontSize', fontscale, 'Interpreter', 'latex')
ylabel('Position Error', 'Fontname', 'Times New Roman', 'FontSize', fontscale, 'Interpreter', 'latex')
l1 = legend([h1, h2], 'Bearing-only', 'Bearing-box', 'Interpreter', 'latex');
xlim([1 max(t)])
ylim([0 7])
set(gca, 'Fontname', 'Times New Roman', 'FontSize', fontscale)

% Velocity estimation subplot
subplot(2, 3, 5)
h1 = plot(t, pred_x_bo(4, :), 'Color', bo_color, 'LineStyle', '-', 'LineWidth', 2); hold on
h2 = plot(t, pred_x_bb(4, :), 'Color', z_color, 'LineStyle', '-', 'LineWidth', 2); hold on
h3 = plot(t, v_o(1, :), 'Color', x_color, 'LineStyle', '-', 'LineWidth', 2);
xlabel('$t(s)$', 'Fontname', 'Times New Roman', 'FontSize', fontscale, 'Interpreter', 'latex')
ylabel('$v_x$ Estimation', 'Fontname', 'Times New Roman', 'FontSize', fontscale, 'Interpreter', 'latex')
l1 = legend([h1, h2, h3], 'BO', 'Ours', '$v_x$', 'Interpreter', 'latex');
xlim([1 max(t)])
ylim([-4 10])
set(gca, 'Fontname', 'Times New Roman', 'FontSize', fontscale)

% Alpha estimation subplot
subplot(2, 3, 3)
h1 = plot(t, pred_x_bb(end, :), 'Color', z_color, 'LineStyle', '-', 'LineWidth', 2); hold on
h2 = plot(t, alphas, 'Color', x_color, 'LineStyle', '-', 'LineWidth', 2);
xlabel('$t(s)$', 'Fontname', 'Times New Roman', 'FontSize', fontscale, 'Interpreter', 'latex')
ylabel('$\alpha$', 'Fontname', 'Times New Roman', 'FontSize', fontscale, 'Interpreter', 'latex')
l1 = legend([h1, h2], '$\hat{\alpha}$', '$\alpha$', 'Interpreter', 'latex');
xlim([1 max(t)])
ylim([0 1.5])
set(gca, 'Fontname', 'Times New Roman', 'FontSize', fontscale)

% NEES subplot
subplot(2, 3, 6)
h1 = plot(t, e_nees_bo, 'Color', bo_color, 'LineStyle', '-', 'LineWidth', 2); hold on
h2 = plot(t, e_nees_bb, 'Color', z_color, 'LineStyle', '-', 'LineWidth', 2);
xlabel('$t(s)$', 'Fontname', 'Times New Roman', 'FontSize', fontscale, 'Interpreter', 'latex')
ylabel('NEES', 'Fontname', 'Times New Roman', 'FontSize', fontscale, 'Interpreter', 'latex')
l1 = legend([h1, h2], 'Bearing-only', 'Bearing-box', 'Interpreter', 'latex');
xlim([1 max(t)])
ylim([0.5 2])
set(gca, 'Fontname', 'Times New Roman', 'FontSize', fontscale)

% Prepare data for export
p_c = p_c';
p_o = p_o';
avg_pred_bo = avg_pred_bo';
avg_pred_bb = avg_pred_bb';
avg_dis_error_bb = avg_dis_error_bb';
avg_dis_error_bo = avg_dis_error_bo';
dis_error_bb = dis_error_bb';
dis_error_bo = dis_error_bo';
t = t';
alphas = alphas';
avg_nees_bo = avg_nees_bo';
avg_nees_bb = avg_nees_bb';
