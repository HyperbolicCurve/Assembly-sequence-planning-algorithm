% 装配序列规划问题的灰狼优化算法（GWO）

% 输入参数
N = 500; % 种群数量
maxIter = 200; % 最大迭代次数
aMax = 2; % a 初始值
minError = 1e-2; % 最小误差

% 读取支撑矩阵 S
S = readmatrix('支撑矩阵.xlsx');

% 读取装配工具信息 T 并转换为字符串数组
T = readtable('装配工具信息.xlsx', 'ReadVariableNames', false);
T = table2array(T);  % 转换为数组
T = string(T);  % 转换为字符串数组

% 读取连接矩阵 C
C = readmatrix('连接矩阵.xlsx');

% 读取干涉矩阵 A_I
A_I_X = readmatrix('干涉矩阵-X.xlsx');
A_I_Y = readmatrix('干涉矩阵-Y.xlsx');
A_I_Z = readmatrix('干涉矩阵-Z.xlsx');

% 假设装配序列有21个零件
m = 21;
% 只允许使用+X和-X两个方向
directionOptions = ["X", "-X"];

% 初始化狼群位置
positions = zeros(N, m);
directions = strings(N, m);
for i = 1:N
    positions(i, 1) = 9;  % 将9号零件固定为第一个零件
    positions(i, 2:end) = setdiff(randperm(m), 9);  % 确保位置是有效的排列且不包括9号零件
    directions(i, 1) = "-X";  % 将9号零件的装配方向固定为-X
    directions(i, 2:end) = directionOptions(randi(length(directionOptions), 1, m-1));  % 随机选择其他零件的装配方向
end

% 初始化 Alpha, Beta, Delta 的位置及其适应度
alpha_pos = zeros(1, m);
alpha_dir = strings(1, m);
alpha_score = inf;
beta_pos = zeros(1, m);
beta_dir = strings(1, m);
beta_score = inf;
delta_pos = zeros(1, m);
delta_dir = strings(1, m);
delta_score = inf;

% 计算初始适应值
for i = 1:N
    fitness = fitness_function(positions(i, :), directions(i, :), S, T, C, A_I_X);
    if fitness < alpha_score
        alpha_score = fitness;
        alpha_pos = positions(i, :);
        alpha_dir = directions(i, :);
    elseif fitness < beta_score
        beta_score = fitness;
        beta_pos = positions(i, :);
        beta_dir = directions(i, :);
    elseif fitness < delta_score
        delta_score = fitness;
        delta_pos = positions(i, :);
        delta_dir = directions(i, :);
    end
end

% 迭代优化
best_fitness_values = zeros(maxIter, 1);
for iter = 1:maxIter
    a = aMax - iter * (aMax / maxIter); % 线性减少 a
    
    for i = 1:N
        % 创建一个临时位置和方向数组来存储更新后的值
        temp_positions = positions(i, :);
        temp_directions = directions(i, :);
        for j = 2:m  % 从第二个零件开始更新位置和方向
            r1 = rand();
            r2 = rand();
            
            A1 = 2 * a * r1 - a;
            C1 = 2 * r2;
            
            D_alpha = abs(C1 * alpha_pos(j) - positions(i, j));
            X1 = alpha_pos(j) - A1 * D_alpha;
            
            r1 = rand();
            r2 = rand();
            
            A2 = 2 * a * r1 - a;
            C2 = 2 * r2;
            
            D_beta = abs(C2 * beta_pos(j) - positions(i, j));
            X2 = beta_pos(j) - A2 * D_beta;
            
            r1 = rand();
            r2 = rand();
            
            A3 = 2 * a * r1 - a;
            C3 = 2 * r2;
            
            D_delta = abs(C3 * delta_pos(j) - positions(i, j));
            X3 = delta_pos(j) - A3 * D_delta;
            
            temp_positions(j) = round((X1 + X2 + X3) / 3); % 四舍五入以获得整数位置
            
            % 随机更新方向
            temp_directions(j) = directionOptions(randi(length(directionOptions)));
        end

        % 确保位置在有效范围内
        temp_positions = max(min(temp_positions, m), 1);
        
        % 确保每个零件只出现一次
        temp_positions = fix_positions(temp_positions);
        
        % 更新主位置数组中的位置和方向
        positions(i, :) = temp_positions;
        directions(i, :) = temp_directions;
        
        % 计算适应值
        fitness = fitness_function(positions(i, :), directions(i, :), S, T, C, A_I_X);
        
        % 更新 Alpha, Beta, Delta
        if fitness < alpha_score
            alpha_score = fitness;
            alpha_pos = positions(i, :);
            alpha_dir = directions(i, :);
        elseif fitness < beta_score
            beta_score = fitness;
            beta_pos = positions(i, :);
            beta_dir = directions(i, :);
        elseif fitness < delta_score
            delta_score = fitness;
            delta_pos = positions(i, :);
            delta_dir = directions(i, :);
        end
    end
    
    % 保存当前迭代的最佳适应度值
    best_fitness_values(iter) = alpha_score;
    
    % 判断是否达到最小误差停止条件
    if alpha_score < minError
        break;
    end
    
    % 尝试减少方向数量
    for i = 1:N
        unique_dirs = unique(directions(i, :));
        if length(unique_dirs) > 2
            % 使用 tabulate 函数统计每个方向的出现次数
            dir_counts = tabulate(directions(i, :));
            % 获取出现次数最多的两个方向
            [~, idx] = sort([dir_counts{:, 2}], 'descend');
            main_dirs = dir_counts(idx(1:2), 1);
            main_dirs = [main_dirs{:}];
            
            % 将其他方向替换为出现次数最多的两个方向之一
            for j = 2:m  % 从第二个零件开始替换方向
                if ~ismember(directions(i, j), main_dirs)
                    directions(i, j) = main_dirs(randi(2));
                end
            end
        end
    end
end

% 输出最优解
best_solution = alpha_pos;
best_directions = alpha_dir;
best_fitness = alpha_score;

% 计算最优解的具体指标
[N_g, N_t, N_d, N_s] = calculate_indicators(best_solution, best_directions, S, T, C, A_I_X);

fprintf('最优装配序列: ');
disp(best_solution);
fprintf('最优适应度: %f\n', best_fitness);
fprintf('装配工具的改变次数: %d\n', N_t);
fprintf('装配方向的改变次数: %d\n', N_d);
fprintf('装配过程中不稳定操作的次数: %d\n', N_s);
fprintf('零件的安装方向: ');
disp(best_directions);

% 绘制目标函数值和迭代次数的关系图
figure;
plot(1:iter, best_fitness_values, '-o');
xlabel('迭代次数');
ylabel('目标函数值');
title('目标函数值与迭代次数的关系');

% 适应度函数
function O = fitness_function(Z_h, D, S, T, C, A_I_X)
    n = length(Z_h);
    N_g = 0;
    N_t = 0;
    N_d = 0;
    N_s = 0;

    % 计算几何干涉次数 N_g
    for i = 1:n-1
        switch D(i)  % 使用装配序列中的索引来访问装配方向
            case 'X'
                A_I = A_I_X;
            case '-X'
                A_I = A_I_X';
            otherwise
                error('未知的装配方向');
        end

        if A_I(Z_h(i), Z_h(i+1)) == 1
            N_g = N_g + 1;
        end
    end

    % 计算工具改变次数 N_t
    for i = 1:n-1
        if T(Z_h(i)) ~= T(Z_h(i+1))
            N_t = N_t + 1;
        end
    end

    % 计算方向改变次数 N_d
    for i = 1:n-1
        if D(i) ~= D(i+1)
            N_d = N_d + 1;
        end
    end

    % 计算不稳定操作次数 N_s
    for i = 2:n
        if C(Z_h(i), Z_h(i-1)) ~= 1 && S(Z_h(i), Z_h(i-1)) ~= 1
            N_s = N_s + 1;
        end
    end

    % 设定权重系数
    omega_1 = 0.3;
    omega_2 = 0.3;
    omega_3 = 0.4;

    % 计算评价函数 E
    if N_g > 0
        O = inf; % 存在几何干涉时，适应度设为无穷大
    else
        O = omega_1 * N_t + omega_2 * N_d + omega_3 * N_s;
    end
end

function [N_g, N_t, N_d, N_s] = calculate_indicators(Z_h, D, S, T, C, A_I_X)
    n = length(Z_h);
    N_g = 0;
    N_t = 0;
    N_d = 0;
    N_s = 0;

    % 计算几何干涉次数 N_g
    for i = 1:n-1
        switch D(i)  % 确保装配方向的索引是合法的
            case 'X'
                A_I = A_I_X;
            case '-X'
                A_I = A_I_X';
            otherwise
                error('未知的装配方向');
        end

        if A_I(Z_h(i), Z_h(i+1)) == 1
            N_g = N_g + 1;
        end
    end

    % 计算工具改变次数 N_t
    for i = 1:n-1
        if T(Z_h(i)) ~= T(Z_h(i+1))
            N_t = N_t + 1;
        end
    end

    % 计算方向改变次数 N_d
    for i = 1:n-1
        if D(i) ~= D(i+1)
            N_d = N_d + 1;
        end
    end

    % 计算不稳定操作次数 N_s
    for i = 2:n
        if C(Z_h(i), Z_h(i-1)) ~= 1 && S(Z_h(i), Z_h(i-1)) ~= 1
            N_s = N_s + 1;
        end
    end
end

function fixedPosition = fix_positions(position)
    m = length(position);
    [~, ia, ~] = unique(position, 'stable');
    duplicate_indices = setdiff(1:m, ia);
    missing_elements = setdiff(1:m, position);

    % 替换重复项
    for i = 1:length(duplicate_indices)
        position(duplicate_indices(i)) = missing_elements(i);
    end
    fixedPosition = position;
end