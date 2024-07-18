% 装配序列规划问题的粒子群算法（PSO）

% 输入参数
N = 200; % 种群数量
maxIter = 200; % 最大迭代次数
w = 0.9; % 惯性因子
c1 = 0.6; % 个体学习因子
c2 = 0.6; % 社会学习因子
vMax = 4; % 最大速度限制
minError = 1e-6; % 算法终止的最小误差

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

% 初始化装配方向矩阵 D
directions = ["+X", "-X", "+Y", "-Y", "+Z", "-Z"];
D = directions(randi(length(directions), N, m)); % 随机初始化每个粒子的装配方向

% 初始化粒子群的位置和速度
positions = zeros(N, m);
velocities = zeros(N, m);
for i = 1:N
    positions(i, :) = randperm(m);
    velocities(i, :) = randn(1, m);
end

% 初始化局部最优位置和全局最优位置
localBestPositions = positions;
localBestDirections = D;
localBestFitness = inf * ones(N, 1);
globalBestPosition = positions(1, :); % 初始化为第一个粒子的初始位置
globalBestDirections = D(1, :); % 初始化为第一个粒子的初始方向
globalBestFitness = inf;

% 计算初始适应值
for i = 1:N
    fitness = fitness_function(positions(i, :), S, T, C, D(i, :), A_I_X, A_I_Y, A_I_Z);
    localBestFitness(i) = fitness;
    if fitness < globalBestFitness
        globalBestFitness = fitness;
        globalBestPosition = positions(i, :);
        globalBestDirections = D(i, :);
    end
end

% 如果初始化的适应度全为无穷大，则生成一个可行解
if isempty(globalBestPosition)
    for i = 1:N
        positions(i, :) = randperm(m);
        D(i, :) = directions(randi(length(directions), 1, m)); % 随机初始化装配方向
        fitness = fitness_function(positions(i, :), S, T, C, D(i, :), A_I_X, A_I_Y, A_I_Z);
        if fitness < globalBestFitness
            globalBestFitness = fitness;
            globalBestPosition = positions(i, :);
            globalBestDirections = D(i, :);
            break;
        end
    end
end

% 迭代优化
best_fitness_values = zeros(maxIter, 1);
for iter = 1:maxIter
    for i = 1:N
        % 更新速度
        r1 = rand(1, m);
        r2 = rand(1, m);
        velocities(i, :) = w * velocities(i, :) ...
            + c1 * r1 .* (localBestPositions(i, :) - positions(i, :)) ...
            + c2 * r2 .* (globalBestPosition - positions(i, :));
        
        % 限制速度
        velocities(i, :) = max(min(velocities(i, :), vMax), -vMax);
        
        % 更新位置
        positions(i, :) = round(positions(i, :) + velocities(i, :)); % 确保位置为整数
        
        % 确保位置在有效范围内
        positions(i, :) = max(min(positions(i, :), m), 1);
        
        % 确保每个零件只出现一次
        positions(i, :) = fix_positions(positions(i, :), m);
        
        % 随机改变部分方向
        D(i, randi(m, 1, floor(m/2))) = directions(randi(length(directions), 1, floor(m/2)));
        
        % 计算适应值
        fitness = fitness_function(positions(i, :), S, T, C, D(i, :), A_I_X, A_I_Y, A_I_Z);
        
        % 更新局部最优
        if fitness < localBestFitness(i)
            localBestFitness(i) = fitness;
            localBestPositions(i, :) = positions(i, :);
            localBestDirections(i, :) = D(i, :);
        end
        
        % 更新全局最优
        if fitness < globalBestFitness
            globalBestFitness = fitness;
            globalBestPosition = positions(i, :);
            globalBestDirections = D(i, :);
        end
    end
    
    % 保存当前迭代的最佳适应度值
    best_fitness_values(iter) = globalBestFitness;
    
    % 判断是否达到最小误差停止条件
    if globalBestFitness < minError
        break;
    end
end

% 输出最优解
best_solution = globalBestPosition;
best_directions = globalBestDirections;
best_fitness = globalBestFitness;

% 计算最优解的具体指标
[N_g, N_t, N_d, N_s] = calculate_indicators(best_solution, S, T, C, best_directions, A_I_X, A_I_Y, A_I_Z);

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



% 计算指标的函数
function [N_g, N_t, N_d, N_s] = calculate_indicators(Z_h, S, T, C, D, A_I_X, A_I_Y, A_I_Z)
    n = length(Z_h);
    N_g = 0;
    N_t = 0;
    N_d = 0;
    N_s = 0;

    % 计算几何干涉次数 N_g
    for i = 1:n-1
        switch D(i)  % 使用装配方向
            case "+X"
                A_I = A_I_X;
            case "-X"
                A_I = A_I_X';
            case "+Y"
                A_I = A_I_Y;
            case "-Y"
                A_I = A_I_Y';
            case "+Z"
                A_I = A_I_Z;
            case "-Z"
                A_I = A_I_Z';
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

% 确保每个零件只出现一次的函数
function unique_positions = fix_positions(positions, m)
    % 确保每个零件只出现一次
    [~, idx] = unique(positions, 'stable');
    if length(idx) < m
        missing_elements = setdiff(1:m, positions);
        positions(setdiff(1:m, idx)) = missing_elements;
    end
    unique_positions = positions;
end
