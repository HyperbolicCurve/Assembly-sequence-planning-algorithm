% 装配序列规划问题的离散布谷鸟算法（DCA）

% 输入参数
N = 100; % 种群数量
Times = 200; % 最大迭代次数
Pa = 0.25; % 寄生概率
NGlobalSwap = 5; % Global-Swap操作次数
NSegments = 3; % Local-Swap分段数

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

% 装配方向矩阵 D
% 这里假设D已经被正确定义为每个零件的装配方向
% 请确保D的长度与装配序列的长度一致
D = ["X", "-X", "Y", "-Y", "Z", "-Z", "X", "-X", "Y", "-Y", "Z", "-Z", "X", "-X", "Y", "-Y", "Z", "-Z", "X", "-X", "Z"];

% 生成初始种群
pop = zeros(N, m);
for i = 1:N
    pop(i, :) = randperm(m);
end

% 初始化保存每次迭代的最佳适应度值
best_fitness_values = zeros(1, Times);

% 迭代优化
for t = 1:Times
    % Lévy飞行产生新解
    for i = 1:N
        new_solution = pop(i, :);
        
        % Global-Swap操作
        for k = 1:NGlobalSwap
            p1 = randi(m);
            p2 = randi(m);
            new_solution([p1 p2]) = new_solution([p2 p1]);
        end
        
        % Cut-and-Paste操作
        l = randi([1 3]);
        p = randi(m - l + 1);
        segment = new_solution(p:p+l-1);
        new_solution(p:p+l-1) = [];
        insert_pos = randi(m - l + 1);
        new_solution = [new_solution(1:insert_pos-1), segment, new_solution(insert_pos:end)];
        
        % 2-Opt操作
        p1 = randi(m);
        p2 = randi(m);
        if p1 > p2
            temp = p1;
            p1 = p2;
            p2 = temp;
        end
        new_solution(p1:p2) = flip(new_solution(p1:p2));
        
        % 确保新解没有重复的零件
        if length(unique(new_solution)) == m
            % 计算新解的适应度
            new_fitness = fitness_function(new_solution, S, T, C, D, A_I_X, A_I_Y, A_I_Z);
            
            % 如果新解更优，则替换
            if new_fitness < fitness_function(pop(i, :), S, T, C, D, A_I_X, A_I_Y, A_I_Z)
                pop(i, :) = new_solution;
            end
        end
    end
    
    % 寄生巢更新
    fitness_values = zeros(1, N);
    for i = 1:N
        fitness_values(i) = fitness_function(pop(i, :), S, T, C, D, A_I_X, A_I_Y, A_I_Z);
    end
    [~, idx] = sort(fitness_values, 'descend');
    num_to_replace = round(Pa * N);
    
    for i = 1:num_to_replace
        worst_idx = idx(i);
        new_solution = pop(worst_idx, :);
        
        % Local-Swap操作
        segment_len = round(m / NSegments);
        for j = 1:NSegments
            segment_start = (j-1) * segment_len + 1;
            segment_end = min(j * segment_len, m);
            if segment_end - segment_start >= 1
                p1 = segment_start + randi(segment_end - segment_start + 1) - 1;
                p2 = segment_start + randi(segment_end - segment_start + 1) - 1;
                new_solution([p1 p2]) = new_solution([p2 p1]);
            end
        end
        
        % Break-Join操作
        direction_list = unique(D); % 获取唯一的装配方向
        for dir = 1:length(direction_list)
            current_direction = direction_list(dir);
            % 找到所有当前方向的零件位置
            indices = find(D == current_direction);
            if length(indices) >= 2
                % 随机选择一个片段
                start_idx = randi(length(indices) - 1);
                end_idx = start_idx + 1;
                segment = new_solution(indices(start_idx):indices(end_idx));
                
                % 在剩余的零件序列中找到相同方向的零件
                remaining_indices = setdiff(1:m, indices(start_idx):indices(end_idx));
                remaining_positions = find(D(remaining_indices) == current_direction);
                if ~isempty(remaining_positions)
                    insert_idx = remaining_indices(randi(length(remaining_positions)));
                    
                    % 插入片段，确保总长度不变
                    new_solution = [new_solution(1:insert_idx), segment, new_solution(insert_idx+1:end)];
                    if length(new_solution) > m
                        new_solution = new_solution(1:m); % 修正长度
                    end
                end
            end
        end
        
        % 确保新解没有重复的零件
        if length(unique(new_solution)) == m
            % 计算新解的适应度
            new_fitness = fitness_function(new_solution, S, T, C, D, A_I_X, A_I_Y, A_I_Z);
            
            % 如果新解更优，则替换
            if new_fitness < fitness_function(pop(worst_idx, :), S, T, C, D, A_I_X, A_I_Y, A_I_Z)
                pop(worst_idx, :) = new_solution;
            end
        end
    end
    
    % 保存当前迭代的最佳适应度值
    best_fitness_values(t) = min(fitness_values);
end

% 输出最优解
best_solution = pop(1, :);
best_fitness = fitness_function(best_solution, S, T, C, D, A_I_X, A_I_Y, A_I_Z);
for i = 2:N
    current_fitness = fitness_function(pop(i, :), S, T, C, D, A_I_X, A_I_Y, A_I_Z);
    if current_fitness < best_fitness
        best_solution = pop(i, :);
        best_fitness = current_fitness;
    end
end

% 计算最优解的具体指标
[N_g, N_t, N_d, N_s] = calculate_indicators(best_solution, S, T, C, D, A_I_X, A_I_Y, A_I_Z);

fprintf('最优装配序列: ');
disp(best_solution);
fprintf('最优适应度: %f\n', best_fitness);
fprintf('装配工具的改变次数: %d\n', N_t);
fprintf('装配方向的改变次数: %d\n', N_d);
fprintf('装配过程中不稳定操作的次数: %d\n', N_s);
fprintf('零件的安装方向: ');
disp(D(best_solution));

% 绘制目标函数值和迭代次数的关系图
figure;
plot(1:Times, best_fitness_values, '-o');
xlabel('迭代次数');
ylabel('目标函数值');
title('目标函数值与迭代次数的关系');

% 适应度函数
function O = fitness_function(Z_h, S, T, C, D, A_I_X, A_I_Y, A_I_Z)
    n = length(Z_h);
    N_g = 0;
    N_t = 0;
    N_d = 0;
    N_s = 0;

    % 计算几何干涉次数 N_g
    for i = 1:n-1
        switch D(Z_h(i))  % 使用装配序列中的索引来访问装配方向
            case 'X'
                A_I = A_I_X;
            case '-X'
                A_I = A_I_X';
            case 'Y'
                A_I = A_I_Y;
            case '-Y'
                A_I = A_I_Y';
            case 'Z'
                A_I = A_I_Z;
            case '-Z'
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
        if D(Z_h(i)) ~= D(Z_h(i+1))
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

% 计算指标的函数
function [N_g, N_t, N_d, N_s] = calculate_indicators(Z_h, S, T, C, D, A_I_X, A_I_Y, A_I_Z)
    n = length(Z_h);
    N_g = 0;
    N_t = 0;
    N_d = 0;
    N_s = 0;

    % 计算几何干涉次数 N_g
    for i = 1:n-1
        switch D(Z_h(i))  % 使用装配序列中的索引来访问装配方向
            case 'X'
                A_I = A_I_X;
            case '-X'
                A_I = A_I_X';
            case 'Y'
                A_I = A_I_Y;
            case '-Y'
                A_I = A_I_Y';
            case 'Z'
                A_I = A_I_Z;
            case '-Z'
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
        if D(Z_h(i)) ~= D(Z_h(i+1))
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
