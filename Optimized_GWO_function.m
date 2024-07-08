function [best_solution, best_directions, bestValues] = Optimized_GWO_function(S, T, C, A_I_X, A_I_Y, A_I_Z, N, maxIter, aMax, minError)
    % 输入参数
    m = 21; % 假设装配序列有21个零件
    directions = ["X", "-X", "Y", "-Y", "Z", "-Z"];

    % 初始化狼群位置
    positions = zeros(N, m);
    directions_idx = randi(length(directions), N, m); % 装配方向索引

    for i = 1:N
        positions(i, :) = randperm(m);
        % 确保9号零件为第一个零件，并且方向为-X
        positions(i, :) = [9, setdiff(positions(i, :), 9, 'stable')];
        directions_idx(i, :) = [2, randi(length(directions), 1, m-1)];
    end

    % 初始化 Alpha, Beta, Delta 的位置及其适应度
    alpha_pos = zeros(1, m);
    alpha_dir = zeros(1, m);
    alpha_score = inf;
    beta_pos = zeros(1, m);
    beta_dir = zeros(1, m);
    beta_score = inf;
    delta_pos = zeros(1, m);
    delta_dir = zeros(1, m);
    delta_score = inf;

    % 计算初始适应值
    for i = 1:N
        fitness = fitness_function(positions(i, :), directions(directions_idx(i, :)), S, T, C, A_I_X, A_I_Y, A_I_Z);
        if fitness < alpha_score
            alpha_score = fitness;
            alpha_pos = positions(i, :);
            alpha_dir = directions_idx(i, :);
        elseif fitness < beta_score
            beta_score = fitness;
            beta_pos = positions(i, :);
            beta_dir = directions_idx(i, :);
        elseif fitness < delta_score
            delta_score = fitness;
            delta_pos = positions(i, :);
            delta_dir = directions_idx(i, :);
        end
    end

    % 迭代优化
    best_fitness_values = zeros(maxIter, 1);
    for iter = 1:maxIter
        a = aMax - iter * (aMax / maxIter); % 线性减少 a

        if iter > maxIter / 2
            directions = ["X", "-X"]; % 后半段限制方向选项
        end
        
        for i = 1:N
            % 创建一个临时位置数组来存储更新后的值
            temp_positions = positions(i, :);
            temp_dirs = directions_idx(i, :);
            for j = 1:m
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
            end
            
            % 确保位置在有效范围内
            temp_positions = max(min(temp_positions, m), 1);
            
            % 确保每个零件只出现一次
            temp_positions = fix_positions(temp_positions);
            
            % 使用Lévy飞行动作生成新位置
            [temp_positions, temp_dirs] = levy_flight(temp_positions, temp_dirs, directions, S, T, C, A_I_X, A_I_Y, A_I_Z);
            
            % 确保9号零件为第一个零件，并且方向为-X
            temp_positions = [9, setdiff(temp_positions, 9, 'stable')];
            temp_dirs(1) = 2; % 2对应于'-X'
            
            % 更新主位置数组中的位置
            positions(i, :) = temp_positions;
            directions_idx(i, :) = temp_dirs;
            
            % 计算适应值
            fitness = fitness_function(positions(i, :), directions(directions_idx(i, :)), S, T, C, A_I_X, A_I_Y, A_I_Z);
            
            % 更新 Alpha, Beta, Delta
            if fitness < alpha_score
                alpha_score = fitness;
                alpha_pos = positions(i, :);
                alpha_dir = directions_idx(i, :);
            elseif fitness < beta_score
                beta_score = fitness;
                beta_pos = positions(i, :);
                beta_dir = directions_idx(i, :);
            elseif fitness < delta_score
                delta_score = fitness;
                delta_pos = positions(i, :);
                delta_dir = directions_idx(i, :);
            end
        end
        
        % 保存当前迭代的最佳适应度值
        best_fitness_values(iter) = alpha_score;
        
        % 判断是否达到最小误差停止条件
        if alpha_score < minError
            break;
        end
    end

    % 输出最优解
    best_solution = alpha_pos;
    best_directions = alpha_dir;
    best_fitness = alpha_score;
    bestValues = best_fitness_values;
    
    % 确保best_directions中的值是有效的索引
    best_directions = max(min(best_directions, length(directions)), 1);

    % 计算最优解的具体指标
    [N_g, N_t, N_d, N_s] = calculate_indicators(best_solution, directions(best_directions), S, T, C, A_I_X, A_I_Y, A_I_Z);
    
    fprintf('最优装配序列: ');
    disp(best_solution);
    fprintf('最优适应度: %f\n', best_fitness);
    fprintf('装配工具的改变次数: %d\n', N_t);
    fprintf('装配方向的改变次数: %d\n', N_d);
    fprintf('装配过程中不稳定操作的次数: %d\n', N_s);
    fprintf('零件的安装方向: ');
    disp(directions(best_directions));
    
end

%%% 辅助函数 %%%
function [new_position, new_directions] = levy_flight(position, directions, all_directions, S, T, C, A_I_X, A_I_Y, A_I_Z)
    % 选择随机一个Global-Swap算子
    if rand < 0.33
        new_position = global_swap(position);
    % 选择随机一个Cut-and-Paste算子
    elseif rand < 0.66
        new_position = cut_and_paste(position);
    % 选择随机一个2-Opt算子
    else
        new_position = two_opt(position);
    end

    % 确保位置在有效范围内
    new_position = max(min(new_position, length(position)), 1);
    
    % 确保每个零件只出现一次
    new_position = fix_positions(new_position);

    % 随机改变装配方向
    new_directions = directions;
    for i = 1:length(directions)
        if rand < 0.1  % 10%的概率改变方向
            new_directions(i) = randi(length(all_directions));
        end
    end

    % 确保new_directions中的索引在有效范围内
    new_directions = max(min(new_directions, length(all_directions)), 1);

    % 如果新位置的目标函数值更好，则接受新位置
    if fitness_function(new_position, all_directions(new_directions), S, T, C, A_I_X, A_I_Y, A_I_Z) < ...
       fitness_function(position, all_directions(directions), S, T, C, A_I_X, A_I_Y, A_I_Z)
        position = new_position;
        directions = new_directions;
    end
    % 返回新的位置和方向
    new_position = position;
    new_directions = directions;
end

% Global-Swap算子
function new_position = global_swap(position)
    idx = randperm(length(position), 2);
    new_position = position;
    new_position(idx(1)) = position(idx(2));
    new_position(idx(2)) = position(idx(1));
end

% Cut-and-Paste算子
function new_position = cut_and_paste(position)
    idx = sort(randperm(length(position), 2));
    new_position = [position(1:idx(1)-1), position(idx(2):end), position(idx(1):idx(2)-1)];
end

% 2-Opt算子
function new_position = two_opt(position)
    idx = sort(randperm(length(position), 2));
    new_position = position;
    new_position(idx(1):idx(2)) = position(idx(2):-1:idx(1));
end

% 适应度函数
function O = fitness_function(Z_h, D, S, T, C, A_I_X, A_I_Y, A_I_Z)
    n = length(Z_h);
    N_g = 0;
    N_t = 0;
    N_d = 0;
    N_s = 0;

    % 计算几何干涉次数 N_g
    for i = 1:n-1
        switch D(i)  % 使用装配序列中的方向
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
    omega_3 = 0.3;

    % 计算评价函数 E
    if N_g > 0
        O = inf; % 存在几何干涉时，适应度设为无穷大
    else
        O = omega_1 * N_t + omega_2 * N_d + omega_3 * N_s;
    end
end

% 计算指标的函数
function [N_g, N_t, N_d, N_s] = calculate_indicators(Z_h, D, S, T, C, A_I_X, A_I_Y, A_I_Z)
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
        if D(i) ~= D(i+1)
            N_d = N_d + 1;
        end
    end

    % 计算不稳定操作次数 N_s
    for i = 2:n
        if C(Z_h(i), Z_h(i-1)) ~= 1 && S(Z_h(i-1), Z_h(i)) ~= 1
            N_s = N_s + 1;
        end
    end
end

% 修复位置的函数，确保每个零件只出现一次
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
