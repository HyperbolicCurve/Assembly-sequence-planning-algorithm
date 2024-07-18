function [best_solution, best_directions, fitness_iter] = GWO_function(N, maxIter, aMax, minError, S, T, C, A_I_X, A_I_Y, A_I_Z)
    % 装配序列规划问题的灰狼优化算法（GWO）
    
    % 假设装配序列有21个零件
    m = 21;
    % 装配方向矩阵 D
    directions_full = ["+X", "-X", "+Y", "-Y", "+Z", "-Z"];
    directions_restricted = ["+X", "-X"];
    directions = directions_full;
    
    % 初始化 Alpha, Beta, Delta 的位置及其适应度
    alpha_pos = zeros(1, m);
    alpha_dir = strings(1, m);
    alpha_score = inf;
    beta_pos = zeros(1, m);
    beta_score = inf;
    delta_pos = zeros(1, m);
    delta_score = inf;

    % 确保初始化的种群有一个非无穷大的适应度
    valid_population_found = false;
    while ~valid_population_found
        positions = zeros(N, m);
        directions_idx = randi(length(directions), N, m); % 装配方向索引
    
        for i = 1:N
            positions(i, :) = randperm(m);
            % 确保9号零件为第一个零件，并且方向为-X
            positions(i, :) = [9, setdiff(positions(i, :), 9, 'stable')];
            directions_idx(i, :) = [2, randi(length(directions), 1, m-1)];
        end
    
        % 计算初始适应值并检查是否有有效的适应度
        for i = 1:N
            fitness = fitness_function(positions(i, :), S, T, C, ...
                                       directions(directions_idx(i, :)), A_I_X, A_I_Y, A_I_Z);
            if fitness < alpha_score
                alpha_pos = positions(i, :);
                alpha_dir = directions_idx(i, :);
                alpha_score = fitness;
            elseif fitness < beta_score
                beta_pos = positions(i, :);
                beta_score = fitness;
            elseif fitness < delta_score
                delta_pos = positions(i, :);
                delta_score = fitness;
            end
            % 如果Alpha的适应度已经有效，则不需要更差的解，可以停止搜索
            if alpha_score ~= inf
                valid_population_found = true;
            end
        end
    end
    
    % 迭代优化
    best_fitness_values = zeros(maxIter, 1);
    for iter = 1:maxIter
        % 动态调整可选择的方向
        if iter > maxIter / 2
            directions = directions_restricted;
            % 修正当前方向索引，确保不会超过新的方向集范围
            directions_idx(directions_idx > length(directions_restricted)) = randi(length(directions_restricted), size(directions_idx(directions_idx > length(directions_restricted))));
        else
            directions = directions_full;
        end
        
        a = aMax - iter * (aMax / maxIter); % 线性减少 a
        step_size = 1 / (1 + exp(-10 * (iter / maxIter - 0.5))); % 动态调整步长
        for i = 1:N
            % 创建一个临时位置和方向数组来存储更新后的值
            temp_positions = positions(i, :);
            temp_directions = directions_idx(i, :);
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
                
                temp_positions(j) = round((X1 + X2 + X3) / 3 * step_size); % 乘以步长进行调整
                % 随机更新方向
                temp_directions(j) = randi(length(directions));
            end
    
            % 确保位置在有效范围内
            temp_positions = max(min(temp_positions, m), 1);
            
            % 确保每个零件只出现一次
            temp_positions = fix_positions(temp_positions);
            
            % 更新主位置数组中的位置和方向
            positions(i, :) = temp_positions;
            directions_idx(i, :) = temp_directions;
            
            % 计算适应值
            fitness = fitness_function(positions(i, :), S, T, C, directions(directions_idx(i, :)), A_I_X, A_I_Y, A_I_Z);
            
            % 更新 Alpha, Beta, Delta
            if fitness < alpha_score
                alpha_score = fitness;
                alpha_pos = positions(i, :);
                alpha_dir = directions_idx(i, :);
            elseif fitness < beta_score
                beta_score = fitness;
                beta_pos = positions(i, :);
            elseif fitness < delta_score
                delta_score = fitness;
                delta_pos = positions(i, :);
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
    fitness_iter = best_fitness_values;

    % 计算最优解的具体指标
    [N_g, N_t, N_d, N_s] = calculate_indicators(best_solution, S, T, C, directions(best_directions), A_I_X, A_I_Y, A_I_Z);
    
    fprintf('最优装配序列: ');
    disp(best_solution);
    fprintf('最优适应度: %f\n', best_fitness);
    fprintf('装配工具的改变次数: %d\n', N_t);
    fprintf('装配方向的改变次数: %d\n', N_d);
    fprintf('装配过程中不稳定操作的次数: %d\n', N_s);
    fprintf('零件的安装方向: ');
    disp(directions(best_directions));

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
