function [omega_1, omega_2, omega_3] = random_topsis(matrix)
    % 提取矩阵的各列
    N_t = matrix(:, 1);
    N_d = matrix(:, 2);
    N_s = matrix(:, 3);
    
    % Step 1: 标准化决策矩阵
    N_t_norm = N_t / norm(N_t);
    N_d_norm = N_d / norm(N_d);
    N_s_norm = N_s / norm(N_s);
    
    % 构建标准化决策矩阵
    norm_matrix = [N_t_norm, N_d_norm, N_s_norm];
    
    % Step 2: 确定理想解和负理想解
    ideal_solution = max(norm_matrix);
    negative_ideal_solution = min(norm_matrix);
    
    % Step 3: 随机生成权重
    num_random_weights = 10000; % 生成 10000 组随机权重
    weights = rand(num_random_weights, 3);
    weights = weights ./ sum(weights, 2); % 归一化权重使得权重之和为1
    
    best_score = -inf;
    best_weights = [0, 0, 0];
    
    % Step 4: 遍历每组随机权重，计算相对接近度并找出最优权重组合
    for i = 1:num_random_weights
        w = weights(i, :);
        
        % 构建加权标准化决策矩阵
        weighted_norm_matrix = norm_matrix .* w;
        
        % 计算与理想解和负理想解的距离
        d_pos = sqrt(sum((weighted_norm_matrix - ideal_solution) .^ 2, 2));
        d_neg = sqrt(sum((weighted_norm_matrix - negative_ideal_solution) .^ 2, 2));
        
        % 计算相对接近度
        relative_closeness = d_neg ./ (d_pos + d_neg);
        
        % 找到最大相对接近度对应的权重
        score = mean(relative_closeness);
        if score > best_score
            best_score = score;
            best_weights = w;
        end
    end
    
    % 返回最优权重
    omega_1 = best_weights(1);
    omega_2 = best_weights(2);
    omega_3 = best_weights(3);
end
