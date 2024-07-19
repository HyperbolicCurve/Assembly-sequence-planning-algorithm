function [best_solution, best_directions, best_fitness_value] = PSO_function(N, maxIter, w, c1, c2, vMax, minError, S, T, C, A_I_X, A_I_Y, A_I_Z)

    % 假设装配序列有21个零件
    m = 21;

    % 初始化装配方向矩阵 D
    directions = ["+X", "-X", "+Y", "-Y", "+Z", "-Z"];
    D = directions(randi(length(directions), N, m)); % 随机初始化每个粒子的装配方向
    
    valid_population_found = false;
    while ~valid_population_found

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
        if globalBestFitness ~= inf 
            valid_population_found = true;
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
    best_fitness_value = best_fitness_values;
    best_fitness = globalBestFitness;

    % 计算最优解的具体指标
    [~, N_t, N_d, N_s] = calculate_indicators(best_solution, S, T, C, best_directions, A_I_X, A_I_Y, A_I_Z);

    fprintf('最优装配序列: ');
    disp(best_solution);
    fprintf('最优适应度: %f\n', best_fitness);
    fprintf('装配工具的改变次数: %d\n', N_t);
    fprintf('装配方向的改变次数: %d\n', N_d);
    fprintf('装配过程中不稳定操作的次数: %d\n', N_s);
    fprintf('零件的安装方向: ');
    disp(best_directions);

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
