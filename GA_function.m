function [bestSequence, bestDirections, bestValues] = GA_function(populationSize, maxGenerations, crossoverRate, mutationRate, minError, S, T, C, A_I_X, A_I_Y, A_I_Z)
    % 定义方向矩阵
    directions = ["+X", "-X", "+Y", "-Y", "+Z", "-Z"];

    % 假设装配序列有21个零件
    m = 21;

    valid_population_found = false;
    while ~valid_population_found
        % 初始化种群
        population = repmat(struct('sequence', [], 'directions', []), populationSize, 1);
        for i = 1:populationSize
            population(i).sequence = randperm(m);
            population(i).directions = directions(randi(length(directions), 1, m));
        end
    
        % 计算初始适应值
        fitnessValues = zeros(populationSize, 1);
        for i = 1:populationSize
            fitnessValues(i) = fitness_function(population(i).sequence, S, T, C, population(i).directions, A_I_X, A_I_Y, A_I_Z);
            if fitnessValues(i) ~= inf
                valid_population_found = true;
            end
        end
    end

    % 初始化最佳适应值
    [bestFitness, bestIndex] = min(fitnessValues);
    bestSolution = population(bestIndex);

    % 进化过程
    bestFitnessValues = zeros(maxGenerations, 1);
    for generation = 1:maxGenerations
        % 选择
        selectedPopulation = selection(population, fitnessValues);

        % 交叉
        newPopulation = crossover(selectedPopulation, crossoverRate);

        % 变异
        newPopulation = mutation(newPopulation, mutationRate, directions);

        % 计算新种群的适应值
        fitnessValues = zeros(populationSize, 1);
        for i = 1:populationSize
            fitnessValues(i) = fitness_function(newPopulation(i).sequence, S, T, C, newPopulation(i).directions, A_I_X, A_I_Y, A_I_Z);
        end

        % 更新种群
        population = newPopulation;

        % 更新最佳适应值
        [currentBestFitness, bestIndex] = min(fitnessValues);
        if currentBestFitness < bestFitness
            bestFitness = currentBestFitness;
            bestSolution = population(bestIndex);
        end

        % 保存当前代的最佳适应值
        bestFitnessValues(generation) = bestFitness;

        % 判断是否达到最小误差停止条件
        if bestFitness < minError
            break;
        end
    end

    % 输出最优解
    bestSequence = bestSolution.sequence;
    bestDirections = bestSolution.directions;
    bestValues = bestFitnessValues;

    % 计算最优解的具体指标
    [N_g, N_t, N_d, N_s] = calculate_indicators(bestSequence, S, T, C, bestDirections, A_I_X, A_I_Y, A_I_Z);

    fprintf('最优装配序列: ');
    disp(bestSequence);
    fprintf('最优适应度: %f\n', bestFitness);
    fprintf('装配工具的改变次数: %d\n', N_t);
    fprintf('装配方向的改变次数: %d\n', N_d);
    fprintf('装配过程中不稳定操作的次数: %d\n', N_s);
    fprintf('零件的安装方向: ');
    disp(bestDirections);

end

% 适应度函数
function O = fitness_function(Z_h, S, T, C, D, A_I_X, A_I_Y, A_I_Z)
    n = length(Z_h);
    N_g = 0;
    N_t = 0;
    N_d = 0;
    N_s = 0;

    % 计算几何干涉次数 N_g
    for i = 1:n-1
        switch D(i)
            case '+X'
                A_I = A_I_X;
            case '-X'
                A_I = A_I_X';
            case '+Y'
                A_I = A_I_Y;
            case '-Y'
                A_I = A_I_Y';
            case '+Z'
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
    omega_3 = 0.4;

    % 计算评价函数 E
    if N_g > 0
        O = inf; % 存在几何干涉时，适应度设为无穷大
    else
        O = omega_1 * N_t + omega_2 * N_d + omega_3 * N_s;
    end
end

% 选择操作（锦标赛选择）
function selectedPopulation = selection(population, fitnessValues)
    % 使用锦标赛选择
    populationSize = length(population);
    selectedPopulation = repmat(struct('sequence', [], 'directions', []), populationSize, 1);
    for i = 1:populationSize
        % 随机选择两个个体
        idx1 = randi(populationSize);
        idx2 = randi(populationSize);
        % 选择适应值更小的个体
        if fitnessValues(idx1) < fitnessValues(idx2)
            selectedPopulation(i) = population(idx1);
        else
            selectedPopulation(i) = population(idx2);
        end
    end
end

% 交叉操作（部分映射交叉 PMX 和单点交叉）
function newPopulation = crossover(population, crossoverRate)
    populationSize = length(population);
    newPopulation = population;
    for i = 1:2:populationSize-1  % 确保处理成对的个体
        if rand < crossoverRate
            parent1 = population(i);
            parent2 = population(i+1);
            % 使用部分映射交叉 (PMX) 处理序列
            [child1Sequence, child2Sequence] = pmx(parent1.sequence, parent2.sequence);
            % 使用单点交叉处理方向
            [child1Directions, child2Directions] = single_point_crossover(parent1.directions, parent2.directions);
            newPopulation(i).sequence = fix_positions(child1Sequence);
            newPopulation(i).directions = child1Directions;
            newPopulation(i+1).sequence = fix_positions(child2Sequence);
            newPopulation(i+1).directions = child2Directions;
        end
    end
end

% PMX 交叉操作
function [child1, child2] = pmx(parent1, parent2)
    n = length(parent1);
    child1 = nan(1, n);
    child2 = nan(1, n);

    % 随机选择两个交叉点
    cp1 = randi(n-1);
    cp2 = randi([cp1+1, n]);

    % 复制交叉片段
    child1(cp1:cp2) = parent1(cp1:cp2);
    child2(cp1:cp2) = parent2(cp1:cp2);

    % 填充子代1的剩余部分
    for i = 1:n
        if ~ismember(parent2(i), child1)
            for j = 1:n
                if isnan(child1(j))
                    child1(j) = parent2(i);
                    break;
                end
            end
        end
    end

    % 填充子代2的剩余部分
    for i = 1:n
        if ~ismember(parent1(i), child2)
            for j = 1:n
                if isnan(child2(j))
                    child2(j) = parent1(i);
                    break;
                end
            end
        end
    end
end

% 单点交叉操作
function [child1, child2] = single_point_crossover(parent1, parent2)
    n = length(parent1);
    cp = randi(n-1);  % 随机选择一个交叉点

    % 生成子代
    child1 = [parent1(1:cp) parent2(cp+1:end)];
    child2 = [parent2(1:cp) parent1(cp+1:end)];
end

% 修正位置函数
function unique_positions = fix_positions(positions)
    m = length(positions);
    missing_elements = setdiff(1:m, positions);
    duplicate_indices = find(histc(positions, unique(positions)) > 1);

    % 替换重复的元素
    for i = 1:length(duplicate_indices)
        positions(duplicate_indices(i)) = missing_elements(i);
    end
    unique_positions = positions;
end

% 变异操作
function mutatedPopulation = mutation(population, mutationRate, directions)
    populationSize = length(population);
    n = length(population(1).sequence);
    mutatedPopulation = population;
    for i = 1:populationSize
        if rand < mutationRate
            % 使用交换变异
            mutationPoints = randperm(n, 2);
            mp1 = mutationPoints(1);
            mp2 = mutationPoints(2);
            temp = mutatedPopulation(i).sequence(mp1);
            mutatedPopulation(i).sequence(mp1) = mutatedPopulation(i).sequence(mp2);
            mutatedPopulation(i).sequence(mp2) = temp;
        end
        if rand < mutationRate
            % 随机改变一个方向
            mp = randi(n);
            mutatedPopulation(i).directions(mp) = directions(randi(length(directions)));
        end
        % 确保唯一性
        mutatedPopulation(i).sequence = fix_positions(mutatedPopulation(i).sequence);
    end
end
