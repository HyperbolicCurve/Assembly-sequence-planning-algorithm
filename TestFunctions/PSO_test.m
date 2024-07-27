% 运行参数设置
N = 30; % 种群大小
maxIter = 500; % 最大迭代次数
aMax = 2; % a的初始值
minError = 1e-6; % 最小误差停止条件
numRuns = 30; % 每个测试函数运行次数

% 测试函数列表
testFunctions = {@F1, @F2, @F3, @F4, @F5, @F6}; 
ranges = {
    [-100, 100],   % F1
    [-100, 100],   % F2
    [-1.28, 1.28], % F3
    [-600, 600],   % F4
    [-32, 32],     % F5
    [-500, 500]    % F6
};
dims = [30, 30, 30, 30, 4, 30]; % 维度设置
results = struct();

% 运行每个测试函数
for funcIdx = 1:length(testFunctions)
    testFunction = testFunctions{funcIdx};
    range = ranges{funcIdx};
    dim = dims(funcIdx);
    fitnessValues = zeros(numRuns, 1);
    runTimes = zeros(numRuns, 1);

    for run = 1:numRuns
        tic;
         [~, bestFitness] = PSO_function(N, maxIter, minError, testFunction, range, dim);
        runTimes(run) = toc;
        fitnessValues(run) = bestFitness;
    end

    % 计算统计值
    bestValue = min(fitnessValues);
    worstValue = max(fitnessValues);
    meanValue = mean(fitnessValues);
    stdValue = std(fitnessValues);
    meanRunTime = mean(runTimes);

    % 存储结果
    results(funcIdx).function = func2str(testFunction);
    results(funcIdx).bestValue = bestValue;
    results(funcIdx).worstValue = worstValue;
    results(funcIdx).meanValue = meanValue;
    results(funcIdx).stdValue = stdValue;
    results(funcIdx).meanRunTime = meanRunTime;
end

% 输出结果
for funcIdx = 1:length(results)
    fprintf('Function: %s\n', results(funcIdx).function);
    fprintf('Best Value: %f\n', results(funcIdx).bestValue);
    fprintf('Worst Value: %f\n', results(funcIdx).worstValue);
    fprintf('Mean Value: %f\n', results(funcIdx).meanValue);
    fprintf('Standard Deviation: %f\n', results(funcIdx).stdValue);
    fprintf('Mean Run Time: %f seconds\n', results(funcIdx).meanRunTime);
    fprintf('\n');
end

function [best_solution, best_fitness_value] = PSO_function(N, maxIter, minError, testFunction, range, dim)
    % PSO参数：
    w = 0.9; % 惯性因子
    c1 = 0.2; % 个体学习因子
    c2 = 0.2; % 社会学习因子
    vMax = 4; % 最大速度限制
    lb = range(1); % 下界
    ub = range(2); % 上界

    % 初始化种群
    positions = lb + (ub - lb) * rand(N, dim);
    velocities = zeros(N, dim); % Initialize velocities as zero
    
    % Initialize local and global best
    localBestPositions = positions;
    localBestFitness = inf * ones(N, 1);
    globalBestPosition = positions(1, :); % Initialize as first particle's position
    globalBestFitness = inf;

    % Calculate initial fitness values
    for i = 1:N
        fitness = testFunction(positions(i, :));
        localBestFitness(i) = fitness;
        if fitness < globalBestFitness
            globalBestFitness = fitness;
            globalBestPosition = positions(i, :);
        end
    end

    % Iterative optimization
    best_fitness_values = zeros(maxIter, 1);
    for iter = 1:maxIter
        for i = 1:N
            % Update velocity
            r1 = rand(1, dim);
            r2 = rand(1, dim);
            velocities(i, :) = w * velocities(i, :) ...
                + c1 * r1 .* (localBestPositions(i, :) - positions(i, :)) ...
                + c2 * r2 .* (globalBestPosition - positions(i, :));

            % Limit velocity
            velocities(i, :) = max(min(velocities(i, :), vMax), -vMax);

            % Update position
            positions(i, :) = positions(i, :) + velocities(i, :);

            % Ensure positions are within bounds
            positions(i, :) = max(min(positions(i, :), 1), 0);

            % Calculate fitness
            fitness = testFunction(positions(i, :));

            % Update local best
            if fitness < localBestFitness(i)
                localBestFitness(i) = fitness;
                localBestPositions(i, :) = positions(i, :);
            end

            % Update global best
            if fitness < globalBestFitness
                globalBestFitness = fitness;
                globalBestPosition = positions(i, :);
            end
        end

        % Save best fitness value for the current iteration
        best_fitness_values(iter) = globalBestFitness;

        % Check termination condition
        if globalBestFitness < minError
            break;
        end
    end

    % Output best solution and details
    best_solution = globalBestPosition;
    best_fitness_value = best_fitness_values(end); % Assuming last fitness value is the best

end

% 测试函数示例
function o = F1(x)
    o = max(abs(x));
end

function o = F2(x)
    dim = size(x, 2);
    o = 0;
    for i = 1:dim
        o = o + sum(x(1:i))^2;
    end
end

function o = F3(x)
    dim = size(x, 2);
    o = sum((1:dim) .* (x.^4)) + rand;
end

function o = F4(x)
    dim = size(x, 2);
    o = sum(x.^2) / 4000 - prod(cos(x ./ sqrt(1:dim))) + 1;
end

function o = F5(x)
    dim = size(x, 2);
    o = -20 * exp(-0.2 * sqrt(sum(x.^2) / dim)) - exp(sum(cos(2 * pi .* x)) / dim) + 20 + exp(1);
end

function o = F6(x)
    o = sum(-x .* sin(sqrt(abs(x))));
end

