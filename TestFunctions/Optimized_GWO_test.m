% 运行参数设置
N = 30; % 种群大小
maxIter = 500; % 最大迭代次数
aMax = 1.6; % a的初始值
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
        [~, bestFitness, fitness_iter] = Optimized_GWO(N, maxIter, aMax, minError, testFunction, range, dim);
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

% 优化灰狼算法主函数
function [best_solution, best_fitness, fitness_iter] = Optimized_GWO(N, maxIter, aMax, minError, testFunction, range, dim)
    lb = range(1); % 下界
    ub = range(2); % 上界

    % 初始化种群
    positions = lb + (ub - lb) * rand(N, dim);
    fitness = arrayfun(@(i) feval(testFunction, positions(i, :)), 1:N);

    % 初始化 Alpha, Beta, Delta
    [sorted_fitness, idx] = sort(fitness);
    alpha_pos = positions(idx(1), :);
    alpha_score = sorted_fitness(1);
    beta_pos = positions(idx(2), :);
    beta_score = sorted_fitness(2);
    delta_pos = positions(idx(3), :);
    delta_score = sorted_fitness(3);

    % 迭代优化
    fitness_iter = zeros(maxIter, 1);
    for iter = 1:maxIter
        a = aMax - iter * (aMax / maxIter); % 线性减少 a
        for i = 1:N
            for j = 1:dim
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

                positions(i, j) = (X1 + X2 + X3) / 3;
            end

            % Lévy飞行操作
            if rand < 0.1
                positions(i, :) = levy_flight(positions(i, :), lb, ub);
            end

            % 确保位置在有效范围内
            positions(i, :) = max(min(positions(i, :), ub), lb);

            % 计算适应值
            fitness(i) = feval(testFunction, positions(i, :));

            % 更新 Alpha, Beta, Delta
            if fitness(i) < alpha_score
                alpha_score = fitness(i);
                alpha_pos = positions(i, :);
            elseif fitness(i) < beta_score
                beta_score = fitness(i);
                beta_pos = positions(i, :);
            elseif fitness(i) < delta_score
                delta_score = fitness(i);
                delta_pos = positions(i, :);
            end
        end

        % 保存当前迭代的最佳适应度值
        fitness_iter(iter) = alpha_score;

        % 判断是否达到最小误差停止条件
        if alpha_score < minError
            break;
        end
    end

    % 输出最优解
    best_solution = alpha_pos;
    best_fitness = alpha_score;
end

% Lévy飞行函数
function new_position = levy_flight(position, lb, ub)
    beta = 1.5;
    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / (gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)))^(1 / beta);
    u = randn(size(position)) * sigma;
    v = randn(size(position));
    step = u ./ abs(v).^(1 / beta);

    step_size = 0.01 * step .* (position - lb);
    new_position = position + step_size .* randn(size(position));

    % 确保位置在有效范围内
    new_position = max(min(new_position, ub), lb);
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
