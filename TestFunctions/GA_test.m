% 运行参数设置
populationSize = 30;
maxGenerations = 500;
crossoverRate = 0.8;
mutationRate = 0.1;
minError = 1e-6;
numRuns = 30;

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
        bestFitnessValues = GA_function(populationSize, maxGenerations, crossoverRate, mutationRate, minError, testFunction, range, dim);
        runTimes(run) = toc;
        fitnessValues(run) = bestFitnessValues(end);
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


function bestFitnessValues = GA_function(populationSize, maxGenerations, crossoverRate, mutationRate, minError, testFunction, range, dim)
    lb = range(1); % 下界
    ub = range(2); % 上界

    % 初始化种群
    population = lb + (ub - lb) * rand(populationSize, dim);
    
    % 计算初始适应值
    fitnessValues = zeros(populationSize, 1);
    for i = 1:populationSize
        fitnessValues(i) = feval(testFunction, population(i, :));
    end

    % 初始化最佳适应值
    [bestFitness, ~] = min(fitnessValues);

    % 进化过程
    bestFitnessValues = zeros(maxGenerations, 1);
    for generation = 1:maxGenerations
        % 选择
        selectedPopulation = selection(population, fitnessValues);

        % 交叉
        newPopulation = crossover(selectedPopulation, crossoverRate);

        % 变异
        newPopulation = mutation(newPopulation, mutationRate);

        % 计算新种群的适应值
        fitnessValues = zeros(populationSize, 1);
        for i = 1:populationSize
            fitnessValues(i) = feval(testFunction, newPopulation(i, :));
        end

        % 更新种群
        population = newPopulation;

        % 更新最佳适应值
        [currentBestFitness, ~] = min(fitnessValues);
        if currentBestFitness < bestFitness
            bestFitness = currentBestFitness;
        end

        % 保存当前代的最佳适应值
        bestFitnessValues(generation) = bestFitness;

        % 判断是否达到最小误差停止条件
        if bestFitness < minError
            break;
        end
    end
end

% 选择操作（锦标赛选择）
function selectedPopulation = selection(population, fitnessValues)
    % 使用锦标赛选择
    populationSize = size(population, 1);
    selectedPopulation = zeros(size(population));
    for i = 1:populationSize
        % 随机选择两个个体
        idx1 = randi(populationSize);
        idx2 = randi(populationSize);
        % 选择适应值更小的个体
        if fitnessValues(idx1) < fitnessValues(idx2)
            selectedPopulation(i, :) = population(idx1, :);
        else
            selectedPopulation(i, :) = population(idx2, :);
        end
    end
end

% 交叉操作（单点交叉）
function newPopulation = crossover(population, crossoverRate)
    populationSize = size(population, 1);
    numGenes = size(population, 2);
    newPopulation = population;
    for i = 1:2:populationSize-1  % 确保处理成对的个体
        if rand < crossoverRate
            parent1 = population(i, :);
            parent2 = population(i+1, :);
            cp = randi(numGenes-1);  % 随机选择一个交叉点
            newPopulation(i, :) = [parent1(1:cp), parent2(cp+1:end)];
            newPopulation(i+1, :) = [parent2(1:cp), parent1(cp+1:end)];
        end
    end
end

% 变异操作
function mutatedPopulation = mutation(population, mutationRate)
    populationSize = size(population, 1);
    numGenes = size(population, 2);
    mutatedPopulation = population;
    for i = 1:populationSize
        for j = 1:numGenes
            if rand < mutationRate
                mutatedPopulation(i, j) = rand * 2 - 1; % 随机变异
            end
        end
    end
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


