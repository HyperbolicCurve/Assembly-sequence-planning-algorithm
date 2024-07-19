% 种群数量
N = 100; 
% 迭代次数
maxIter = 300;

% PSO参数：
w = 0.9; % 惯性因子
c1 = 0.2; % 个体学习因子
c2 = 0.2; % 社会学习因子
vMax = 4; % 最大速度限制
minError_1 = 1e-2; % 算法终止的最小误差

% GA参数
crossoverRate = 0.8; % 交叉概率
mutationRate = 0.1; % 变异概率
minError_2 = 1e-2; % 算法终止的最小误差

% GWO&Optimized_GWO参数
aMax = 1.6;
minError_3 = 1e-2;

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

% PSO
 fprintf("PSO algorithm:\n")
 [pso_best_solution, pso_best_directions, pso_fitness_values] = PSO_function(N, maxIter, w, c1, c2, vMax, minError_1, S, T, C, A_I_X, A_I_Y, A_I_Z);

% GA
 fprintf("GA algorithm:\n")
 [ga_bestSequence, ga_bestDirections, ga_fitness_values] = GA_function(N, maxIter, crossoverRate, mutationRate, minError_2, S, T, C, A_I_X, A_I_Y, A_I_Z);

% GWO
 fprintf("GWO algorithm:\n")
 [gwo_best_solution, gwo_best_directions, gwo_fitness_values] = GWO_function(N, maxIter, aMax, minError_3, S, T, C, A_I_X, A_I_Y, A_I_Z);

% Optimized_GWO
fprintf("Optimized GWO algorithm:\n")
[opt_gwo_best_solution, opt_gwo_best_directions, opt_gwo_fitness_values] = Optimized_GWO_function(S, T, C, A_I_X, A_I_Y, A_I_Z, N, maxIter, aMax, minError_3);

% Calculate N
N_Matrix = CalculateN(pso_best_solution, ga_bestSequence, gwo_best_solution, opt_gwo_best_solution,  S, T, C, pso_best_directions, ga_bestDirections, gwo_best_directions, opt_gwo_best_directions);

% 随机TOPSIS法，On
[omega_1, omega_2, omega_3] = random_topsis(N_Matrix);
disp(['omega_1: ', num2str(omega_1)]);
disp(['omega_2: ', num2str(omega_2)]);
disp(['omega_3: ', num2str(omega_3)]);

% 绘制适应度变化图 (GWO, PSO, GA)
figure;
plot(1:maxIter, gwo_fitness_values, 'r', 'LineWidth', 2);
hold on;
plot(1:maxIter, pso_fitness_values, 'g', 'LineWidth', 2);
plot(1:maxIter, ga_fitness_values, 'b', 'LineWidth', 2);
plot(1:maxIter, opt_gwo_fitness_values, 'm', 'LineWidth', 2)

xlabel('迭代次数');
ylabel('适应度值');
legend('GWO', 'PSO', 'GA', 'Optimized-GWO');
title('GWO, PSO, GA, Optimized-GWO适应度变化图');
hold off;


