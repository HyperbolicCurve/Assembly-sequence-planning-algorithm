function O = CalculateN(pso_best_solution, ga_bestSequence, gwo_best_solution, opt_gwo_best_solution, S, T, C, PSO_D, GA_D, GWO_D, Optimized_D)
    direction_arrays = {PSO_D, GA_D, GWO_D, Optimized_D};
    position_arrays = {pso_best_solution, ga_bestSequence, gwo_best_solution, opt_gwo_best_solution};
    num_directions = length(direction_arrays); 
    n = length(pso_best_solution);
    O = zeros(num_directions , 3); % 初始化输出矩阵，行数为 num_directions * n，列数为 3

    for k = 1:num_directions
        current_D = direction_arrays{k}; % 获取当前方向数组
        current_P = position_arrays{k};

        N_t = 0;
        N_d = 0;
        N_s = 0;

        % 计算工具改变次数 N_t
        for i = 1:n-1
            if T(current_P(i)) ~= T(current_P(i+1))
                N_t = N_t + 1;
            end
        end

        % 计算方向改变次数 N_d
        for i = 1:n-1
            if current_D(i) ~= current_D(i + 1)
                N_d = N_d + 1;
            end
        end

        % 计算不稳定操作次数 N_s
        for i = 2:n
            if C(current_P(i), current_P(i-1)) ~= 1 && S(current_P(i), current_P(i-1)) ~= 1
                N_s = N_s + 1;
            end
        end
        % 将结果存储在输出矩阵中
        O(k, :) = [N_t, N_d, N_s];
    end
end
