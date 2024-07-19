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
            case "+X"
                A_I = A_I_X;
            case "-X"
                A_I = A_I_X';
            case "+Y"
                A_I = A_I_Y;
            case "Y"
                A_I = A_I_Y;
            case "-Y"
                A_I = A_I_Y;
            case "+Z"
                A_I = A_I_Z;
            case "Z"
                A_I = A_I_Z;
            case "-Z"
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

    omega_1 = 0.2; % 权重因子
    omega_2 = 0.4; % 权重因子
    omega_3 = 0.3; % 权重因子
    if N_g > 0
        O = inf; % 如果几何干涉次数大于0，适应度设为无穷大
    else
        O = omega_1 * N_t + omega_2 * N_d + omega_3 * N_s;
    end
end