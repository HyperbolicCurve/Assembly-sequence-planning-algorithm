% Clear the workspace and close all figures
clear;
close all;

% Function 1: F1(x) = max_i(|x_i|), 1 <= i <= n
n = 30;
[x1, x2] = meshgrid(-100:1:100, -100:1:100);
F1 = max(abs(x1), abs(x2));

figure;
mesh(x1, x2, F1);
title('Function F_1');
xlabel('x_1');
ylabel('x_2');
zlabel('F_1(x)');
colorbar;

% Function 2: F2(x) = sum_{i=1}^{n}((sum_{j=1}^{i} x_{j})^2)
F2 = zeros(size(x1));
for i = 1:n
    F2 = F2 + (x1 + x2).^2;
end

figure;
mesh(x1, x2, F2);
title('Function F_2');
xlabel('x_1');
ylabel('x_2');
zlabel('F_2(x)');
colorbar;

% Function 3: F3(x) = sum_{i=1}^{n} i * x_i^4 + random[0,1)
[x1, x2] = meshgrid(-1.28:0.01:1.28, -1.28:0.01:1.28);
F3 = x1.^4 + x2.^4 + rand(size(x1));

figure;
mesh(x1, x2, F3);
title('Function F_3');
xlabel('x_1');
ylabel('x_2');
zlabel('F_3(x)');
colorbar;

% Function 4: F4(x) = 1 + sum_{i=1}^{n} x_i^2 - prod(cos(x_i/sqrt(i))) + 1
[x1, x2] = meshgrid(-600:10:600, -600:10:600);
F4 = 1 + 1/4000 * (x1.^2 + x2.^2) - cos(x1 ./ sqrt(1)) .* cos(x2 ./ sqrt(2)) + 1;

figure;
mesh(x1, x2, F4);
title('Function F_4');
xlabel('x_1');
ylabel('x_2');
zlabel('F_4(x)');
colorbar;

% Function 5: F5(x) = sum_{i=1}^{7}[(x-a_i)(x-a_i)^T + c_i]^{-1}
a = [-32, -16, 0, 16, 32];
c = [10, 20, 30, 40, 50];
[x1, x2] = meshgrid(-32:1:32, -32:1:32);
F5 = zeros(size(x1));
for i = 1:length(a)
    F5 = F5 + 1 ./ ((x1 - a(i)).^2 + (x2 - a(i)).^2 + c(i));
end

figure;
mesh(x1, x2, F5);
title('Function F_5');
xlabel('x_1');
ylabel('x_2');
zlabel('F_5(x)');
colorbar;

% Function 6: F6(x) = sum_{i=1}^{n} -x_i \sin(\sqrt{|x_i|})
[x1, x2] = meshgrid(-500:10:500, -500:10:500);
F6 = -x1 .* sin(sqrt(abs(x1))) - x2 .* sin(sqrt(abs(x2)));

figure;
mesh(x1, x2, F6);
title('Function F_6');
xlabel('x_1');
ylabel('x_2');
zlabel('F_6(x)');
colorbar;
