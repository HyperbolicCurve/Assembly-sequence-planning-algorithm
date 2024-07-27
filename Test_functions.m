% Define search ranges
ranges = {
    [-100, 100],   % F1
    [-100, 100],   % F2
    [-1.28, 1.28], % F3
    [-600, 600],   % F4
    [-32, 32],     % F5
    [-500, 500]    % F6
};

% Define dimensions
dims = [30, 30, 30, 30, 4, 30];

% Loop over each function and plot
for i = 1:6
    range = ranges{i};
    dim = dims(i);
    
    % Generate grid data within the range for 2D plot
    [x1, x2] = meshgrid(linspace(range(1), range(2), 100), linspace(range(1), range(2), 100));
    
    % Calculate function output
    switch i
        case 1
            z = arrayfun(@(a, b) F1([a, b]), x1, x2);
        case 2
            z = arrayfun(@(a, b) F2([a, b]), x1, x2);
        case 3
            z = arrayfun(@(a, b) F3([a, b]), x1, x2);
        case 4
            z = arrayfun(@(a, b) F4([a, b]), x1, x2);
        case 5
            z = arrayfun(@(a, b) F5([a, b]), x1, x2);
        case 6
            z = arrayfun(@(a, b) F6([a, b]), x1, x2);
    end
    
    % Plot function output in a new figure
    figure;
    surf(x1, x2, z);
    title(['Function F', num2str(i)]);
    xlabel('x_1');
    ylabel('x_2');
    zlabel(['F', num2str(i), '(x_1, x_2)']);
end

% Test functions
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
    o = sum([1:dim] .* (x.^4)) + rand;
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
