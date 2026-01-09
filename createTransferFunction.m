function G = generateTransferFunction(num, den)
% 生成传递函数G
% 输入参数:
%   num - 传递函数分子多项式系数向量
%   den - 传递函数分母多项式系数向量
% 输出参数:
%   G   - 传递函数对象

% 如果没有提供参数，使用默认的二阶系统
if nargin < 2
    % 默认传递函数: G(s) = 1/(s^2 + 2ζωs + ω^2)
    % 其中 ω = 2, ζ = 0.7
    num = [1];
    den = [1, 2 * 0.7 * 2, 2^2];
    fprintf('使用默认传递函数: G(s) = 1/(s^2 + %.2fs + %.2f)\n', 2 * 0.7 * 2, 2^2);
else
    fprintf('使用自定义传递函数: G(s) = (');
    for i = 1:length(num)
        if i < length(num)
            fprintf('%.2fs^%d + ', num(i), length(num)-i);
        else
            fprintf('%.2f)/(', num(i));
        end
    end

    for i = 1:length(den)
        if i < length(den)
            fprintf('%.2fs^%d + ', den(i), length(den)-i);
        else
            fprintf('%.2f)\n', den(i));
        end
    end
end

% 创建传递函数
G = tf(num, den);

% 显示传递函数
disp('生成的传递函数:');
G

% 绘制阶跃响应
figure;
step(G);
title('传递函数阶跃响应');
grid on;

% 绘制波特图
figure;
bode(G);
title('传递函数波特图');
grid on;
end