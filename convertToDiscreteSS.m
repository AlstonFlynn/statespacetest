function [Ad, Bd, Cd, Dd, sysd] = convertToDiscreteSS(G, Ts, method)
% 将连续传递函数转换为离散状态空间模型
% 输入参数:
%   G      - 连续传递函数对象
%   Ts     - 采样时间（秒），默认为0.1
%   method - 离散化方法: 'zoh'（零阶保持器）, 'foh'（一阶保持器）,
%           'tustin'（双线性变换）, 默认为'zoh'
% 输出参数:
%   Ad, Bd, Cd, Dd - 离散状态空间矩阵
%   sysd - 离散状态空间系统对象

% 设置默认值
if nargin < 3
    method = 'zoh';
end
if nargin < 2
    Ts = 0.1;
end

fprintf('离散化参数:\n');
fprintf('  采样时间 Ts = %.3f 秒\n', Ts);
fprintf('  离散化方法: %s\n', method);

% 转换为连续状态空间
sysc = ss(G);
Ac = sysc.A;
Bc = sysc.B;
Cc = sysc.C;
Dc = sysc.D;

fprintf('连续状态空间模型:\n');
disp('A矩阵:');
disp(Ac);
disp('B矩阵:');
disp(Bc);
disp('C矩阵:');
disp(Cc);
disp('D矩阵:');
disp(Dc);

% 离散化状态空间模型
sysd = c2d(sysc, Ts, method);

% 获取离散状态空间矩阵
Ad = sysd.A;
Bd = sysd.B;
Cd = sysd.C;
Dd = sysd.D;

fprintf('\n离散状态空间模型:\n');
disp('A矩阵:');
disp(Ad);
disp('B矩阵:');
disp(Bd);
disp('C矩阵:');
disp(Cd);
disp('D矩阵:');
disp(Dd);

% 计算和显示系统特性
fprintf('系统特性:\n');
fprintf('  系统阶数: %d\n', size(Ad,1));
fprintf('  采样频率: %.2f Hz\n', 1/Ts);
fprintf('  奈奎斯特频率: %.2f Hz\n', 1/(2*Ts));

% 检查稳定性
eig_cont = eig(Ac);
eig_disc = eig(Ad);

fprintf('  连续系统特征值: ');
fprintf('%.3f ', real(eig_cont));
if all(real(eig_cont) < 0)
    fprintf('(稳定)\n');
else
    fprintf('(不稳定)\n');
end

fprintf('  离散系统特征值: ');
fprintf('%.3f ', abs(eig_disc));
if all(abs(eig_disc) < 1)
    fprintf('(稳定)\n');
else
    fprintf('(不稳定)\n');
end

% 绘制离散系统的阶跃响应
figure;
step(sysd, 5); % 显示5秒的响应
title(sprintf('离散系统阶跃响应 (Ts=%.3f秒)', Ts));
grid on;

% 绘制零极点图
figure;
pzmap(sysd);
title(sprintf('离散系统零极点图 (Ts=%.3f秒)', Ts));
grid on;
end