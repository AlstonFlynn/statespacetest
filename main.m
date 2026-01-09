%% 主程序：生成传递函数并转换为离散状态空间模型
clear;
close all;
clc;

fprintf('========================================\n');
fprintf('传递函数转换为离散状态空间模型\n');
fprintf('========================================\n\n');

%% 第一部分：生成传递函数
fprintf('1. 生成传递函数\n');
fprintf('----------------------------------------\n');

% 可以选择使用默认传递函数或自定义传递函数
useDefault = input('是否使用默认传递函数? (y/n, 默认: y): ', 's');
if isempty(useDefault) || lower(useDefault) == 'y'
    % 使用默认传递函数
    G = generateTransferFunction();
else
    % 自定义传递函数参数
    fprintf('\n输入传递函数参数:\n');

    % 输入分子多项式系数
    num_str = input('输入分子系数向量 (例如: [1, 2] 表示 s+2): ', 's');
    num = str2num(num_str); %#ok<ST2NM>

    % 输入分母多项式系数
    den_str = input('输入分母系数向量 (例如: [1, 3, 2] 表示 s^2+3s+2): ', 's');
    den = str2num(den_str); %#ok<ST2NM>

    % 生成传递函数
    G = generateTransferFunction(num, den);
end

fprintf('\n按任意键继续...\n');
pause;

%% 第二部分：转换为离散状态空间模型
fprintf('\n\n2. 转换为离散状态空间模型\n');
fprintf('----------------------------------------\n');

% 设置离散化参数
Ts = input('输入采样时间 Ts (秒) (默认: 0.1): ');
if isempty(Ts)
    Ts = 0.1;
end

fprintf('可选离散化方法:\n');
fprintf('  1. zoh - 零阶保持器\n');
fprintf('  2. foh - 一阶保持器\n');
fprintf('  3. tustin - 双线性变换 (Tustin)\n');
method_choice = input('选择离散化方法 (1/2/3, 默认: 1): ', 's');

if isempty(method_choice)
    method = 'zoh';
else
    switch method_choice
        case '1'
            method = 'zoh';
        case '2'
            method = 'foh';
        case '3'
            method = 'tustin';
        otherwise
            method = 'zoh';
            fprintf('无效选择，使用默认方法: zoh\n');
    end
end

% 转换为离散状态空间模型
[Ad, Bd, Cd, Dd, sysd] = convertToDiscreteSS(G, Ts, method);

%% 第三部分：系统分析
fprintf('\n\n3. 系统分析\n');
fprintf('----------------------------------------\n');

% 比较连续和离散系统的阶跃响应
figure;
t_cont = 0:0.01:5;  % 连续时间
t_disc = 0:Ts:5;     % 离散时间

% 计算连续系统响应
y_cont = step(G, t_cont);

% 计算离散系统响应
y_disc = step(sysd, t_disc);

% 绘制比较图
subplot(2,1,1);
plot(t_cont, y_cont, 'b-', 'LineWidth', 1.5);
hold on;
stairs(t_disc, y_disc, 'r--', 'LineWidth', 1.5);
xlabel('时间 (秒)');
ylabel('幅值');
title('连续与离散系统阶跃响应比较');
legend('连续系统', '离散系统', 'Location', 'best');
grid on;

subplot(2,1,2);
% 计算误差
t_common = 0:0.01:5;
y_cont_interp = interp1(t_cont, y_cont, t_common);
y_disc_interp = interp1(t_disc, y_disc, t_common, 'previous');
error = y_cont_interp - y_disc_interp;
plot(t_common, error, 'g-', 'LineWidth', 1.5);
xlabel('时间 (秒)');
ylabel('误差');
title('离散化误差');
grid on;

% 计算性能指标
fprintf('性能指标:\n');
% 连续系统
info_cont = stepinfo(G);
fprintf('连续系统:\n');
fprintf('  上升时间: %.3f 秒\n', info_cont.RiseTime);
fprintf('  超调量: %.2f%%\n', info_cont.Overshoot);
fprintf('  调节时间: %.3f 秒\n', info_cont.SettlingTime);

% 离散系统
info_disc = stepinfo(sysd);
fprintf('离散系统:\n');
fprintf('  上升时间: %.3f 秒\n', info_disc.RiseTime);
fprintf('  超调量: %.2f%%\n', info_disc.Overshoot);
fprintf('  调节时间: %.3f 秒\n', info_disc.SettlingTime);

%% 第四部分：保存结果
fprintf('\n\n4. 保存结果\n');
fprintf('----------------------------------------\n');

save_results = input('是否保存结果到文件? (y/n, 默认: n): ', 's');
if lower(save_results) == 'y'
    filename = input('输入文件名 (不带扩展名, 默认: system_data): ', 's');
    if isempty(filename)
        filename = 'system_data';
    end

    % 保存数据
    save(filename, 'G', 'Ad', 'Bd', 'Cd', 'Dd', 'sysd', 'Ts', 'method');
    fprintf('结果已保存到 %s.mat\n', filename);

    % 生成报告文件
    report_filename = [filename, '_report.txt'];
    fid = fopen(report_filename, 'w');
    fprintf(fid, '系统离散化报告\n');
    fprintf(fid, '================\n\n');
    fprintf(fid, '1. 传递函数\n');
    fprintf(fid, '   G(s) = %s\n\n', evalc('G'));
    fprintf(fid, '2. 离散化参数\n');
    fprintf(fid, '   采样时间: %.4f 秒\n', Ts);
    fprintf(fid, '   离散化方法: %s\n\n', method);
    fprintf(fid, '3. 离散状态空间模型\n');
    fprintf(fid, '   A = \n');
    for i = 1:size(Ad,1)
        fprintf(fid, '     ');
        fprintf(fid, '%.6f ', Ad(i,:));
        fprintf(fid, '\n');
    end
    fprintf(fid, '\n   B = \n');
    fprintf(fid, '     ');
    fprintf(fid, '%.6f ', Bd');
    fprintf(fid, '\n\n');
    fprintf(fid, '   C = \n');
    fprintf(fid, '     ');
    fprintf(fid, '%.6f ', Cd);
    fprintf(fid, '\n\n');
    fprintf(fid, '   D = %.6f\n\n', Dd);
    fprintf(fid, '4. 系统特性\n');
    fprintf(fid, '   系统阶数: %d\n', size(Ad,1));
    fprintf(fid, '   采样频率: %.2f Hz\n', 1/Ts);
    fprintf(fid, '   奈奎斯特频率: %.2f Hz\n', 1/(2*Ts));
    fclose(fid);
    fprintf('报告已保存到 %s\n', report_filename);
end

fprintf('\n程序执行完成！\n');