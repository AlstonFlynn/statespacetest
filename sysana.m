classdef sysana
    % 线性系统分析器
    % 提供系统稳定性、能控性、能观性和模态分析的静态方法

    methods(Static)

        function stability = stabilityAnalysis(A, isContinuous, varargin)
            % 系统稳定性分析
            % 输入参数:
            %   A - 系统矩阵
            %   isContinuous - 是否为连续系统 (true:连续, false:离散)
            %   varargin - 可选参数: 采样时间Ts (离散系统需要)
            % 返回值:
            %   stability - 结构体，包含稳定性分析结果

            % 解析输入参数
            Ts = [];
            if ~isContinuous && nargin > 2
                Ts = varargin{1}; % 离散系统采样时间
            end

            % 计算特征值
            eigVals = eig(A);

            % 初始化返回结构
            stability = struct();
            stability.eigenvalues = eigVals;
            stability.isContinuous = isContinuous;

            if isContinuous
                % 连续系统：检查实部是否都小于0
                realParts = real(eigVals);
                stability.isStable = all(realParts < 0);
                stability.maxRealPart = max(realParts);
                stability.stablePoles = eigVals(realParts < 0);
                stability.unstablePoles = eigVals(realParts >= 0);
                stability.analysis = sprintf('连续系统稳定性: %s\n最大实部: %.4f', ...
                    iif(stability.isStable, '稳定', '不稳定'), stability.maxRealPart);
            else
                % 离散系统：检查模是否都小于1
                magnitudes = abs(eigVals);
                stability.isStable = all(magnitudes < 1);
                stability.maxMagnitude = max(magnitudes);
                stability.stablePoles = eigVals(magnitudes < 1);
                stability.unstablePoles = eigVals(magnitudes >= 1);
                stability.analysis = sprintf('离散系统稳定性: %s\n最大模: %.4f', ...
                    iif(stability.isStable, '稳定', '不稳定'), stability.maxMagnitude);

                % 如果有采样时间，计算等效连续系统特征值
                if ~isempty(Ts)
                    stability.continuousEquivalent = log(eigVals)/Ts;
                end
            end

            % 显示结果
            disp('=== 稳定性分析结果 ===');
            disp(stability.analysis);
            fprintf('特征值:\n');
            disp(eigVals);
        end

        function controllability = controllabilityAnalysis(A, B)
            % 系统能控性分析
            % 输入参数:
            %   A - 系统矩阵 (n×n)
            %   B - 输入矩阵 (n×m)
            % 返回值:
            %   controllability - 结构体，包含能控性分析结果

            n = size(A, 1); % 状态维度
            m = size(B, 2); % 输入维度

            % 构造能控性矩阵
            Co = B;
            for i = 1:n-1
                Co = [Co, A^i * B];
            end

            % 判断能控性
            rankCo = rank(Co);
            isControllable = (rankCo == n);

            % 计算能控性格拉姆矩阵（可选）
            % 注意：对于大系统，可能需要数值稳定性考虑

            % 返回结果
            controllability = struct();
            controllability.matrix = Co;
            controllability.rank = rankCo;
            controllability.isControllable = isControllable;
            controllability.stateDimension = n;

            if isControllable
                controllability.analysis = sprintf('系统完全能控\n能控性矩阵秩: %d (等于状态维度: %d)', ...
                    rankCo, n);
            else
                controllability.analysis = sprintf('系统不完全能控\n能控性矩阵秩: %d (小于状态维度: %d)', ...
                    rankCo, n);
                % 计算能控子空间
                [U, S, V] = svd(Co);
                controllability.controllableSubspace = U(:, 1:rankCo);
            end

            % 显示结果
            disp('=== 能控性分析结果 ===');
            disp(controllability.analysis);
            fprintf('能控性矩阵维度: %dx%d\n', size(Co, 1), size(Co, 2));
        end

        function observability = observabilityAnalysis(A, C)
            % 系统能观性分析
            % 输入参数:
            %   A - 系统矩阵 (n×n)
            %   C - 输出矩阵 (p×n)
            % 返回值:
            %   observability - 结构体，包含能观性分析结果

            n = size(A, 1); % 状态维度
            p = size(C, 1); % 输出维度

            % 构造能观性矩阵
            Ob = C;
            for i = 1:n-1
                Ob = [Ob; C * A^i];
            end

            % 判断能观性
            rankOb = rank(Ob);
            isObservable = (rankOb == n);

            % 返回结果
            observability = struct();
            observability.matrix = Ob;
            observability.rank = rankOb;
            observability.isObservable = isObservable;
            observability.stateDimension = n;

            if isObservable
                observability.analysis = sprintf('系统完全能观\n能观性矩阵秩: %d (等于状态维度: %d)', ...
                    rankOb, n);
            else
                observability.analysis = sprintf('系统不完全能观\n能观性矩阵秩: %d (小于状态维度: %d)', ...
                    rankOb, n);
                % 计算能观子空间的正交补
                [U, S, V] = svd(Ob');
                observability.unobservableSubspace = U(:, rankOb+1:end);
            end

            % 显示结果
            disp('=== 能观性分析结果 ===');
            disp(observability.analysis);
            fprintf('能观性矩阵维度: %dx%d\n', size(Ob, 1), size(Ob, 2));
        end

        function modal = modalAnalysis(A, isContinuous, varargin)
            % 系统模态分析
            % 输入参数:
            %   A - 系统矩阵
            %   isContinuous - 是否为连续系统
            %   varargin - 可选参数: 采样时间Ts (离散系统需要)
            % 返回值:
            %   modal - 结构体，包含模态分析结果

            % 解析输入参数
            Ts = [];
            if ~isContinuous && nargin > 2
                Ts = varargin{1};
            end

            % 特征值分解
            [V, D] = eig(A);
            eigenvalues = diag(D);

            % 初始化返回结构
            modal = struct();
            modal.eigenvalues = eigenvalues;
            modal.eigenvectors = V;
            modal.isContinuous = isContinuous;

            if isContinuous
                % 连续系统模态分析
                modal.naturalFrequencies = abs(eigenvalues);
                modal.dampingRatios = -real(eigenvalues) ./ abs(eigenvalues);

                % 对于实特征值，阻尼比为1或-1
                realIdx = imag(eigenvalues) == 0;
                modal.dampingRatios(realIdx) = sign(-real(eigenvalues(realIdx)));

                % 计算时间常数
                modal.timeConstants = -1 ./ real(eigenvalues);
                modal.timeConstants(real(eigenvalues) >= 0) = Inf;

                % 计算振荡频率（仅对复共轭对）
                modal.oscillationFreqs = abs(imag(eigenvalues)) / (2*pi);

                % 分类极点
                modal.stablePoles = eigenvalues(real(eigenvalues) < 0);
                modal.unstablePoles = eigenvalues(real(eigenvalues) >= 0);
                modal.realPoles = eigenvalues(imag(eigenvalues) == 0);
                modal.complexPoles = eigenvalues(imag(eigenvalues) ~= 0);

            else
                % 离散系统模态分析
                modal.magnitudes = abs(eigenvalues);
                modal.angles = angle(eigenvalues);

                % 如果有采样时间，计算等效连续系统参数
                if ~isempty(Ts)
                    s = log(eigenvalues)/Ts;
                    modal.equivalentContinuous = s;
                    modal.equivalentNaturalFreqs = abs(s);
                    modal.equivalentDampingRatios = -real(s) ./ abs(s);
                end

                % 分类极点
                modal.stablePoles = eigenvalues(abs(eigenvalues) < 1);
                modal.unstablePoles = eigenvalues(abs(eigenvalues) >= 1);
                modal.unitCirclePoles = eigenvalues(abs(abs(eigenvalues) - 1) < 1e-10);
            end

            % 计算模态参与因子（可选）
            modal.participationFactors = abs(V) .^ 2;

            % 显示结果
            disp('=== 模态分析结果 ===');
            fprintf('系统矩阵A的特征值:\n');
            for i = 1:length(eigenvalues)
                fprintf('  λ%d = ', i);
                if imag(eigenvalues(i)) >= 0
                    fprintf('%.4f + %.4fj\n', real(eigenvalues(i)), imag(eigenvalues(i)));
                else
                    fprintf('%.4f - %.4fj\n', real(eigenvalues(i)), -imag(eigenvalues(i)));
                end
            end

            if isContinuous
                disp('连续系统参数:');
                for i = 1:length(eigenvalues)
                    if imag(eigenvalues(i)) ~= 0
                        fprintf('  极点%d: 自然频率=%.4f rad/s, 阻尼比=%.4f', ...
                            i, modal.naturalFrequencies(i), modal.dampingRatios(i));
                        if modal.dampingRatios(i) > 0 && modal.dampingRatios(i) < 1
                            fprintf(', 振荡频率=%.4f Hz\n', modal.oscillationFreqs(i));
                        else
                            fprintf('\n');
                        end
                    end
                end
            end
        end

        function analyzeSystem(A, B, C, isContinuous, varargin)
            % 综合系统分析
            % 输入参数:
            %   A, B, C - 系统矩阵
            %   isContinuous - 是否为连续系统
            %   varargin - 可选参数

            disp('========================================');
            disp('        线性系统综合分析');
            disp('========================================');

            % 1. 稳定性分析
            if ~isempty(varargin) && ~isContinuous
                stability = sysana.stabilityAnalysis(A, isContinuous, varargin{1});
            else
                stability = sysana.stabilityAnalysis(A, isContinuous);
            end

            % 2. 能控性分析
            if ~isempty(B)
                controllability = sysana.controllabilityAnalysis(A, B);
            end

            % 3. 能观性分析
            if ~isempty(C)
                observability = sysana.observabilityAnalysis(A, C);
            end

            % 4. 模态分析
            if ~isempty(varargin) && ~isContinuous
                modal = sysana.modalAnalysis(A, isContinuous, varargin{1});
            else
                modal = sysana.modalAnalysis(A, isContinuous);
            end

            % 综合判断
            disp('========================================');
            disp('        系统综合属性');
            disp('========================================');

            isStable = stability.isStable;
            isCtrl = (~isempty(B) && controllability.isControllable);
            isObs = (~isempty(C) && observability.isObservable);

            fprintf('稳定性: %s\n', iif(isStable, '稳定', '不稳定'));
            if ~isempty(B)
                fprintf('能控性: %s\n', iif(isCtrl, '完全能控', '不完全能控'));
            end
            if ~isempty(C)
                fprintf('能观性: %s\n', iif(isObs, '完全能观', '不完全能观'));
            end

            if isStable && isCtrl && isObs
                disp('系统是稳定、能控、能观的');
            end
        end
    end
end

% 辅助函数：条件表达式
function result = iif(condition, trueValue, falseValue)
if condition
    result = trueValue;
else
    result = falseValue;
end
end