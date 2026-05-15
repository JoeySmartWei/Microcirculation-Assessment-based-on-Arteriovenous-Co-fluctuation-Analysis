function [CorrResults, PhaseDiffResults] = MixBigPixelSignalPackage(ROISignal)

% 配置基础参数
[Long, PixelNums, types] = size(ROISignal);

% 定义混合比例梯度（5% 递进，从5%到95%）
MixingRatios = 0.00:0.02:1.00;  % 静脉信号占比，动脉占比为1-该值
numRatios = length(MixingRatios);

% 定义分组大小序列（与Purify程序保持一致）
GroupSizes = [1, 2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100];
GroupSizes = GroupSizes(GroupSizes <= PixelNums);  % 过滤无效分组大小
if isempty(GroupSizes)
    error('所有分组大小均超过像素总数，请调整GroupSizes');
end
numGroupSizes = length(GroupSizes);

% 获取动脉参考信号
ReferenceSignal = mean(ROISignal(:,:,1),2); 

% 初始化结果矩阵（维度：组数×混合比例数×分组大小数）
CorrResults = cell(numGroupSizes, 1);       % 每个元素为"组数×混合比例"的相关性矩阵
PhaseDiffResults = cell(numGroupSizes, 1);  % 每个元素为"组数×混合比例"的相位差矩阵

% 按不同混合比例处理
for ratioIdx = 1:numRatios
    veinRatio = MixingRatios(ratioIdx);  % 静脉信号占比
    arteryRatio = 1 - veinRatio;        % 动脉信号占比

    % 生成当前比例的混合信号（临时变量，不长期存储）
    MixingSignal = generateMixedSignal(ROISignal, PixelNums, arteryRatio);

    % 对当前混合信号进行分组处理并计算相关性和相位差
    for g = 1:numGroupSizes
        groupSize = GroupSizes(g);
        groupCount = floor(PixelNums / groupSize);  % 计算有效组数
        if groupCount < 1
            warning('混合比例%.0f%%对应的分组大小%d有效组数不足，已跳过', veinRatio*100, groupSize);
            if ratioIdx == 1  % 初始化空矩阵
                CorrResults{g} = NaN(groupCount, numRatios);
                PhaseDiffResults{g} = NaN(groupCount, numRatios);
            end
            continue;
        end

        % 随机打乱像素并分配到各组
        shuffledIndices = randperm(PixelNums);
        groupSignals = zeros(Long, groupCount);  % 临时存储当前分组的均值
        for i = 1:groupCount
            startIdx = (i-1)*groupSize + 1;
            endIdx = i*groupSize;
            currentIndices = shuffledIndices(startIdx:endIdx);
            groupSignals(:, i) = mean(MixingSignal(:, currentIndices), 2);  % 计算组均值
        end

        % 计算并存储相关性和相位差（按分组大小和混合比例索引）
        if ratioIdx == 1  % 初始化矩阵
            CorrResults{g} = zeros(groupCount, numRatios);
            PhaseDiffResults{g} = zeros(groupCount, numRatios);
        end
        CorrResults{g}(:, ratioIdx) = corr(ReferenceSignal, groupSignals);
        PhaseDiffResults{g}(:, ratioIdx) = PhaseDifference(ReferenceSignal, groupSignals);
        
        % 对当前混合比例下的结果排序（保持原逻辑）
        CorrResults{g}(:, ratioIdx) = sort(CorrResults{g}(:, ratioIdx));
        PhaseDiffResults{g}(:, ratioIdx) = sort(PhaseDiffResults{g}(:, ratioIdx));
    end

%     fprintf('已完成混合比例 %.0f%%（静脉）的处理\n', veinRatio*100);
end

%% 子函数：生成指定比例的混合信号
function MixingSignal = generateMixedSignal(ROISignal, PixelNums, arteryRatio)
    FixedROISignal = ROISignal(:,:,1);  % 动脉信号
    MovedROISignal = ROISignal(:,:,2);  % 静脉信号
    
    % 计算每种信号需要的像素数量
    arteryPixels = round(PixelNums * arteryRatio);
    veinPixels = PixelNums - arteryPixels;  % 确保总像素数不变
    
    % 随机选择动脉和静脉像素
    arteryIndices = randperm(PixelNums, arteryPixels);
    veinIndices = setdiff(1:PixelNums, arteryIndices);  % 避免重复选择
    
    % 按比例混合信号
    arteryMixed = FixedROISignal(:, arteryIndices);
    veinMixed = MovedROISignal(:, veinIndices);
    
    % 扩展为矩阵形式并随机打乱像素顺序
    MixingSignal = [arteryMixed,veinMixed];
    MixingSignal = MixingSignal(:,randperm(size(FixedROISignal, 2)));
end

function [average_phase] = PhaseDifference(signal1,signal2)

if nargin == 2
    % 1. 希尔伯特变换得到瞬时相位
    hSignal1 = hilbert(signal1);
    hSignal2 = hilbert(signal2);
    phase1 = unwrap(angle(hSignal1)); % 瞬时相位
    phase2 = unwrap(angle(hSignal2)); % 瞬时相位
    
    % 计算相位差
    phase_diff = phase1 - phase2;
else
    phase_diff = unwrap(angle(hilbert(signal1)));
end

% 将相位转换到复数域
complex_sequence = exp(1i * phase_diff);

% 计算平均
average_complex = mean(complex_sequence);

% 转换回相位
average_phase = abs(angle(average_complex));

end
end