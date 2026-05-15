function [OrganizedMixSignal] = MixPixelSignal(ROISignal1,ROISignal2,ROISignal3)

% 配置基础参数
[Long PixelNums types] = size(ROISignal1);
RefSignal(:,1) = mean(ROISignal1(:,:,1),2);
RefSignal(:,2) = mean(ROISignal2(:,:,1),2);
RefSignal(:,3) = mean(ROISignal3(:,:,1),2);

Repeats = 12;
Colorss = ['r','b','r'];
% 提取混合信号的计算参数
MixedNum = [1:PixelNums-1];
OrganizedMixSignal1 = MixingROIsMask(ROISignal1,MixedNum,Repeats); % 生成随机混合的动脉+脑组织Mask
OrganizedMixSignal2 = MixingROIsMask(ROISignal2,MixedNum,Repeats); % 生成随机混合的动脉+静脉Mask
OrganizedMixSignal3 = MixingROIsMask(ROISignal3,MixedNum,Repeats); % 生成随机混合的脑组织+静脉Mask
OrganizedMixSignal = cat(4,OrganizedMixSignal1,OrganizedMixSignal2,OrganizedMixSignal3);

% 2. 计算i与所有j的相位差，集合存储到PhaseDiff(:, i, chooseN, ch)
PhaseDiff = ones([Long*(Repeats-1), Repeats, length(MixedNum), size(OrganizedMixSignal,4)]);
% PhaseDiff = ones([Long, Repeats, length(MixedNum), size(OrganizedMixSignal,4)]);

for ch = 1:size(OrganizedMixSignal,4)
    for chooseN = 1:length(MixedNum)
        for i = 1:Repeats
            sig_i = OrganizedMixSignal(:, i, chooseN,ch);  % 第i个信号
            % 创建临时矩阵存储当前i对应的所有j的相位差（每列一个j）
            temp_phase = zeros(Long, Repeats-1) .* NaN;  % 维度[Long, Repeats]
            k = 1;
            for j = 1:Repeats
                if i ~= j
                    sig_j = OrganizedMixSignal(:, j, chooseN,ch);  % 第j个信号
                    temp_phase(:, k) = PhaseDifference(sig_i, sig_j);  % 第j列存储i与j的相位差
%                     temp_phase(:, k) = PhaseDifference(RefSignal(:,ch), sig_j);  % 第j列存储i与j的相位差
                    k = k + 1;
                end
            end
            % 将所有j的相位差整体赋值给PhaseDiff的第i个位置
            PhaseDiff(:, i, chooseN,ch) = temp_phase(:);
        end
    end
end

%% subfunctions
% 生成混合信号的掩膜
function [OrganizedMixSignal] = MixingROIsMask(ROISignal,ChooseNum,Repeats)
    
    FixedROISignal = ROISignal(:,:,1);
    MovedROISignal = ROISignal(:,:,2);
    
    PixelNums = size(FixedROISignal,2);
    
    OrganizedMixSignal = zeros([size(FixedROISignal,1) Repeats length(ChooseNum)]);
    for chooseN = 1:length(ChooseNum)
        for rdm = 1:Repeats
            MovedSelectIndex = transpose(randperm(PixelNums,ChooseNum(chooseN)));
            FixedSelectIndex = transpose(randperm(PixelNums,PixelNums - ChooseNum(chooseN)));
            MixedSignal = [FixedROISignal(:,FixedSelectIndex),MovedROISignal(:,MovedSelectIndex)];
            OrganizedMixSignal(:,rdm,chooseN) = mean(MixedSignal(:,randperm(PixelNums,12)),2);
%             OrganizedMixSignal(:,rdm,chooseN) = mean(MixedSignal,2);
        end
    end
end

% 计算信号之间的相位差
function [phase_diff] = PhaseDifference(signal1,signal2)
    
    % 1. 希尔伯特变换得到瞬时相位
    hSignal1 = hilbert(signal1);
    hSignal2 = hilbert(signal2);
    phase1 = unwrap(angle(hSignal1)); % 瞬时相位
    phase2 = unwrap(angle(hSignal2)); % 瞬时相位
    
    % 计算相位差
    phase_diff = phase1 - phase2;

end

%% 绘图程序
Seq = [2,4,6,9,13];
% Seq = [2,5,9,14,24];

ks = 3;
rows = 5;  % 行数 
cols = ks;  % 列数
margin = 0.001;  % 边缘压缩量（固定值）

figure(10); set(gcf,'position',[144,66,931,1246]);
for chooseN = 1:length(Seq)
    for ch = 1
        subplot(5,ks,[chooseN*ks-ks+1:chooseN*ks-ks+3]);
        plot(OrganizedMixSignal(:,:,Seq(chooseN),ch),'Color',Colorss(ch),'LineWidth',0.8);
        alpha(0.5);
        hold on;
        xlim([301 1100]);
        set(gca, 'YTick', [],'XTick', []);  % 移除刻度线
        set(gca, 'XColor', 'none', 'YColor', 'none');  % 可选：隐藏坐标轴线条（若需要）
    end
    pos = get(gca, 'Position');
    new_pos = [pos(1) + margin, pos(2) + margin, pos(3)+ 50*margin, pos(4)+ 50*margin];
    set(gca, 'Position', new_pos);
end

ks = 1;
for ch = 1:size(OrganizedMixSignal,4)
    figure(11+ch); set(gcf,'position',[144,66,347,1246]);
    for chooseN = 1:length(Seq)
        subplot(5,ks,chooseN);
        PhaseDiffX = PhaseDiff(:,:,chooseN,ch);
        PhaseDiffX1 = PhaseDiffX(:);
        
        % 1. 设置角度范围（例如0~120度）和目标Bin数量（例如30个，比默认更多）
        angle_min = -40;       % 角度范围最小值（度）
        angle_max = 40;     % 角度范围最大值（度）
        num_bins = 60;       % 增加Bin数量（值越大，分箱越细）
        
        % 2. 绘制极坐标直方图，指定Bin数量
        % circ_plot的第3个参数用于指定Bin边缘，第4个参数为Bin数量
        circ_plot(PhaseDiffX1, 'hist', [], num_bins, true, true, ...
            'linewidth', 3, 'color', 'b');
        
        % 3. 获取极坐标轴句柄并修改角度范围
        pax = gca;
        pax.ThetaAxisUnits = 'degrees';  % 确保角度单位为度
        pax.ThetaLim = [angle_min angle_max];  % 设置角度范围为0~120度
        
        % 4. 可选：调整角度刻度以匹配新范围
        pax.ThetaTick = linspace(angle_min, angle_max, 5);  % 5个刻度点（0,30,60,90,120）
        
        pax.ThetaColor = 'k';
        pax.RColor = 'b';
        pax.LineWidth = 3;
        pax.GridColor = 'r';
        pax.RLim = [0, 0.12];  % 径向范围留20%余量
        pax.ThetaTickLabel = {};  % 替换原有的 pax.XTickLabel = {};
        drawnow;
        set(gca,'FontSize',20,'FontName','Arial','FontWeight','Bold');
        
        pos = get(gca, 'Position');
%         new_pos = [pos(1) + margin, pos(2) + margin, pos(3) - 2*margin, pos(4) - 2*margin];
        new_pos = [pos(1) + margin, pos(2) + margin, pos(3) + 20*margin, pos(4) + 20*margin];
        set(gca, 'Position', new_pos);
    end
end
end