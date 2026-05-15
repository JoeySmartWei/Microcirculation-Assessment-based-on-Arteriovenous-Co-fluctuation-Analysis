function [SelectSignal,PhaseDiff] = PurifyPixelSignalPhase(PixelSignal)
% 从像素组成的信号中随机提取信号，并对其均值、方差、与参考信号的相似性进行计算
[Long PixelNums Channels] = size(PixelSignal);
Repeats = 10;
ChooseNum = [1:PixelNums-1];
RefSignal = squeeze(mean(PixelSignal,2));
Colorss = ['r','g','b'];
Vessels = {'Artery','Tissue','Vein','ArteryTissue','VeinTissue'};

% 初始化SelectSignal和PhaseDiff（维度：[信号长度, 重复次数, 抽样数量, 通道]）
SelectSignal = ones([Long, Repeats, length(ChooseNum), Channels]);
PhaseDiff = ones([Long*(Repeats-1), Repeats, length(ChooseNum), Channels]);

for ch = 1:Channels
    for chooseN = 1:length(ChooseNum)
        % 1. 提取当前参数下的所有Repeats个信号
        for rdm = 1:Repeats
            SelectIndex = transpose(randperm(PixelNums, ChooseNum(chooseN)));
            SelectSignal(:, rdm, chooseN, ch) = mean(PixelSignal(:, SelectIndex, ch), 2);
        end
        
        % 2. 计算i与所有j的相位差，集合存储到PhaseDiff(:, i, chooseN, ch)
        for i = 1:Repeats
            sig_i = SelectSignal(:, i, chooseN, ch);  % 第i个信号
            % 创建临时矩阵存储当前i对应的所有j的相位差（每列一个j）
            temp_phase = zeros(Long, Repeats-1) .* NaN;  % 维度[Long, Repeats]
            k = 1;
            for j = 1:Repeats
                if i ~= j
                    sig_j = SelectSignal(:, j, chooseN, ch);  % 第j个信号
                    temp_phase(:, k) = PhaseDifference(sig_i, sig_j);  % 第j列存储i与j的相位差
                    k = k + 1;
                end
            end
            % 将所有j的相位差整体赋值给PhaseDiff的第i个位置
            PhaseDiff(:, i, chooseN, ch) = temp_phase(:);
        end
    end
end


    function [phase_diff] = PhaseDifference(signal1,signal2)
        
        % 1. 希尔伯特变换得到瞬时相位
        hSignal1 = hilbert(signal1);
        hSignal2 = hilbert(signal2);
        phase1 = unwrap(angle(hSignal1)); % 瞬时相位
        phase2 = unwrap(angle(hSignal2)); % 瞬时相位
        
        % 计算相位差
        phase_diff = phase1 - phase2;
        
    end

%% 绘图

Seq = [2,5,9,12,20];

ks = 3;
rows = 5;  % 行数 
cols = ks;  % 列数
margin = 0.001;  % 边缘压缩量（固定值）

figure(10); set(gcf,'position',[144,66,931,1246]);
for chooseN = 1:length(Seq)
    for ch = 1:Channels
        subplot(5,ks,[chooseN*ks-ks+1:chooseN*ks-ks+3]);
        plot(SelectSignal(:,:,Seq(chooseN),ch),'Color',Colorss(ch),'LineWidth',0.8);
        alpha(0.5);
        hold on;
        xlim([1001 1800]);
        set(gca, 'YTick', [],'XTick', []);  % 移除刻度线
        set(gca, 'XColor', 'none', 'YColor', 'none');  % 可选：隐藏坐标轴线条（若需要）
    end
    pos = get(gca, 'Position');
    new_pos = [pos(1) + margin, pos(2) + margin, pos(3)+ 50*margin, pos(4)+ 50*margin];
    set(gca, 'Position', new_pos);
end

ks = 1;
for ch = 1:Channels
    figure(11+ch); set(gcf,'position',[144,66,347,1246]);
    for chooseN = 1:length(Seq)
        subplot(5,ks,chooseN);
        PhaseDiffX = PhaseDiff(:,:,chooseN,ch);
        PhaseDiffX1 = PhaseDiffX(:);
        circ_plot(PhaseDiffX1,'hist',[],20,true,true,'linewidth',3,'color','b');
        pax = gca;
        pax.ThetaColor = 'k';
        pax.RColor = 'b';
        pax.LineWidth = 3;
        pax.GridColor = 'r';
        pax.RLim = [0, 0.14];  % 径向范围留20%余量
        pax.ThetaTickLabel = {};  % 替换原有的 pax.XTickLabel = {};
        drawnow;
        set(gca,'FontSize',20,'FontName','Arial','FontWeight','Bold');
        
        pos = get(gca, 'Position');
%         new_pos = [pos(1) + margin, pos(2) + margin, pos(3) - 2*margin, pos(4) - 2*margin];
        new_pos = [pos(1) + margin, pos(2) + margin, pos(3)+ 40*margin, pos(4)+ 40*margin];
        set(gca, 'Position', new_pos);
    end
end

end