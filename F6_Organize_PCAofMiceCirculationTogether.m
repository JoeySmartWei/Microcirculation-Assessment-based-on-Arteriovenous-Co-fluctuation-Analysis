%% PCA Analysis for Mice Circulation Data
% Author: Bochao Niu
% Description: This script performs PCA analysis on mice CBF (Cerebral Blood Flow) data,
% extracts regional brain signals, and generates statistical plots/activity maps.

%% ========================== Configuration (Constants) ==========================
clearvars; close all; clc;

% Experimental groups (Stroke/Control)
EXP_GROUPS = {'Stroke','Control'};
% Mouse number sequences for each group
MOUSE_SEQUENCES = {[1:19], [1:12,14:31]};
% PCA parameter categories
PCA_PARAMS = {'CirculationParas','CBFSignalParas'};
% PCA analysis types (General/Top/Bottom)
ANALYSIS_TYPES = {'General','Top','Bottom'};

% File path constants (MODIFY THESE PATHS ACCORDING TO YOUR ENVIRONMENT)
ATLAS_PATH = 'E:\CBF-data\GongData\VesselsSegmentation\AtlasROIsNetwork\Mouse_Atlas_RoisNew.mat';
CBF_FILE_PATH = 'E:\CBF-data\GongData\OutPut\CBFPackageAffine\';
PCA_RESULT_PATH = 'E:\CBF-data\GongData\OutPut\PCAResultsAffine\';
FIGURE_OUTPUT_PATH = fullfile(PCA_RESULT_PATH, 'Figure');
mkdir(FIGURE_OUTPUT_PATH); % Ensure output directory exists

% Sampling frequency for time series analysis
SAMPLING_FREQ = 10; % Hz
% Gaussian filter sigma for activity map smoothing
GAUSS_FILTER_SIGMA = 2.5;

%% ========================== Brain Region Definition ==========================
% Load brain atlas ROI masks
fprintf('Loading brain atlas ROI masks...\n');
load(ATLAS_PATH);

% Full brain region list (left/right hemispheres + whole brain)
brainRegions = {
    'L_Olfactory','L_Prelimbic','L_Cingulate','L_FrontalAssociation','L_M1','L_M2',...
    'L_FrontalArea3','L_Trunk','L_Hindlimb','L_Shoulder','L_Forelimb','L_Barrel','L_Head',...
    'L_V2L','L_V1','L_V2M','L_Retrosplenial','L_Medial','L_Lateral','L_Posterior',...
    'L_Auditory','L_Association','L_Superior','L_Inferior',...
    'R_Olfactory','R_Prelimbic','R_Cingulate','R_FrontalAssociation','R_M1','R_M2',...
    'R_FrontalArea3','R_Trunk','R_Hindlimb','R_Shoulder','R_Forelimb','R_Barrel','R_Head',...
    'R_V2L','R_V1','R_V2M','R_Retrosplenial','R_Medial','R_Lateral','R_Posterior',...
    'R_Auditory','R_Association','R_Superior','R_Inferior',...
    'L_Limbic','L_Motor','L_Parietal','L_Visual','L_Colliculi','L_Temporal','L_Somatosensory',...
    'R_Limbic','R_Motor','R_Parietal','R_Visual','R_Colliculi','R_Temporal','R_Somatosensory',...
    'L_Brain','R_Brain','Whole_Brain'
};

% Bilateral region groups (combined left+right)
bilateralRegions1 = {
    'Olfactory','Prelimbic','Cingulate','FrontalAssociation','M1','M2',...
    'FrontalArea3','Trunk','Hindlimb','Shoulder','Forelimb','Barrel','Head','Auditory','Association',...
    'V2L','V1','V2M','Retrosplenial','Medial','Lateral','Posterior','Superior','Inferior','Motor','Somatosensory','Visual'
};

bilateralRegions2 = {
    'Forelimb','Barrel','Head','M1','M2','Parietal','V2L','V1','V2M','Colliculi','Temporal','Limbic'
};

%% ========================== Preprocess ROI Masks ==========================
% Resize masks and remove empty rows (optimize memory)
fprintf('Preprocessing brain ROI masks...\n');
processedMasks = preprocessROIMasks(brainRegions);

% Combine left/right hemispheres for bilateral regions
bilateralMasks = combineBilateralRegions(processedMasks, bilateralRegions1, bilateralRegions2);

%% ========================== Main PCA Analysis Pipeline ==========================
% Process each experimental group (start from index 2: Control group in original code)
for groupIdx = 2:length(EXP_GROUPS)
    currentGroup = EXP_GROUPS{groupIdx};
    mouseNumbers = MOUSE_SEQUENCES{groupIdx};
    
    % Initialize result storage structure
    groupResults = initGroupResults(currentGroup, PCA_PARAMS, ANALYSIS_TYPES, ...
        length(mouseNumbers), bilateralRegions1, bilateralRegions2);
    
    % Process each mouse in the group
    for mouseIdx = 1:length(mouseNumbers)
        mouseNum = mouseNumbers(mouseIdx);
        fprintf('Processing %s Group: Mouse %d CBF Image ......\n', currentGroup, mouseNum);
        tic;

        % Load PCA results and time series data
        [pcaData, timeSeriesData] = loadMouseData(PCA_RESULT_PATH, currentGroup, mouseNum);
        
        % Calculate vascular percentage (Artery/Vein/Mixed)
        vascularPercent = calculateVascularPercentage(timeSeriesData, processedMasks.Whole_Brain);
        groupResults.VesselsPercentage(mouseIdx, :) = vascularPercent;
        
        % Extract PCA features for each parameter type
        for paramIdx = 1:length(PCA_PARAMS)
            paramName = PCA_PARAMS{paramIdx};
            pcaParamData = pcaData.(['Mice', num2str(mouseNum), '_', currentGroup, paramName]);
            
            % Extract features for each analysis type (General/Top/Bottom)
            groupResults = extractPCAFEatures(groupResults, paramName, pcaParamData, ...
                ANALYSIS_TYPES, mouseIdx, processedMasks, bilateralMasks, ...
                bilateralRegions1, bilateralRegions2);
            
            % Calculate correlations between analysis types
            groupResults = calculateTypeCorrelations(groupResults, paramName, mouseIdx, ANALYSIS_TYPES);
        end
        
        % Clean up temporary variables
        clearvars -regexp ^Mice\d+;
        fprintf('Mouse %d processed in %.2f seconds\n', mouseNum, toc);
    end

    %% ========================== Generate Plots ==========================
    % 1. Boxplot for correlation/explained variance
    generateCorrelationBoxplots(groupResults, currentGroup, PCA_PARAMS, FIGURE_OUTPUT_PATH);
    
    % 2. Boxplot for regional brain signals
    generateRegionalSignalBoxplots(groupResults, currentGroup, PCA_PARAMS, ANALYSIS_TYPES, ...
        bilateralRegions2, FIGURE_OUTPUT_PATH);
    
    % 3. Activity maps (PCA Mean)
    generateActivityMaps(groupResults, currentGroup, PCA_PARAMS, ANALYSIS_TYPES, ...
        CBF_FILE_PATH, mouseNum, FIGURE_OUTPUT_PATH);
    
    % 4. Circulation time series plot + histogram
    generateTimeSeriesPlots(groupResults, currentGroup, SAMPLING_FREQ, FIGURE_OUTPUT_PATH);
end

fprintf('All analysis completed successfully!\n');

%% ========================== Helper Functions ==========================
function processedMasks = preprocessROIMasks(brainRegions)
    % Preprocess ROI masks: remove empty rows + resize to 0.5 scale
    % Input: brainRegions - cell array of ROI names
    % Output: processedMasks - struct of resized masks
    
    processedMasks = struct();
    for roiIdx = 1:length(brainRegions)
        roiName = brainRegions{roiIdx};
        mask = eval(roiName); % Load mask from workspace
        mask(1:209, :) = [];  % Remove empty rows
        mask = double(imresize(mask, 0.5)); % Resize
        mask(mask == 0) = NaN; % Set background to NaN
        processedMasks.(roiName) = mask;
    end
end

function bilateralMasks = combineBilateralRegions(processedMasks, regions1, regions2)
    % Combine left/right hemispheres for bilateral regions
    bilateralMasks = struct();
    
    % Process first set of bilateral regions
    for roiIdx = 1:length(regions1)
        roiName = regions1{roiIdx};
        leftMask = processedMasks.(['L_', roiName]);
        rightMask = processedMasks.(['R_', roiName]);
        bilateralMasks.(roiName) = leftMask + rightMask;
        bilateralMasks.(roiName)(bilateralMasks.(roiName) == 0) = NaN;
    end
    
    % Process second set of bilateral regions
    for roiIdx = 1:length(regions2)
        roiName = regions2{roiIdx};
        leftMask = processedMasks.(['L_', roiName]);
        rightMask = processedMasks.(['R_', roiName]);
        bilateralMasks.(roiName) = leftMask + rightMask;
        bilateralMasks.(roiName)(bilateralMasks.(roiName) == 0) = NaN;
    end
end

function groupResults = initGroupResults(groupName, paramList, typeList, mouseCount, regions1, regions2)
    % Initialize result structure for experimental group
    groupResults = struct();
    groupResults.VesselsPercentage = zeros(mouseCount, 3); % Artery/Vein/Mixed
    
    for paramIdx = 1:length(paramList)
        paramName = paramList{paramIdx};
        groupResults.(paramName) = struct();
        
        % Initialize PCA feature fields
        for typeIdx = 1:length(typeList)
            typeName = typeList{typeIdx};
            % PCA components (top 3)
            groupResults.(paramName).(typeName).PC1 = zeros(mouseCount, numel(processedMasks.Whole_Brain));
            groupResults.(paramName).(typeName).PC2 = zeros(mouseCount, numel(processedMasks.Whole_Brain));
            groupResults.(paramName).(typeName).PC3 = zeros(mouseCount, numel(processedMasks.Whole_Brain));
            % Mean/Variance/Explained variance
            groupResults.(paramName).(typeName).Mean = zeros(mouseCount, size(processedMasks.Whole_Brain,1), size(processedMasks.Whole_Brain,2));
            groupResults.(paramName).(typeName).Variance = zeros(mouseCount, size(processedMasks.Whole_Brain,1), size(processedMasks.Whole_Brain,2));
            groupResults.(paramName).(typeName).Eexplained = zeros(mouseCount, 1);
            % Regional signals
            groupResults.(paramName).(typeName).Regions1 = zeros(mouseCount, length(regions1));
            groupResults.(paramName).(typeName).Regions2 = zeros(mouseCount, length(regions2));
        end
        
        % Correlation fields
        groupResults.(paramName).PCsCorrelation = zeros(mouseCount, 3);
        groupResults.(paramName).MeanCorrelation = zeros(mouseCount, 3);
        groupResults.(paramName).Eexplained = zeros(mouseCount, 3);
        groupResults.(paramName).VascularPercentage = zeros(mouseCount, 3);
    end
end

function [pcaData, timeSeriesData] = loadMouseData(pcaPath, groupName, mouseNum)
    % Load PCA results and time series data for a single mouse
    pcaData = load(fullfile(pcaPath, ['Mice', num2str(mouseNum), '_', groupName, 'CirculationParas.mat']));
    timeSeriesData = load(fullfile(pcaPath, [groupName, num2str(mouseNum), 'CirculationTimeSeries.mat']));
end

function vascularPercent = calculateVascularPercentage(timeSeriesData, wholeBrainMask)
    % Calculate percentage of artery/vein/mixed vessels in whole brain
    arteryMask = timeSeriesData.([groupName, num2str(mouseNum), 'ArteryMask']);
    veinMask = timeSeriesData.([groupName, num2str(mouseNum), 'VeinMask']);
    mixedMask = timeSeriesData.([groupName, num2str(mouseNum), 'MixedMask']);
    
    totalBrainPixels = sum(wholeBrainMask(:) == 1);
    arteryPercent = sum(arteryMask(:) == 1) / totalBrainPixels;
    veinPercent = sum(veinMask(:) == 1) / totalBrainPixels;
    mixedPercent = sum(mixedMask(:) == 1) / totalBrainPixels;
    
    vascularPercent = [arteryPercent, veinPercent, mixedPercent];
end

function groupResults = extractPCAFEatures(groupResults, paramName, pcaData, typeList, mouseIdx, ...
    roiMasks, bilateralMasks, regions1, regions2)
    % Extract PCA features (PCs, Mean, Variance, regional signals)
    
    for typeIdx = 1:length(typeList)
        typeName = typeList{typeIdx};
        % Extract top 3 PCs
        for pcIdx = 1:3
            pcData = pcaData.(typeName).PCs(:,:,pcIdx);
            groupResults.(paramName).(typeName).(['PC', num2str(pcIdx)])(mouseIdx, :) = pcData(:)';
        end
        
        % Mean/Variance/Explained variance
        groupResults.(paramName).(typeName).Mean(mouseIdx, :, :) = pcaData.(typeName).Mean(:,:,1);
        groupResults.(paramName).(typeName).Variance(mouseIdx, :, :) = pcaData.(typeName).Variance;
        groupResults.(paramName).(typeName).Eexplained(mouseIdx) = pcaData.(typeName).Eexplained(1);
        
        % Extract regional signals (bilateral regions 1)
        for roiIdx = 1:length(regions1)
            roiName = regions1{roiIdx};
            roiMask = bilateralMasks.(roiName);
            pc1Data = pcaData.(typeName).PCs(:,:,1) .* roiMask;
            regionalSignal = nanmean(nanmean(pc1Data));
            groupResults.(paramName).(typeName).Regions1(mouseIdx, roiIdx) = regionalSignal;
        end
        
        % Extract regional signals (bilateral regions 2)
        for roiIdx = 1:length(regions2)
            roiName = regions2{roiIdx};
            roiMask = bilateralMasks.(roiName);
            pc1Data = pcaData.(typeName).PCs(:,:,1) .* roiMask;
            regionalSignal = nanmean(nanmean(pc1Data));
            groupResults.(paramName).(typeName).Regions2(mouseIdx, roiIdx) = regionalSignal;
        end
    end
    return groupResults;
end

function groupResults = calculateTypeCorrelations(groupResults, paramName, mouseIdx, typeList)
    % Calculate correlations between General/Top/Bottom analysis types
    
    % Extract PC1 data for correlation
    generalPC1 = groupResults.(paramName).General.PC1(mouseIdx, :);
    topPC1 = groupResults.(paramName).Top.PC1(mouseIdx, :);
    bottomPC1 = groupResults.(paramName).Bottom.PC1(mouseIdx, :);
    
    % Calculate PC correlations
    groupResults.(paramName).PCsCorrelation(mouseIdx, :) = [
        corr(generalPC1', topPC1'), ...
        corr(generalPC1', bottomPC1'), ...
        corr(topPC1', bottomPC1')
    ];
    
    % Extract Mean data for correlation
    generalMean = groupResults.(paramName).General.Mean(mouseIdx, :, :);
    topMean = groupResults.(paramName).Top.Mean(mouseIdx, :, :);
    bottomMean = groupResults.(paramName).Bottom.Mean(mouseIdx, :, :);
    
    % Calculate Mean correlations
    groupResults.(paramName).MeanCorrelation(mouseIdx, :) = [
        corr(generalMean(:), topMean(:)), ...
        corr(generalMean(:), bottomMean(:)), ...
        corr(topMean(:), bottomMean(:))
    ];
    
    % Combine explained variance
    groupResults.(paramName).Eexplained(mouseIdx, :) = [
        groupResults.(paramName).General.Eexplained(mouseIdx), ...
        groupResults.(paramName).Top.Eexplained(mouseIdx), ...
        groupResults.(paramName).Bottom.Eexplained(mouseIdx)
    ];
    return groupResults;
end

function generateCorrelationBoxplots(groupResults, groupName, paramList, figPath)
    % Generate boxplots for PCsCorrelation/MeanCorrelation/Eexplained
    plotParams = {'PCsCorrelation', 'MeanCorrelation', 'Eexplained'};
    
    for paramIdx = 1:length(paramList)
        paramName = paramList{paramIdx};
        for plotIdx = 1:length(plotParams)
            plotParam = plotParams{plotIdx};
            plotData = groupResults.(paramName).(plotParam);
            plotDataCell = {plotData(:,1), plotData(:,2), plotData(:,3)};
            
            % Generate boxplot with scatter
            boxplot_scatter(plotDataCell, {'General','Top','Bottom'});
            saveas(gcf, fullfile(figPath, [groupName, paramName, plotParam, 'Comparison.png']), 'png');
            close all;
        end
    end
end

function generateRegionalSignalBoxplots(groupResults, groupName, paramList, typeList, regions2, figPath)
    % Generate boxplots for regional brain signals (Regions2)
    for paramIdx = 1:length(paramList)
        paramName = paramList{paramIdx};
        for typeIdx = 1:length(typeList)
            typeName = typeList{typeIdx};
            plotData = groupResults.(paramName).(typeName).Regions2;
            plotDataCell = num2cell(plotData, 1);
            
            boxplot_scatter2(plotDataCell, regions2);
            saveas(gcf, fullfile(figPath, [groupName, paramName, typeName, 'RegionalSignals.png']), 'png');
            close all;
        end
    end
end

function generateActivityMaps(groupResults, groupName, paramList, typeList, cbfPath, mouseNum, figPath)
    % Generate smoothed activity maps for PCA Mean values
    addpath('E:\CBF-data\GongData\GongScript\PurifySignals\BasedPCAAnalysisMask\');
    
    % Load CBF mask
    load(fullfile(cbfPath, [groupName, num2str(mouseNum), 'CBFimg.mat']));
    cbfMask = CBFMeanImg(:,:,1) > 0;
    cbfMask = FiltMask(cbfMask, 600, 4); % Filter mask
    cbfMask(cbfMask == 0) = NaN;
    
    for paramIdx = 1:length(paramList)
        paramName = paramList{paramIdx};
        for typeIdx = 1:length(typeList)
            typeName = typeList{typeIdx};
            % Extract Mean data and apply mask
            meanData = groupResults.(paramName).(typeName).Mean(mouseNum, :, :) .* cbfMask;
            meanData = imgaussfilt(meanData, GAUSS_FILTER_SIGMA); % Smooth with Gaussian filter
            
            % Plot activity map
            figure('Position', [1089,383,633,693]);
            imagesc(meanData);
            axis image;
            colormap jet;
            set(gca, 'XTick', [], 'YTick', [], 'Visible', 'off');
            % Remove white edges
            set(gca,'Position',get(gca,'OuterPosition')-get(gca,'TightInset').*[0 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 1]);
            
            saveas(gcf, fullfile(figPath, [groupName, paramName, typeName, 'MeanActivity.png']), 'png');
            close all;
        end
    end
end

function generateTimeSeriesPlots(groupResults, groupName, fs, figPath)
    % Generate circulation time series plot and histogram
    
    % Time vector
    timeVec = (1:length(groupResults.CirculationSignal{2}))/fs;
    circSignal = groupResults.CirculationSignal{1};
    
    % 1. Time series plot
    figure('Position', [307,295,1459,639]);
    hold on;
    plot(timeVec, circSignal, 'Color', [0.3 0.3 0.3], 'LineWidth', 2.8);
    xlim([0 600]);
    set(gca, 'Visible', 'off');
    set(gca,'Position',get(gca,'OuterPosition')-get(gca,'TightInset').*[0 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 1]);
    saveas(gcf, fullfile(figPath, [groupName, 'CirculationTimeSeries.png']), 'png');
    close all;
    
    % 2. Histogram (probability density)
    [counts, binCenters, bins] = histcounts(circSignal, 500);
    binCenters(1) = [];
    binWidth = bins(2) - bins(1);
    pdfVals = abs(counts / (sum(counts) * binWidth));
    
    figure('Position', [417,218,656,639]);
    area(binCenters, pdfVals, 'FaceColor', [0.7 0.9 0.9], 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 3.0);
    xlim([-0.0001 0.003]);
    ylim([0 0.014]);
    set(gca, 'Visible', 'off');
    set(gca,'Position',get(gca,'OuterPosition')-get(gca,'TightInset').*[0 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 1]);
    saveas(gcf, fullfile(figPath, [groupName, 'CirculationHistogram.png']), 'png');
    close all;
end