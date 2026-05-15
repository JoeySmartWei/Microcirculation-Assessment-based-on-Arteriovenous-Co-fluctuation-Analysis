% Map_VesselsPurifiedMixedFeature.m
% Purpose: Process CBF (Cerebral Blood Flow) signals for Stroke/Control mouse groups,
%          extract purified/mixed vascular signals (Artery/Tissue/Vein) and plot time series/phase difference
% Author: Bochao Niu

clearvars; close all; clc;

%% ====================== Configuration (User Modifiable) ======================
% Group definition (Stroke/Control)
groups = {'Stroke','Control'};
% Mouse number sequence for each group (row 1:Stroke, row 2:Control)
mouseSequences = {[1:19]; [1:12,14:31]};
% Vascular ROI types
vesselTypes = {'Artery','Tissue','Vein','ArteryTissue','VeinTissue'};
% Signal analysis types
analysisTypes = {'Correlation','PhaseDiff'};

% File path configuration
mainFilePath = 'E:\CBF-data\GongData\OutPut\CBFPackage\';
roisFilePath = 'E:\CBF-data\GongData\OutPut\VesselsROIs\';
outputPath = 'E:\CBF-data\GongData\OutPut\VesselsROISignal\Figures\';

% Constant parameters
numROIs = 20;                % Number of ROIs to process
plotDPI = 300;               % Resolution for saved figures
purifyROIRange = 11:20;      % ROI range for purified signal plotting

%% ====================== Main Processing Loop ======================
% Process Control group (index 2 in groups, adjust if need to process Stroke)
for groupIdx = 2:length(groups)
    currentGroup = groups{groupIdx};
    currentMiceSeq = mouseSequences{groupIdx};
    
    % Only process mouse 10 (original hardcoded sq=10, adjust to 1:length(currentMiceSeq) for full sequence)
    mouseIdx = 10;             
    if mouseIdx > length(currentMiceSeq)
        warning('Mouse index %d out of range for %s group', mouseIdx, currentGroup);
        continue;
    end
    
    mouseNum = currentMiceSeq(mouseIdx);
    fprintf('Processing %s Group: Mouse %d CBF Image ......\n', currentGroup, mouseNum);
    
    % Create output directory if not exist
    createDirectoryIfNotExist(outputPath);
    
    %% Load raw data
    % Load CBF signal and vascular ROIs
    cbfSignal = loadCBFSignal(mainFilePath, currentGroup, mouseNum);
    vesselROIs = loadVesselROIs(roisFilePath, currentGroup, mouseNum);
    
    %% Preprocess CBF signal
    % Reshape: [width x length x time] -> [time x (width*length)]
    [width, length, timeLength] = size(cbfSignal);
    cbfSignal = reshape(cbfSignal, width*length, timeLength)';
    % Preprocess: detrend -> fill outliers -> z-score normalization
    cbfSignal = preprocessCBFSignal(cbfSignal);
    
    %% Extract ROI signals for each vessel type
    % Get index count from Artery ROI (reference for all vessel types)
    arteryROIFirst = vesselROIs.Artery(:,:,1);
    indexCount = length(find(arteryROIFirst == 1));
    vesselSignals = extractVesselSignals(vesselROIs, cbfSignal, vesselTypes, timeLength, indexCount, numROIs);
    
    % Clear large raw signal variable to save memory
    clear cbfSignal;
    
    %% Plot Purified Signals (Time Series + Phase Difference)
    plotPurifiedSignals(vesselSignals, vesselTypes, purifyROIRange, outputPath, mouseNum, currentGroup, plotDPI);
    
    %% Plot Mixed Signals (Artery-Tissue/Artery-Vein/Tissue-Vein combinations)
    plotMixedSignals(vesselSignals, outputPath, mouseNum, currentGroup, plotDPI);
end

%% ====================== Helper Functions ======================
function createDirectoryIfNotExist(dirPath)
    % Create directory if it does not exist
    % Input: dirPath - Full path of target directory
    if ~exist(dirPath, 'dir')
        mkdir(dirPath);
        fprintf('Created directory: %s\n', dirPath);
    end
end

function cbfSignal = loadCBFSignal(filePath, group, mouseNum)
    % Load CBF signal mat file
    fileName = [filePath, group, num2str(mouseNum), 'CBFSignal.mat'];
    data = load(fileName);
    varName = [group, num2str(mouseNum), 'CBFSignal'];
    cbfSignal = data.(varName);
end

function vesselROIs = loadVesselROIs(roisPath, group, mouseNum)
    % Load vascular ROIs and organize into struct
    % Output: vesselROIs - Struct with fields: Artery, Tissue, Vein, ArteryTissue, VeinTissue
    fileName = [roisPath, group, 'Mice', num2str(mouseNum), 'VascularROIs.mat'];
    data = load(fileName);
    vesselTypes = {'Artery','Tissue','Vein','ArteryTissue','VeinTissue'};
    
    vesselROIs = struct();
    for i = 1:length(vesselTypes)
        roiVarName = ['Mice', num2str(mouseNum), '_', vesselTypes{i}, 'ROIs'];
        vesselROIs.(vesselTypes{i}) = data.(roiVarName);
    end
end

function processedSignal = preprocessCBFSignal(rawSignal)
    % Preprocess CBF signal: detrend -> fill outliers -> z-score
    % Input: rawSignal - Raw CBF signal matrix [time x pixels]
    % Output: processedSignal - Normalized signal
    signalDouble = double(rawSignal);
    signalDetrend = detrend(signalDouble);
    signalFilled = filloutliers(signalDetrend, 'nearest', 'mean');
    processedSignal = zscore(signalFilled);
end

function vesselSignals = extractVesselSignals(vesselROIs, cbfSignal, vesselTypes, timeLen, idxCount, numROIs)
    % Extract signal for each vessel ROI type
    % Output: vesselSignals - Struct with vessel type as fields, each is [time x index x ROI]
    
    vesselSignals = struct();
    for i = 1:length(vesselTypes)
        vesselType = vesselTypes{i};
        roiData = vesselROIs.(vesselType);
        % Initialize signal matrix: [time x index count x number of ROIs]
        vesselSignals.(vesselType) = zeros(timeLen, idxCount, numROIs);
        
        for roiNum = 1:numROIs
            % Find ROI pixels (value = 1)
            roiPixels = find(roiData(:,:,roiNum) == 1);
            % Extract signal for current ROI
            vesselSignals.(vesselType)(:,:,roiNum) = cbfSignal(:, roiPixels);
        end
    end
end

function plotPurifiedSignals(vesselSignals, vesselTypes, roiRange, outPath, mouseNum, group, dpi)
    % Plot purified signals (Artery/Tissue/Vein) and save figures
    for roiNum = roiRange
        % Combine Artery/Tissue/Vein signals (first 3 vessel types)
        pixelTimeseriesNew = [];
        for i = 1:3
            vesselType = vesselTypes{i};
            pixelData = vesselSignals.(vesselType)(:,:,roiNum);
            pixelTimeseriesNew = cat(3, pixelTimeseriesNew, pixelData);
        end
        
        close all;
        % Plot phase-based purified signal
        [selectSignal, phaseDiff] = PurifyPixelSignalPhase(pixelTimeseriesNew);
        savePurifiedFigs(outPath, mouseNum, group, roiNum, dpi);
        
        % Plot single-channel purified signal
        [selectSignal, phaseDiff] = PurifyPixelSignal(pixelTimeseriesNew);
        print(figure(10), [outPath, 'Mice', num2str(mouseNum), '_', group, 'TimeseriesSingle'], '-dpng', ['-r', num2str(dpi)]);
    end
end

function savePurifiedFigs(outPath, mouseNum, group, roiNum, dpi)
    % Save purified signal figures (helper for plotPurifiedSignals)
    baseName = [outPath, 'Mice', num2str(mouseNum), '_', group];
    print(figure(10), [baseName, 'Timeseries'], '-dpng', ['-r', num2str(dpi)]);
    print(figure(12), [baseName, 'ArteryPhaseDifferenceNumber', num2str(roiNum)], '-dpng', ['-r', num2str(dpi)]);
    print(figure(13), [baseName, 'TissuePhaseDifferenceNumber', num2str(roiNum)], '-dpng', ['-r', num2str(dpi)]);
    print(figure(14), [baseName, 'VeinPhaseDifferenceNumber', num2str(roiNum)], '-dpng', ['-r', num2str(dpi)]);
end

function plotMixedSignals(vesselSignals, outPath, mouseNum, group, dpi)
    % Plot mixed signals (mean of ROI signals) and save figures
    % Calculate mean across ROIs for each vessel type
    arteryMean = mean(vesselSignals.Artery, 3);
    tissueMean = mean(vesselSignals.Tissue, 3);
    veinMean = mean(vesselSignals.Vein, 3);
    
    % Combine mixed signal combinations
    newSignalA = cat(3, arteryMean, tissueMean);  % Artery-Tissue
    newSignalB = cat(3, arteryMean, veinMean);    % Artery-Vein
    newSignalC = cat(3, tissueMean, veinMean);    % Tissue-Vein
    
    close all;
    % Plot phase-based mixed signal
    [mixedSignalParas] = MixPixelSignalPhase(newSignalA, newSignalB, newSignalC);
    saveMixedFigs(outPath, mouseNum, group, dpi);
    
    % Plot single-channel mixed signal
    [mixedSignalParas] = MixPixelSignal(newSignalA, newSignalB, newSignalC);
    print(figure(10), [outPath, 'Mice', num2str(mouseNum), '_', group, 'TimeseriesMixedSingle'], '-dpng', ['-r', num2str(dpi)]);
end

function saveMixedFigs(outPath, mouseNum, group, dpi)
    % Save mixed signal figures (helper for plotMixedSignals)
    baseName = [outPath, 'Mice', num2str(mouseNum), '_', group];
    print(figure(10), [baseName, 'TimeseriesMixed'], '-dpng', ['-r', num2str(dpi)]);
    print(figure(12), [baseName, 'MixedArteryTissuePhaseDifference'], '-dpng', ['-r', num2str(dpi)]);
    print(figure(13), [baseName, 'MixedArteryVeinPhaseDifference'], '-dpng', ['-r', num2str(dpi)]);
    print(figure(14), [baseName, 'MixedTissueVeinPhaseDifference'], '-dpng', ['-r', num2str(dpi)]);
end