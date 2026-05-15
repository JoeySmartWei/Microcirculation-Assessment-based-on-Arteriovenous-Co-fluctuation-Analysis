%% Calculate Vascular Purification and Mixing Properties for Mice CBF Data
% Author: Anonymous
% Function: Compute multi-dimensional vascular properties (Corr/PhaseDiff/Variance etc.)
% Input: CBF signal mat files, Vascular ROIs mat files
% Output: Saved mat files with purified/mixed vessel signal properties

% ========================== Configuration Parameters ==========================
clear all; close all; clc;

% Experimental groups (Stroke/Control)
expGroups = {'Stroke','Control'};
% Mouse number sequences for each group (Stroke:[1:2], Control:[1:12,14:31])
mouseSequences = {[1:19]; [1:12,14:31]};
% Vascular ROI types
vesselTypes = {'Artery','Tissue','Vein','ArteryTissue','VeinTissue'};
% Analysis types (extended properties)
analysisTypes = {'Corr','PhaseDiff','Variance','InnerVariance','InnerPhaseDiff'};
% ROI combination matrix (all pairwise combinations of 5 vessel types)
vesselCombinations = [1 2; 1 3; 1 4; 1 5; 2 3; 2 4; 2 5; 3 4; 3 5; 4 5];

% File path configuration (modify these paths according to your environment)
config.filePaths = struct(...
    'CBFSignal', 'F:\CBF-data\GongData\OutPut\CBFPackage\', ...
    'ROIs',      'F:\CBF-data\GongData\OutPut\VesselsROIs\', ...
    'Output',    'F:\CBF-data\GongData\OutPut\VesselsROISignal\');

% ROI constants
config.ROINums = 20;  % Total number of ROIs per vessel type
config.startMouseIdx = 21;  % Start index for mouse sequence (Control group)

% ========================== Main Processing Loop ==========================
% Process only Control group (index 2) - keep consistent with original logic
for groupIdx = 2  % Fixed to Control group (index 2)
    currentGroup = expGroups{groupIdx};
    currentMiceSeq = mouseSequences{groupIdx};
    
    % Create output directory if not exists
    if ~exist(config.filePaths.Output, 'dir')
        mkdir(config.filePaths.Output);
        fprintf('Created output directory: %s\n', config.filePaths.Output);
    end
    
    tic; % Start timer for group processing
    fprintf('\n=============================================\n');
    fprintf('Starting processing for %s Group (Mice %d to %d)\n', ...
        currentGroup, currentMiceSeq(config.startMouseIdx), currentMiceSeq(end));
    
    % Process mice from start index to end
    for mouseIdx = config.startMouseIdx:length(currentMiceSeq)
        mouseNum = currentMiceSeq(mouseIdx);
        fprintf('\nProcessing %s Group: Mouse %d CBF Image ......\n', currentGroup, mouseNum);
        
        % ========================== Load Data ==========================
        % Load CBF signal and vascular ROIs
        CBFSignalFile = fullfile(config.filePaths.CBFSignal, [currentGroup, num2str(mouseNum), 'CBFSignal.mat']);
        roisFile = fullfile(config.filePaths.ROIs, [currentGroup, 'Mice', num2str(mouseNum), 'VascularROIs.mat']);
        
        % Check file existence
        if ~exist(CBFSignalFile, 'file') || ~exist(roisFile, 'file')
            warning('Missing files for Mouse %d - skipping!', mouseNum);
            continue;
        end
        
        load(CBFSignalFile);
        load(roisFile);
        
        % Extract CBF signal and ROIs (replace eval with direct access)
        cbfSignalVarName = [currentGroup, num2str(mouseNum), 'CBFSignal'];
        CBFSignal = evalin('base', cbfSignalVarName);
        
        % Store ROIs in struct for easy access
        vesselROIs = struct();
        for vesIdx = 1:length(vesselTypes)
            roiVarName = ['Mice', num2str(mouseNum), '_', vesselTypes{vesIdx}, 'ROIs'];
            vesselROIs.(vesselTypes{vesIdx}) = evalin('base', roiVarName);
        end
        
        % Clear raw data variables to save memory
        clearVarsByPattern([ '^', currentGroup, num2str(mouseNum) ], [ '^Mice', num2str(mouseNum) ]);
        
        % ========================== Preprocess CBF Signal ==========================
        [width, length, timePoints] = size(CBFSignal);
        % Reshape: [width*length, timePoints] -> transpose to [timePoints, width*length]
        CBFSignal = transpose(reshape(CBFSignal, [width*length, timePoints]));
        % Preprocess: fill outliers -> detrend -> z-score normalization
        CBFSignal = zscore(detrend(filloutliers(double(CBFSignal), 'nearest', 'mean')));
        
        % ========================== Extract Vessel Signals ==========================
        % Get pixel count from Artery ROI (reference for all vessel types)
        refPixelCount = length(find(vesselROIs.Artery(:,:,1) == 1));
        vesselSignals = struct();  % Store signals for each vessel type
        
        for vesIdx = 1:length(vesselTypes)
            vesselName = vesselTypes{vesIdx};
            % Initialize signal matrix: [timePoints, pixelCount, ROI count]
            vesselSignals.(vesselName) = zeros(timePoints, refPixelCount, config.ROINums);
            
            % Extract signal for each ROI
            for roiNum = 1:config.ROINums
                % Find pixel indices for current ROI
                roiPixelIdx = find(vesselROIs.(vesselName)(:,:,roiNum) == 1);
                % Extract corresponding CBF signals
                vesselSignals.(vesselName)(:, :, roiNum) = CBFSignal(:, roiPixelIdx);
            end
        end
        
        % ========================== Purified Signal Analysis ==========================
        pixelsSignal = struct();   % Pixel-level purified properties
        roisSignal = struct();     % ROI-level purified properties
        
        for vesIdx = 1:length(vesselTypes)
            vesselName = vesselTypes{vesIdx};
            vesselSignal = vesselSignals.(vesselName);
            
            % Compute mean timeseries (pixel: average over ROIs; ROI: average over pixels)
            pixelTimeseries = mean(vesselSignal, 3);  % Average across ROIs
            roiTimeseries = squeeze(mean(vesselSignal, 2));  % Average across pixels
            
            % Purify signal and compute all properties
            [pixelCorr, pixelPhaseDiff, pixelVariance, pixelInnerVariance, pixelInnerPhaseDiff] = ...
                PurifyPixelSignal(pixelTimeseries);
            [roiCorr, roiPhaseDiff, roiVariance, roiInnerVariance, roiInnerPhaseDiff] = ...
                PurifyPixelSignal(roiTimeseries);
            
            % Store results in struct (dynamic field names)
            pixelsSignal.([vesselName, '_Corr']) = pixelCorr;
            pixelsSignal.([vesselName, '_PhaseDiff']) = pixelPhaseDiff;
            pixelsSignal.([vesselName, '_Variance']) = pixelVariance;
            pixelsSignal.([vesselName, '_InnerVariance']) = pixelInnerVariance;
            pixelsSignal.([vesselName, '_InnerPhaseDiff']) = pixelInnerPhaseDiff;
            
            roisSignal.([vesselName, '_Corr']) = roiCorr;
            roisSignal.([vesselName, '_PhaseDiff']) = roiPhaseDiff;
            roisSignal.([vesselName, '_Variance']) = roiVariance;
            roisSignal.([vesselName, '_InnerVariance']) = roiInnerVariance;
            roisSignal.([vesselName, '_InnerPhaseDiff']) = roiInnerPhaseDiff;
        end
        
        % ========================== Mixed Vessel Signal Analysis ==========================
        mixedPixelParams = struct();  % Mixed pixel-level parameters
        mixedROIParams = struct();    % Mixed ROI-level parameters
        
        for combIdx = 1:size(vesselCombinations, 1)
            % Get paired vessel types
            ves1Idx = vesselCombinations(combIdx, 1);
            ves2Idx = vesselCombinations(combIdx, 2);
            ves1Name = vesselTypes{ves1Idx};
            ves2Name = vesselTypes{ves2Idx};
            combName = [ves1Name, '2', ves2Name];
            
            % ---------------- Pixel-level mixing ----------------
            fixPixelTS = mean(vesselSignals.(ves1Name), 3);
            movePixelTS = mean(vesselSignals.(ves2Name), 3);
            pixelMixSignal = cat(3, fixPixelTS, movePixelTS);
            mixedPixelParas = MixPixelSignal(pixelMixSignal);
            
            % ---------------- ROI-level mixing ----------------
            roiSignalFix = squeeze(mean(vesselSignals.(ves1Name), 2));
            roiSignalMove = squeeze(mean(vesselSignals.(ves2Name), 2));
            roiMixSignal = cat(3, roiSignalFix, roiSignalMove);
            mixedROIParas = MixPixelSignal(roiMixSignal);
            
            % Store mixed parameters
            mixedPixelParams.([combName, 'Paras']) = mixedPixelParas;
            mixedROIParams.([combName, 'Paras']) = mixedROIParas;
        end
        
        % ========================== Save Results ==========================
        % Create result variable names
        resVarNames = struct(...
            'pixelsSignal', ['Mice', num2str(mouseNum), '_', currentGroup, 'PixelsSignal'], ...
            'roisSignal',   ['Mice', num2str(mouseNum), '_', currentGroup, 'ROIsSignal'], ...
            'mixedPixels',  ['Mice', num2str(mouseNum), '_', currentGroup, 'MixedPixelParas'], ...
            'mixedROIs',    ['Mice', num2str(mouseNum), '_', currentGroup, 'MixedROIParas']);
        
        % Assign results to variables
        assignin('base', resVarNames.pixelsSignal, pixelsSignal);
        assignin('base', resVarNames.roisSignal, roisSignal);
        assignin('base', resVarNames.mixedPixels, mixedPixelParams);
        assignin('base', resVarNames.mixedROIs, mixedROIParams);
        
        % Save to mat file
        savePath = fullfile(config.filePaths.Output, ...
            ['Mice', num2str(mouseNum), '_', currentGroup, 'VascularPixelsROIsSignal.mat']);
        save(savePath, resVarNames.pixelsSignal, resVarNames.roisSignal, ...
            resVarNames.mixedPixels, resVarNames.mixedROIs);
        
        % Clear temporary variables
        clear CBFSignal;
        clearVarsByPattern(['^Mice', num2str(mouseNum)], 'mixedPixelParas', 'mixedROIParas');
        
        fprintf('Completed processing Mouse %d\n', mouseNum);
    end
    
    % ========================== Group Processing Summary ==========================
    fprintf('\n=============================================\n');
    fprintf('Completed processing for %s Group\n', currentGroup);
    fprintf('End time: %s\n', datestr(clock));
    fprintf('Total elapsed time for group: %.2f minutes\n', toc/60);
end

fprintf('\nAll processing completed!\n');

%% Helper Function: Clear Variables by Pattern
function clearVarsByPattern(varargin)
    % CLEARVARSBYPATTERN Clear variables matching regex patterns
    % Input: Variable number of regex patterns
    for i = 1:nargin
        pattern = varargin{i};
        vars = who('-regexp', pattern);
        if ~isempty(vars)
            clear(vars{:});
        end
    end
end