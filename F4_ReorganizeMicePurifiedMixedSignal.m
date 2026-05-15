% ReorganizeMicePurifiedMixedSignal.m
% Purpose: Process CBF (Cerebral Blood Flow) signals for stroke/control mouse groups
%          Extract artery/vein signals, calculate correlation/phase difference matrices
% Author: Bochao Niu

clear all; close all; clc;

%% ====================== Configuration Parameters ======================
% Group settings: 'Stroke' and 'Control' groups
groupList = {'Stroke','Control'};
% Sequence of mouse numbers for each group (Control: 1-12,14-31; Stroke:1-19)
mouseSequence = {[1:19], [1:12,14:31]};  
% Vascular types
vesselTypes = {'Artery','Vein'};

% File path configuration (MODIFY THESE PATHS ACCORDING TO YOUR ENVIRONMENT)
filePath = 'F:\CBF-data\GongData\OutPut\CBFPackage\';
pcaPath = 'F:\CBF-data\GongData\OutPut\PCAResults\';
outputPath = 'F:\CBF-data\GongData\OutPut\PCAVascularSignal\';

%% ====================== Main Processing Loop ======================
% Create output directory if not exists
if ~exist(outputPath, 'dir')
    mkdir(outputPath);
    fprintf('Created output directory: %s\n', outputPath);
end

% Process each group (start from index 2 to skip header)
for groupIdx = 2:length(groupList)
    currentGroup = groupList{groupIdx};
    currentMouseSeq = mouseSequence{groupIdx};
    
    % Process each mouse in the group (start from index 2 to skip header)
    for seqIdx = 2:length(currentMouseSeq)
        tic; % Start timer for processing time calculation
        mouseNum = currentMouseSeq(seqIdx);
        fprintf('\n=============================================\n');
        fprintf('Processing %s Group: Mouse %d CBF Image ......\n', currentGroup, mouseNum);
        fprintf('Start time: %s\n', datestr(clock));

        %% 1. Load input data
        % Load CBF signal data
        cbfSignalFile = [filePath, currentGroup, num2str(mouseNum), 'CBFSignal.mat'];
        pcaResultFile = [pcaPath, currentGroup, num2str(mouseNum), 'PCACirculation.mat'];
        
        % Check if files exist
        if ~exist(cbfSignalFile, 'file')
            error('CBF Signal file not found: %s', cbfSignalFile);
        end
        if ~exist(pcaResultFile, 'file')
            error('PCA Result file not found: %s', pcaResultFile);
        end
        
        load(cbfSignalFile);
        load(pcaResultFile);
        
        % Get CBF signal and PCA results (avoid eval by direct variable assignment)
        cbfSignalVarName = [currentGroup, num2str(mouseNum), 'CBFSignal'];
        pcaResultVarName = [currentGroup, num2str(mouseNum), 'PCACirculation'];
        CBFSignal = evalin('base', cbfSignalVarName);
        PCAResults = evalin('base', pcaResultVarName);

        %% 2. Extract artery/vein signal masks
        % Get PCA time series and sort by correlation sum
        pcaTimeSeries = PCAResults.RawClusterSignal;
        [~, sortedIdx] = sort(sum(corr(pcaTimeSeries), 2));
        
        % Find vein/artery indices from PCA cluster map
        veinIndex = find(PCAResults.RawClusterMap == sortedIdx(1));
        arteryIndex = find(PCAResults.RawClusterMap == sortedIdx(2));

        %% 3. Preprocess CBF signal
        % Clear redundant variables to save memory
        clearvars -regexp ^Control\d+ ^Stroke\d+ ^Mice\d+;
        
        % Reshape CBF signal: [width x length x time] -> [time x (width*length)]
        [width, length, timeLen] = size(CBFSignal);
        CBFSignal = transpose(reshape(CBFSignal, [width*length, timeLen]));
        
        % Signal normalization: detrend -> fill outliers -> z-score
        CBFSignal = double(CBFSignal);
        CBFSignal = detrend(CBFSignal);
        CBFSignal = filloutliers(CBFSignal, 'nearest', 'mean');
        CBFSignal = zscore(CBFSignal);

        %% 4. Extract artery/vein signals (sample first 5000 points)
        sampleSize = 5000;
        arterySignal = CBFSignal(:, arteryIndex(1:min(sampleSize, length(arteryIndex))));
        veinSignal = CBFSignal(:, veinIndex(1:min(sampleSize, length(veinIndex))));
        
        % Combine artery and vein signals into 3D matrix
        mixedSignal = cat(3, arterySignal, veinSignal);

        %% 5. Calculate correlation and phase difference matrices
        [mixedCorrMatrix, mixedPhaseDiffMatrix] = MixBigPixelSignalPackage(mixedSignal);

        %% 6. Save results (use structured variable naming)
        % Define result variable names
        corrMatrixVar = ['Mice', num2str(mouseNum), '_', currentGroup, 'MixedCorrMatrix'];
        phaseDiffMatrixVar = ['Mice', num2str(mouseNum), '_', currentGroup, 'MixedPhaseDiffMatrix'];
        
        % Assign variables and save
        assignin('base', corrMatrixVar, mixedCorrMatrix);
        assignin('base', phaseDiffMatrixVar, mixedPhaseDiffMatrix);
        
        save([outputPath, 'Mice', num2str(mouseNum), '_', currentGroup, 'VascularMixedMatrixParas.mat'], ...
            corrMatrixVar, phaseDiffMatrixVar);

        %% 7. Cleanup and log
        % Clear temporary variables to save memory
        clearvars CBFSignal mixedSignal arterySignal veinSignal;
        clearvars -regexp ^Mice\d+ MixedPixelParas MixedROIParas;
        
        % Log processing time and completion
        processingTime = toc;
        fprintf('Completed processing Mouse %d (Group: %s)\n', mouseNum, currentGroup);
        fprintf('End time: %s\n', datestr(clock));
        fprintf('Processing time: %.2f seconds\n', processingTime);
        fprintf('=============================================\n');
    end
end

fprintf('\nAll processing completed!\n');

