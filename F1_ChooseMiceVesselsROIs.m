% ChooseMiceVesselsROIs.m
% Purpose: Automatically select and save vascular ROIs (Region of Interest) for stroke/control mice CBF images
% Author: [Your Name]
% Date: [Date]
% Version: 1.0

%% Clear workspace and set parameters
clc; clearvars; close all;

% Experimental groups (Stroke/Control)
expGroups = {'Stroke','Control'};
% Sequence of mouse numbers for each group
mouseSequences = {[1:2]; [1:12,14:31]};
% Vascular ROI types to analyze
vesselTypes = {'Artery','Tissue','Vein','ArteryTissue','VeinTissue'};
% Base path for data storage
baseFilePath = 'F:\CBF-data\GongData\';

% ROI configuration
roiNumber = 20; % Number of ROIs to select for each vessel type
imageResizeFactor = 0.5; % Image resize scale factor
gaussianFilterSigma = 1.5; % Sigma for Gaussian filtering
gaussianFilterSize = 5; % Filter size for Gaussian smoothing

%% Process each experimental group
for groupIdx = 2:length(expGroups)
    currentGroup = expGroups{groupIdx};
    currentMouseSeq = mouseSequences{groupIdx,1};
    
    % Process each mouse in the sequence (starting from index 21)
    for seqIdx = 21:length(currentMouseSeq)
        mouseNumber = currentMouseSeq(seqIdx);
        fprintf('Processing %s Group: Mouse %d CBF Image ......\n', currentGroup, mouseNumber);
        
        %% 1. Set up file paths
        % Input path for CBF images
        cbfImagePath = fullfile(baseFilePath, 'CBFImages', [currentGroup, num2str(mouseNumber)]);
        % Output path for ROI results
        outputPath = fullfile(baseFilePath, 'OutPut', 'VesselsROIs');
        
        % Create output directory if not exists
        if ~exist(outputPath, 'dir')
            mkdir(outputPath);
            fprintf('Created output directory: %s\n', outputPath);
        end
        
        %% 2. Load and preprocess CBF images (1001-1200 frames)
        % Get reference image size
        refImagePath = fullfile(cbfImagePath, '1.jpg');
        if ~exist(refImagePath, 'file')
            warning('Reference image not found: %s. Skipping mouse %d', refImagePath, mouseNumber);
            continue;
        end
        
        refImage = imread(refImagePath);
        refImage = imresize(fliplr(permute(refImage, [2 1 3])), imageResizeFactor);
        [imgWidth, imgHeight, imgChannels] = size(refImage);
        
        % Preallocate memory for CBF image sequence (1001-1200 frames)
        cbfImageSequence = uint8(zeros(imgWidth, imgHeight, imgChannels, 200));
        frameIdx = 1;
        
        % Load and process each frame
        for frameNumber = 1001:1200
            framePath = fullfile(cbfImagePath, num2str(frameNumber) + '.jpg');
            if ~exist(framePath, 'file')
                warning('Frame %d not found: %s. Using empty frame', frameNumber, framePath);
                continue;
            end
            
            % Read and preprocess frame
            cbfFrame = imread(framePath);
            % Gaussian filtering for noise reduction and signal smoothing
            cbfFrame = imgaussfilt(cbfFrame, gaussianFilterSigma, 'FilterSize', gaussianFilterSize);
            cbfFrame = imgaussfilt(cbfFrame, gaussianFilterSigma);
            % Resize and flip to standard orientation
            cbfFrame = imresize(fliplr(permute(cbfFrame, [2 1 3])), imageResizeFactor);
            
            cbfImageSequence(:,:,:,frameIdx) = cbfFrame;
            frameIdx = frameIdx + 1;
        end
        
        % Create mean CBF image and resize again
        meanCBFImage = uint8(mean(cbfImageSequence, 4));
        meanCBFImage = imresize(meanCBFImage, imageResizeFactor);
        
        %% 3. Generate ROIs for each vessel type
        roiData = struct(); % Struct to store ROI data (avoids eval)
        
        for vesselIdx = 1:length(vesselTypes)
            currentVessel = vesselTypes{vesselIdx};
            % Generate ROI mask using FormMask function
            [vesselROI] = FormMask(meanCBFImage, currentVessel, roiNumber);
            roiData.(currentVessel) = vesselROI;
        end
        
        %% 4. Visualize and save ROI results
        % Create figure for ROI visualization
        fig = figure('Position', [150, 150, 1578, 967]);
        subplot(2, 3, 1);
        imshow(meanCBFImage);
        set(gca, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', []);
        title('CBF Image (Mean of Frames 1001-1200)', 'FontSize', 20);
        
        % Plot each vessel ROI
        for vesselIdx = 1:length(vesselTypes)
            currentVessel = vesselTypes{vesselIdx};
            roiMask = roiData.(currentVessel);
            roiMaskFinal = roiMask(:,:,end); % Get final ROI mask
            
            subplot(2, 3, vesselIdx + 1);
            imshow(meanCBFImage .* uint8(roiMaskFinal));
            set(gca, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', []);
            title([currentVessel, ' ROIs'], 'FontSize', 20);
        end
        
        % Add super title and save figure
        sgtitle([currentGroup, ' Mice ', num2str(mouseNumber), ' Vascular ROIs Coordinates'], 'FontSize', 30);
        figSavePath = fullfile(outputPath, [currentGroup, 'Mice', num2str(mouseNumber), '_Vascular_ROIs_Coordinates.jpg']);
        saveas(fig, figSavePath);
        fprintf('Saved ROI visualization: %s\n', figSavePath);
        close(fig);
        
        %% 5. Save ROI and CBF image data to MAT files
        % Prepare ROI data for saving
        saveData = struct();
        saveData.CBFImage = meanCBFImage;
        
        for vesselIdx = 1:length(vesselTypes)
            currentVessel = vesselTypes{vesselIdx};
            saveData.([currentVessel, 'ROIs']) = roiData.(currentVessel);
        end
        
        % Save ROI data
        roiMatPath = fullfile(outputPath, [currentGroup, 'Mice', num2str(mouseNumber), 'VascularROIs.mat']);
        save(roiMatPath, '-struct', 'saveData', 'ArteryROIs', 'VeinROIs', 'TissueROIs', ...
            'ArteryTissueROIs', 'VeinTissueROIs');
        
        % Save CBF image data
        cbfMatPath = fullfile(outputPath, [currentGroup, 'Mice', num2str(mouseNumber), 'CBFimg.mat']);
        save(cbfMatPath, '-struct', 'saveData', 'CBFImage');
        
        fprintf('Saved ROI data: %s\n', roiMatPath);
        fprintf('Saved CBF image: %s\n', cbfMatPath);
    end
end

fprintf('All mice processing completed!\n');

%% Helper Function (if FormMask is not a separate function)
% function [mask] = FormMask(image, vesselType, roiNum)
%     % FormMask: Generate ROI mask for specific vessel type
%     % Input:
%     %   image - Input CBF image (mean image)
%     %   vesselType - Type of vessel (Artery/Tissue/Vein etc.)
%     %   roiNum - Number of ROIs to generate
%     % Output:
%     %   mask - Binary ROI mask (same size as input image)
%     % TODO: Implement your ROI selection logic here
% end