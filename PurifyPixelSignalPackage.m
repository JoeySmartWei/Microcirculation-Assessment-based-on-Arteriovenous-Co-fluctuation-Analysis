function [correlationMatrix, phaseDiffMatrix] = PurifyPixelSignalPackage(pixelSignal)
% PURIFYPIXELSIGNALPACKAGE - Extract and analyze pure pixel signals by random sampling
%
% Syntax:
%   [correlationMatrix, phaseDiffMatrix] = PurifyPixelSignalPackage(pixelSignal)
%
% Inputs:
%   pixelSignal - 2D matrix (timeLength × pixelNumber), original pixel signal matrix
%
% Outputs:
%   correlationMatrix - 3D matrix (repeats × repeats × chooseNum), correlation coefficients between sampled signals
%   phaseDiffMatrix   - 3D matrix (repeats × repeats × chooseNum), phase differences between sampled signals
%
% Description:
%   Randomly extract pixel signals, calculate mean/variance/similarity with reference signal,
%   and output correlation and phase difference matrices across different sampling scales.

    % Get basic dimensions of input signal
    [timeLength, pixelNumber] = size(pixelSignal);
    
    % Configuration parameters
    numRepeats = 20;                  % Number of random sampling repetitions
    chooseNum = 1:(pixelNumber - 1);  % Range of pixel numbers to select each time
    refSignal = mean(pixelSignal, 2); % Reference signal (mean of all pixel signals)
    
    % Set sampling gap (optimize computation for large pixel numbers)
    if pixelNumber > 50
        gap = pixelNumber / 100;
    else
        gap = 1;
    end
    
    % Initialize output matrices with NaN (for uncalculated values)
    correlationMatrix = nan(numRepeats, numRepeats, length(chooseNum));
    phaseDiffMatrix = nan(numRepeats, numRepeats, length(chooseNum));
    
    % Main loop: sample signals with different pixel numbers
    for chooseIdx = 1:gap:length(chooseNum)
        currChooseNum = chooseNum(chooseIdx);
        selectSignal = zeros(timeLength, numRepeats);
        
        % Random sampling loop
        for rdmIdx = 1:numRepeats
            % Randomly select pixel indices and calculate mean signal
            selectIndex = randperm(pixelNumber, currChooseNum)';
            selectSignal(:, rdmIdx) = mean(pixelSignal(:, selectIndex), 2);
        end
        
        % Calculate correlation and phase difference
        [corrMat, phaseDiffMat] = calculateAssociation(selectSignal);
        correlationMatrix(:, :, chooseIdx) = corrMat;
        phaseDiffMatrix(:, :, chooseIdx) = phaseDiffMat;
    end

end

% -------------------------------------------------------------------------
function [corrMat, phaseDiffMat] = calculateAssociation(signalMatrix)
% CALCULATEASSOCIATION - Calculate correlation and phase difference between signal pairs
%
% Inputs:
%   signalMatrix - 2D matrix (timeLength × signalNumber), matrix of signals to compare
%
% Outputs:
%   corrMat      - 2D matrix (signalNumber × signalNumber), correlation coefficients
%   phaseDiffMat - 2D matrix (signalNumber × signalNumber), phase differences

    numSignals = size(signalMatrix, 2);
    corrMat = nan(numSignals, numSignals);
    phaseDiffMat = nan(numSignals, numSignals);
    
    % Calculate pairwise correlation and phase difference
    for i = 1:numSignals
        for j = (i + 1):numSignals
            corrMat(j, i) = corr(signalMatrix(:, i), signalMatrix(:, j));
            phaseDiffMat(i, j) = calculatePhaseDifference(signalMatrix(:, i), signalMatrix(:, j));
        end
    end

end

% -------------------------------------------------------------------------
function avgPhase = calculatePhaseDifference(signal1, signal2)
% CALCULATEPHASEDIFfERENCE - Compute average phase difference between two signals using Hilbert transform
%
% Inputs:
%   signal1 - 1D vector, first input signal
%   signal2 - 1D vector, second input signal
%
% Outputs:
%   avgPhase - Scalar, average phase difference (absolute value)

    % Step 1: Hilbert transform to get analytic signals
    hSignal1 = hilbert(signal1);
    hSignal2 = hilbert(signal2);
    
    % Step 2: Unwrap instantaneous phases
    phase1 = unwrap(angle(hSignal1));
    phase2 = unwrap(angle(hSignal2));
    
    % Step 3: Calculate phase difference
    phaseDiff = phase1 - phase2;
    
    % Step 4: Convert to complex domain and compute average phase
    complexSeq = exp(1i * phaseDiff);
    avgComplex = mean(complexSeq);
    avgPhase = abs(angle(avgComplex));

end

%% Example plotting code (commented out)
% figure(10);
% subplot(4,6,chooseIdx); imagesc(correlationMatrix(:,:,chooseIdx));
% caxis([0.4 1.1]); colormap jet; title(sprintf('Correlation - %d Pixels', chooseNum(chooseIdx)));
% 
% figure(11);
% subplot(4,6,chooseIdx); imagesc(phaseDiffMatrix(:,:,chooseIdx));
% caxis([0 0.2]); colormap jet; title(sprintf('Phase Difference - %d Pixels', chooseNum(chooseIdx)));