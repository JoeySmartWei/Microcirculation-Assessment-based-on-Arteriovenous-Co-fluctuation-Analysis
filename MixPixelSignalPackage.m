function [mixedCorrelationMatrix, mixedPhaseDiffMatrix] = MixPixelSignalPackage(roiSignal)
% MIXPIXEL SIGNALPACKAGE - Generate mixed arterial/venous pixel signals and analyze their association
%
% Syntax:
%   [mixedCorrelationMatrix, mixedPhaseDiffMatrix] = MixPixelSignalPackage(roiSignal)
%
% Inputs:
%   roiSignal - 3D matrix (timeLength × pixelNumber × typeNumber), ROI signals (type 1: artery, type 2: vein)
%
% Outputs:
%   mixedCorrelationMatrix - 3D matrix (repeats × repeats × mixNum), correlation of mixed signals
%   mixedPhaseDiffMatrix   - 3D matrix (repeats × repeats × mixNum), phase difference of mixed signals
%
% Description:
%   Randomly mix arterial and venous pixel signals, calculate correlation and phase difference
%   across different mixing ratios, and output the association matrices.

    % Get basic dimensions
    [timeLength, pixelNumber, numTypes] = size(roiSignal);
    
    % Configuration parameters
    numRepeats = 20;                  % Number of random mixing repetitions
    chooseNum = 1:(pixelNumber - 1);  % Range of pixel numbers for mixing
    mixNum = 1:length(chooseNum);     % Mixing ratio range
    
    % Reference signal (mean of arterial pixels, type 1)
    refSignal = mean(roiSignal(:, :, 1), 2);
    
    % Generate mixed ROI signals with random masks
    organizedMixSignal = generateMixROIMask(roiSignal, mixNum, numRepeats);
    
    % Initialize output matrices
    mixedCorrelationMatrix = nan(numRepeats, numRepeats, length(mixNum));
    mixedPhaseDiffMatrix = nan(numRepeats, numRepeats, length(mixNum));
    
    % Calculate association for each mixing ratio
    for mixIdx = 1:length(mixNum)
        currMixSignal = organizedMixSignal(:, :, mixIdx);
        [corrMat, phaseDiffMat] = calculateAssociation(currMixSignal);
        mixedCorrelationMatrix(:, :, mixIdx) = corrMat;
        mixedPhaseDiffMatrix(:, :, mixIdx) = phaseDiffMat;
    end

end

% -------------------------------------------------------------------------
function [correlation, phaseDiff, variance, innerVariance, innerPhaseDiff] = randomPixel(mixedSignal, chooseNum, numRepeats, refSignal)
% RANDOMPIXEL - Randomly sample mixed signals and calculate statistical metrics
%
% Inputs:
%   mixedSignal - 2D matrix (timeLength × indexNumber), mixed signal matrix
%   chooseNum   - Scalar, number of pixels to select
%   numRepeats  - Scalar, number of sampling repetitions
%   refSignal   - 1D vector, reference signal for comparison
%
% Outputs:
%   correlation   - Mean correlation with reference signal
%   phaseDiff     - Mean phase difference with reference signal
%   variance      - Mean variance of (signal - reference)
%   innerVariance - Mean internal variance of sampled signals
%   innerPhaseDiff- Mean internal phase difference of sampled signals

    [timeLength, indexNumber] = size(mixedSignal);
    
    % Set sampling gap (optimize computation)
    gap = indexNumber > 50 ? indexNumber / 100 : 1;
    
    % Only process the last chooseNum (consistent with original logic)
    currChooseNum = chooseNum(end);
    selectSignal = zeros(timeLength, numRepeats);
    
    % Random sampling
    for rdmIdx = 1:numRepeats
        selectIndex = randperm(indexNumber, currChooseNum)';
        selectSignal(:, rdmIdx) = mean(mixedSignal(:, selectIndex), 2);
    end
    
    % Calculate statistical metrics
    correlation = mean(corr(selectSignal, refSignal));
    phaseDiff = mean(arrayfun(@(x) calculatePhaseDifference(selectSignal(:, x), refSignal), 1:numRepeats));
    variance = mean(var(selectSignal - refSignal, [], 1));
    innerVariance = mean(var(selectSignal, [], 1));
    innerPhaseDiff = mean(arrayfun(@(x) calculatePhaseDifference(selectSignal(:, x)), 1:numRepeats));

end

% -------------------------------------------------------------------------
function organizedMixSignal = generateMixROIMask(roiSignal, mixNum, numRepeats)
% GENERATEMIXROIMASK - Generate random mixed ROI signals (artery + vein)
%
% Inputs:
%   roiSignal   - 3D matrix (timeLength × pixelNumber × typeNumber), original ROI signals
%   mixNum      - Array, range of mixing ratios
%   numRepeats  - Scalar, number of repetitions for each mixing ratio
%
% Outputs:
%   organizedMixSignal - 3D matrix (timeLength × numRepeats × length(mixNum)), mixed signals

    % Separate arterial/venous signals (type 1: fixed/artery, type 2: moved/vein)
    fixedROISignal = roiSignal(:, :, 1);
    movedROISignal = roiSignal(:, :, 2);
    [timeLength, pixelNumber] = size(fixedROISignal);
    
    % Set sampling gap
    gap = pixelNumber > 50 ? pixelNumber / 100 : 1;
    
    % Initialize mixed signal matrix
    organizedMixSignal = zeros(timeLength, numRepeats, length(mixNum));
    
    % Generate mixed signals for each ratio
    for mixIdx = 1:gap:length(mixNum)
        currMixNum = mixNum(mixIdx);
        
        for rdmIdx = 1:numRepeats
            % Randomly select venous/arterial pixels
            movedIdx = randperm(pixelNumber, currMixNum)';
            fixedIdx = randperm(pixelNumber, pixelNumber - currMixNum)';
            
            % Mix signals and randomly select 12 pixels for mean calculation
            mixedSignal = [fixedROISignal(:, fixedIdx), movedROISignal(:, movedIdx)];
            random12Idx = randperm(size(mixedSignal, 2), 12);
            organizedMixSignal(:, rdmIdx, mixIdx) = mean(mixedSignal(:, random12Idx), 2);
        end
    end

end

% -------------------------------------------------------------------------
function [corrMat, phaseDiffMat] = calculateAssociation(signalMatrix)
% CALCULATEASSOCIATION - Calculate pairwise correlation and phase difference (reusable function)
%
% Inputs:
%   signalMatrix - 2D matrix (timeLength × signalNumber), signals to compare
%
% Outputs:
%   corrMat      - 2D matrix (signalNumber × signalNumber), correlation coefficients
%   phaseDiffMat - 2D matrix (signalNumber × signalNumber), phase differences

    numSignals = size(signalMatrix, 2);
    corrMat = nan(numSignals, numSignals);
    phaseDiffMat = nan(numSignals, numSignals);
    
    % Pairwise calculation (upper triangle for phase, lower for correlation)
    for i = 1:numSignals
        for j = (i + 1):numSignals
            corrMat(j, i) = corr(signalMatrix(:, i), signalMatrix(:, j));
            phaseDiffMat(i, j) = calculatePhaseDifference(signalMatrix(:, i), signalMatrix(:, j));
        end
    end

end

% -------------------------------------------------------------------------
function avgPhase = calculatePhaseDifference(signal1, signal2)
% CALCULATEPHASEDIFfERENCE - Compute average phase difference (supports 1/2 input signals)
%
% Inputs:
%   signal1 - 1D vector, primary signal
%   signal2 - 1D vector (optional), secondary signal (if omitted, use signal1's phase)
%
% Outputs:
%   avgPhase - Scalar, average phase difference (absolute value)

    if nargin == 2
        % Two-input mode: phase difference between two signals
        hSignal1 = hilbert(signal1);
        hSignal2 = hilbert(signal2);
        phase1 = unwrap(angle(hSignal1));
        phase2 = unwrap(angle(hSignal2));
        phaseDiff = phase1 - phase2;
    else
        % Single-input mode: phase of single signal
        phaseDiff = unwrap(angle(hilbert(signal1)));
    end
    
    % Convert to complex domain and compute average
    complexSeq = exp(1i * phaseDiff);
    avgComplex = mean(complexSeq);
    avgPhase = abs(angle(avgComplex));

end