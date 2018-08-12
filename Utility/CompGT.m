function [x_start, x_end, y_start, y_end, GT_patch, Plots] = CompGT( DataIn, GT, SRF, Pos, CorrFact, PlotFlag )
%COMPGT - Compare diffraction limited data / reconstructed HR image with ground truth.
%
% Syntax:
% -------
% [x_start, x_end, y_start, y_end] = CompGT( DataIn, GT, SRF, Pos, CorrFact, PlotFlag )
%
% Inputs:
% -------
% DataIn       - Input movie (low resolution). Size of [M, N, K]
% GT           - Ground truth image
% SRF          - Super resolution factor (by how much to enlarge DataIn such that the dimensions of DataIn and GT coincide)
% Pos          - Pos(1) starting row position for desired patch from DataIn.
%                Pos(2) starting column position for desired patch from DataIn.
%                Pos(3) size of slice (square).
% CorrFact     - By how much to enlarge Pos(3). For the purpose of correlation calculation only
% PlotFlag     - 1 - show plots, 0 - don't show plots
%
% Outputs:
% --------
% x_start      - Beginning of aligned patch (row)
% x_end        - end of aligned patch (row)
% y_start      - Beginning of aligned patch (column)
% y_end        - end of aligned patch (column)
% GT_patch     - Ground truth aligned patch
% Plots        - Cell array of different plots of the function (only if PlotFlag is TRUE) - not including the correlation image
%
% Ver 1: Written by Oren solomon (with some help from MATLAB's help & examples), Technion I.I.T
%

% Interpolation safety margin
Marg = 0;

%% DataIn is a diffraction limited (low resolution) movie
% -----------------------------------------------------------------------------------------------------------------------------------------------------------
% Sum over the 3rd dimension - relevant only if this is a diffraction limited movie (Raw data)
ImSum = sum(DataIn, 3);

% Data interpolated to the size of GT
EnlargedData = imresize(ImSum, SRF);

% Take chunk location and size in LR movie
Xp      = Pos(1);
Yp      = Pos(2);
BlkSize = Pos(3)*CorrFact + Marg;                            % Add additional pixels in order to avoid interpolation artifacts at the edges

% Cut the chunk and resize it according to SRF (and subtract mean)
Chunk = imresize(ImSum(Xp:Xp + BlkSize - 1, Yp:Yp + BlkSize - 1), SRF);
Chunk = Chunk(1: end - Marg*SRF, 1:end - Marg*SRF);       % Due to interpolation artifacts at the borders

% Chunk used for correlation (big)
Plots{1} = Chunk;
if PlotFlag; figure; colormap gray; imagesc(Plots{1}); title(['Selected Block after interpolation by factor ' num2str(SRF) ' - large (X' num2str(CorrFact) ')']); end;

% Desired patch (small)
Patch = imresize(ImSum(Xp:Xp + Pos(3) + Marg - 1, Yp:Yp + Pos(3) + Marg - 1), SRF);
Patch = Patch(1: end - Marg*SRF, 1:end - Marg*SRF);
Plots{2} = Patch;
if PlotFlag; figure; colormap gray; imagesc(Plots{2}); title(['Selected Block after interpolation by factor ' num2str(SRF) ' - small']); end;

% Find mean of the entire enlarged image
MeanTot = mean(mean(EnlargedData));

% Subtract mean
EnlargedData = EnlargedData - MeanTot;
Chunk        = Chunk - MeanTot;

% Perform correlation with enlarged image to find the relevant locations in GT
crr = xcorr2(EnlargedData, Chunk);

% Find maximum
[ssr, snd] = max(crr(:));                            % ssr - maximum value, snd - location of ssr in the vector
[ij, ji]   = ind2sub(size(crr), snd);

% Display correlation graph
if PlotFlag
    figure; plot(crr(:)); grid on;
    title('Cross-Correlation - High resolution (interpolation)');
    hold on; plot(snd,ssr,'or');
    hold off; text(snd*1.05,ssr,'Maximum');
end

% Beginning & end of the chunk - in the enlarged dimensions
x = ij - (BlkSize - Marg)*SRF + 1;
X = ij;
y = ji - (BlkSize - Marg)*SRF + 1;
Y = ji;
X2 = x + ((BlkSize - Marg)/CorrFact)*SRF - 1;
Y2 = y + ((BlkSize - Marg)/CorrFact)*SRF - 1;

% Display results - sanity check
Plots{3} = imfuse(EnlargedData, GT);
if PlotFlag
    figure; imagesc(Plots{3});
    colormap gray; title('Super position of interpolated LR and ground truth');
    hold on; plot([y y Y Y y],[x X X x x],'r'); hold off;
    hold on; plot([y y Y2 Y2 y],[x X2 X2 x x],'y'); hold off;
end

%% Arrange output
% -----------------------------------------------------------------------------------------------------------------------------------------------------------
x_start = x;
x_end = X2;
y_start = y;
y_end = Y2;

GT_patch = GT(x_start:x_end, y_start:y_end);

% Plot diffractino limited patch with ground truth
Plots{4} = imfuse(Patch, GT_patch);
if PlotFlag; PlotsHandles{4} = figure; colormap gray; imagesc(Plots{4}); title(['Selected Block after interpolation by factor ' num2str(SRF) ' with GT']); end;




