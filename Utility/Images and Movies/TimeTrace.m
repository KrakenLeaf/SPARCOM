function TimeTrace( MovieIn, varargin )
%TIMETRACE Summary of this function goes here
%   Detailed explanation goes here
%
% Syntax:
% -------
% TimeTrace( MovieIn, PixelsToShow, CalcBias )
%
% Inputs:
% -------
% MovieIn      - Input movie to display the time traces
% PixelsToShow - A matrix if size [K, 2] - which pixels' time trace to display. Leave empty in
%                order to go over all pixels in the movie (optional)
% CalcBias     - 1 means calculate bias and subtract from the pixels, 0 - no. Relevant only when PixelsToShow is not empty (optional)
%
% Ver 1: Written by Oren Solomon, Technion I.I.T, 05-10-2015
%


% Determine size of input movie
[M, N, K] = size(MovieIn);

figure;hold on;grid on;
if nargin == 1
% Show all pixles' time traces
% ----------------------------------------------------------------------------------------------
    ColUse = jet(M*M);
    kk = 1;
    for ii = 1:M
        for jj = 1:N
            plot(1:K, squeeze(MovieIn(ii, jj, :)), 'color', ColUse(kk, :));
            kk = kk + 1;
        end
    end
else
% Show selected pixles' time traces
% ----------------------------------------------------------------------------------------------
    [Mv, Nv] = size(varargin{1});
    
    % Calculate bias flag
    try
        BiasFlag = varargin{2};
    catch
        BiasFlag = 0;
    end
        
    ColUse = jet(Mv);
    LegendUse = {};
    for ii = 1:Mv
        if BiasFlag == 1
        % Calculate bias and subtract it for each pixel
        % --------------------------------------------------------------------------------------
            TmpBias = mean(MovieIn(varargin{1}(ii, 1), varargin{1}(ii, 2), :));
            plot(1:K, squeeze(MovieIn(varargin{1}(ii, 1), varargin{1}(ii, 2), :)) - TmpBias, 'color', ColUse(ii, :));
            TitleAdd = [', Bias subtracted'];
        else
        % Do not calculate bias
        % --------------------------------------------------------------------------------------
            plot(1:K, squeeze(MovieIn(varargin{1}(ii, 1), varargin{1}(ii, 2), :)), 'color', ColUse(ii, :));
            TitleAdd = '';
        end
        LegendUse{ii} = ['Pixel [' num2str(varargin{1}(ii, 1)) ' ,' num2str(varargin{1}(ii, 2)) ']'];
    end
    legend(LegendUse);
end
xlabel('Time [AU]');ylabel('Intensity [AU]');title(['Time traces of chosen pixels in the Movie' TitleAdd]);