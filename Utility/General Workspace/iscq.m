function iscq( ImageIn, varargin )
% ISCQ - Quick display of grayscale figures (and other options)
%
% Syntax:
% -------
% iscq( ImageIn, title, colormap )
%
% Inputs:
% -------
% ImageIn  - Arbitrary image to display
% title    - (Optional) Image title
% colormap - (Optional) Change colormap (default is 'gray')
%
% Ver 1: Written by Oren Solomon, Technion, I.I.T, 12-10-2015
%

%% Handle inputs
% --------------------------------------------------------
if nargin > 1
    % Title
    TitleUse = varargin{1};
    
    % Colormap
    try
        ColormapUse = varargin{2};
    catch
        ColormapUse = (((0:2^10-1)/(2^10-1))'*[1 1 1]).^.4;
    end
else
    % Default values
    TitleUse = [];
    ColormapUse = (((0:2^10-1)/(2^10-1))'*[1 1 1]).^.4;
end

%% Display figure
% --------------------------------------------------------
figure;
try
    colormap(ColormapUse);
catch
    eval(['colormap' ColormapUse ';']);
end
imagesc(ImageIn);axis image;
title(TitleUse);
