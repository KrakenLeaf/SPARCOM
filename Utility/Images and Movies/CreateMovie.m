function CreateMovie(MatIn, MovieName, varargin)
% CREATEMOVIE - Generate an AVI movie from a given cube (each 2D slice is assume to be a single frame)
% 
% Syntax:
% -------
% CreateMovie(MatIn, MovieName, SaveFolder, ColormapUse, DupFramesNum)
%
% Inputs:
% -------
% MatIn        - Input cube [M, N, K]
% MovieName    - Name of the output movie
% SaveFolder   - (optional) Location of the output move. Default is current directory
% ColormapUse  - (optional) Which colormap to use. Default is 'hot'
% DupFramesNum - (optional) How many times each frame will be duplicated, for slower frame-rate. Default is 2
%
% Ver 1. Written by Oren Solomon, Technion I.I.T. 24-10-2015
%

%% Initialization
[~, ~, K] = size(MatIn);

% Extract save folder destination
if nargin > 2
    SaveFolder = varargin{1};
else
    SaveFolder = '.';
end

% Extract colormap 
if nargin > 3
    ColormapUse = varargin{2};
else
    ColormapUse = 'hot';
end

% Extract number of duplicate frames
if nargin > 4
    DupFramesNum = varargin{3};
else
    DupFramesNum = 2;
end

%% Capture each frame
disp('Playing movie...');
h1 = figure;
for ii = 1:K
    eval(['colormap ' ColormapUse ';']);
    imagesc(MatIn(:, :, ii));
    set(gca,'XTickLabel','');set(gca,'YTickLabel','');
    axis image; drawnow;
    CapMov(ii) = getframe;
end

% Expand indices by factor of 2 - play movie slower
ind = ones(DupFramesNum,1)*(1:K);
ind = [ind(:);ind(:)];
CapMov = CapMov(ind(:));

%% Generate the movie
disp('Converting movie to AVI...');
movie2avi(CapMov, fullfile(SaveFolder, MovieName), 'compression', 'None');
close(h1);

disp('Done.');