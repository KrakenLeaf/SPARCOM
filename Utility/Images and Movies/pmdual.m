function pmdual( MovieMat1, MovieMat2, varargin )
%PLAYMOVIE - Play two movies stored in MovieMat1 and MovieMat2 simultaneously with specific rate and colormap
%
% Syntax:
% -------
% pm( MovieMat1, MovieMat2, varargin )
%
% Inputs:
% -------
% MovieMat1 - M X N X L stack. Each l = 1,...,L is an M X N frame.
% MovieMat2 - M X N X L stack. Each l = 1,...,L is an M X N frame.
% MovSpeed  - Movie speed to play [sec]. Default value is 1/24 sec.
% ColMap    - Colormap for the movie. Default is 'hot'.
%
% Written by Oren Solomon, Technion I.I.T., Ver 1.
% 

% Handle inputs
if nargin == 2                                              % Only MovieMat as input, set other values to default values
    MovSpeed = 1/24;
    ColMap = 'hot';
elseif nargin == 3                                          % Set only colormap to default
    MovSpeed = varargin{1};
    ColMap = 'hot';
elseif nargin == 4                                          % User specified frame-rate and colormap
    MovSpeed = varargin{1};
    ColMap = varargin{2};
end

% Movie length
[~, ~, MovLength1] = size(MovieMat1(1, 1, :));
[~, ~, MovLength2] = size(MovieMat2(1, 1, :));
MovLength = min(MovLength1, MovLength2);

% Play both movies
figure;
eval(['colormap ' ColMap ';']);                             % Set colormap
for ii = 1:MovLength
    subplot(1,2,1);
    imagesc(MovieMat1(:, :, ii));                            % Show frame
    axis image;
    title(['Frame ' num2str(ii) '/' num2str(MovLength)]);   % Frame number
    subplot(1,2,2);
    imagesc(MovieMat2(:, :, ii));                            % Show frame
    axis image;
    title(['Frame ' num2str(ii) '/' num2str(MovLength)]);   % Frame number
    pause(MovSpeed);                                        % Movie speed
end