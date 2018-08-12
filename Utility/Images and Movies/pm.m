function pm( MovieMat, varargin )
%PLAYMOVIE - Play movie stored in MovieMat with specific rate and colormap
%
% Syntax:
% -------
% pm( MovieMat, MovSpeed, ColMap )
%
% Inputs:
% -------
% MovieMat - M X N X L stack. Each l = 1,...,L is an M X N frame.
% MovSpeed - Movie speed to play [sec]. Default value is 1/24 sec.
% ColMap   - Colormap for the movie. Default is 'hot'.
%
% Written by Oren Solomon, Technion I.I.T., Ver 1.
% 

% Handle inputs
if nargin == 1                                              % Only MovieMat as input, set other values to default values
    MovSpeed = 1/24;
    ColMap = 'hot';
elseif nargin == 2                                          % Set only colormap to default
    MovSpeed = varargin{1};
    ColMap = 'hot';
elseif nargin == 3                                          % User specified frame-rate and colormap
    MovSpeed = varargin{1};
    ColMap = varargin{2};
end

% Movie length
[~, ~, MovLength] = size(MovieMat(1, 1, :));

% Play the movie
figure;
eval(['colormap ' ColMap ';']);                             % Set colormap
for ii = 1:MovLength
    imagesc(MovieMat(:, :, ii));                            % Show frame
    axis square;
    title(['Frame ' num2str(ii) '/' num2str(MovLength)]);   % Frame number
    pause(MovSpeed);                                        % Movie speed
end