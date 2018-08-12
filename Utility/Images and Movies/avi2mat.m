function [ MovieMat ] = avi2mat( AVIname, varargin )
% AVIREAD - Read AVI file and convert to a mat file. Play the video if needed
%
% Syntax:
% -------
% [ MovieMat ] = avi2mat( AVIname, FileFolder, MaxFrames, PlaybackFlag )
%
% Inputs:
% -------
% AVIname      - Name of the input AVI file
% FileFolder   - Location of the AVI file (optional). Default is current directory
% MaxFrames    - Maximum nuber of frames to read (optional). Defualt is the entire length of the movie
% PlaybackFlag - Play the movie flag (optional). Default is 0 - no
%
% Outputs:
% --------
% MovieMat     - A mat file containing the movie. Each slice is a single grayscale frame
%
% Ver 1. WRitten by Oren Solomon with the MATLAB help. 08-02-2016
%

% Input checks
if nargin == 2
    FileFolder    = varargin{1};
    MaxFrames     = [];
    PlaybackFlag  = 0;
elseif nargin == 3
    FileFolder    = varargin{1};
    MaxFrames     = varargin{2};
    PlaybackFlag  = 0;
elseif nargin == 4
    FileFolder    = varargin{1};
    MaxFrames     = varargin{2};
    PlaybackFlag  = varargin{3};
else
    FileFolder    = '.';
    MaxFrames     = [];
    PlaybackFlag  = 0;
end

% Read AVI file
xyloObj = VideoReader( fullfile(FileFolder,AVIname) );

% Determine size of the movie
nFrames = xyloObj.NumberOfFrames;
vidHeight = xyloObj.Height;
vidWidth = xyloObj.Width;

% Determine whether number of frames to read was specified by the user or if it
% exceeds the number of frames in the file
if isempty(MaxFrames) || MaxFrames > nFrames
    MaxFrames = nFrames;
end

% Preallocate the movie structure
if PlaybackFlag
    mov(1:MaxFrames) = struct('cdata', zeros(vidHeight,vidWidth, 3, 'uint8'), 'colormap',[]);
end
MovieMat = zeros(vidHeight, vidWidth, MaxFrames);

% Read one frame at a time.
for kk = 1 : MaxFrames
    tmp = read(xyloObj,kk);
    if PlaybackFlag
        mov(kk).cdata = tmp;
    end
    
    % Read the movie to a mat file
    MovieMat(:, :, kk) = tmp(:, :, 1);
end

if PlaybackFlag
    % Size a figure based on the video's width and height.
    hf = figure;
    set(hf, 'position', [150 150 vidWidth vidHeight])
    
    % Play back the movie once at the video's frame rate.
    movie(hf, mov, 1, xyloObj.FrameRate);
end