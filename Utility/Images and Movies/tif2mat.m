function [ Movie, Info ] = tif2mat( fname, varargin )
%TIF2MAT - Read a TIFF file containing multiple frames and convert to a mat file
%   
% Syntax:
% -------
% [ Movie, info ] = tif2mat( fname, MaxFrames )
%
% Input:
% ------
% fname     - Name of the TIFF file. Must contain the .tif suffix
% MaxFrames - Maximum number of frames to read from the first frame
%
% Outputs:
% --------
% Movie - Mat file containing the movie
% Info  - Movie's information from the TIFF file
%
% Ver 1. Written by Oren Solomon, Technion I.I.T.
% Ver 2. Added ability to read a pre-specified number of frames
%

% Initialization
if strcmp(fname(end-2:end), 'tif')
    UseName = fname;
else
    UseName = [fname '.tif'];
end
Info = imfinfo(UseName);

if nargin > 1
    num_images = varargin{1};
else
    num_images = numel(Info);
end

% Allocate memory
Movie = zeros(Info(1).Height, Info(1).Width, num_images);

% Conver each frame
for kk = 1:num_images
    Movie(:, :, kk) = imread(UseName, kk);
end

% output
% Movie = imread(UseName);