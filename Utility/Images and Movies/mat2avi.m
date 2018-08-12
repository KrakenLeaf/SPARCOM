function mat2avi( InputMovie, FileName, varargin )
% MAT2AVI creates an AVI file from a fiven 3D matrix
%
% Syntax:
% -------
% mat2avi( InputMovie, FileName, Colormap, NumberOfFrames )
%
% Inputs:
% -------
% InputMovie     - 3D stack to record as an AVI movie
% FileName       - Name of AVI movie 
% Colormap       - (optional) specific colormap to use. Default is 'hot'
% NumberOfFrames - (optional) number of frames to write
%
% V.1 Written by Oren Solomon, Technion I.I.T. 21-11-2016
%

% Input parsing
try 
    Colormap = varargin{1};
catch
    Colormap = 'hot';
end
try 
    NumberOfFrames = varargin{2};
catch
    NumberOfFrames = size(InputMovie, 3);
end

% Convert movie to double
InputMovie = double(InputMovie); 

% Normalization factor - values must be in the range [0, 255]
MaxVal = max(max(max(InputMovie)));
InputMovie = InputMovie/MaxVal*255;

% Open AVI object
v = VideoWriter([FileName '.avi'], 'Indexed AVI');
eval(['v.Colormap = ' Colormap '(256);']);

% Open movie object
open(v);

% Write each frame to the movie
for ii = 1:NumberOfFrames
    writeVideo(v, uint8(InputMovie(:, :, ii)));
end

% Close movie object
close(v);




