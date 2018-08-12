function mat2tif( MatIn, FileName, varargin )
%MAT2TIFF - converts a movie stored as a 3D matrix to a 16bit tiff format.
%
% Syntax:
% -------
% mat2tif( MatIn, FileName, LocFolder )
%
% Inputs:
% -------
% MatIn     - Input movie
% FileName  - Desired output file name
% LocFolder - Location in which to create the tiff (optional)
%
% Ver 1. WRitten by Oren Solomon, Technion I.I.T. 18-01-2016 
%

if nargin == 3
    LocFolder = varargin{1};
else
    LocFolder = '.';
end

% Determine size
[m, n, k] = size(MatIn);

% Convert to unit16 for 16 bit resolution
% MatIn = uint16(MatIn);

% header
header = [];

% Write first image
% imwrite(MatIn(:, :, 1), [LocFolder '\' FileName '.tif'], 'Compression', 'none', 'Resolution', [m, n]);
imwrite2tif(MatIn(:, :, 1), header, [LocFolder '\' FileName '.tif'], 'uint64');


% Append rest of images
for ii = 2:k
%     imwrite(MatIn(:, :, ii), [LocFolder '\' FileName '.tif'], 'WriteMode', 'append', 'Compression', 'none', 'Resolution', [m, n]);
    imwrite2tif(MatIn(:, :, ii), header, [LocFolder '\' FileName '.tif'], 'uint64');
%     pause(0.5);
end
