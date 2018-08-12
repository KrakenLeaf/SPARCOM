function mat2txt( MatIn, TXTname, varargin )
% MAT2TXT - Writes a movie in MatIn (each slice is a frame) to a text file.
%
% Syntax;
% -------
% mat2txt( MatIn, TXTname, TXTfolder, PrecType )
%
% Inputs:
% -------
% MatIn     - Input 3D matrix of arbitrary dimensions. Each slice is a frame, starting from 1.
% TXTname   - Desired name for the output text file.
% TXTfolder - Desired location for the text file (optional). Default is current directory.
% PrecType  - Desired writing precision (optional). Default is 'uint8'.
%
% Ver 1. Written by Oren Solomon, Technion I.I.T. 03-02-2016
%

% Set save directory for the text file and desired precision
if nargin == 3
    TXTfolder = varargin{1};                                    % Save destination folder
    PrecType = 'uint8';                                         % Writing precision
elseif nargin == 4
    TXTfolder = varargin{1};
    PrecType = varargin{2};
else
    TXTfolder = '.';
    PrecType = 'uint8';
end

% Dimensions of MatIn
k = size(MatIn, 3);

% Open new file for writing
% fid = fopen(fullfile(TXTfolder, [TXTname '.txt']), 'w');
fid = fopen(fullfile(TXTfolder, TXTname), 'w');

% Read each frame in the 3D matrix and convert it to a line in the text file
for ii = 1:k
    switch PrecType
        case 'uint8'
            fwrite(fid, uint8(MatIn(:, :, ii)), PrecType);
        case 'uint16'
            fwrite(fid, uint16(MatIn(:, :, ii)), PrecType);
        case 'uint32'
            fwrite(fid, uint32(MatIn(:, :, ii)), PrecType);
        case 'uint64'
            fwrite(fid, uint64(MatIn(:, :, ii)), PrecType);
        otherwise
            error('MAT2TXT: Unknown precision type.');
    end
end

% Close the text file
fclose(fid);