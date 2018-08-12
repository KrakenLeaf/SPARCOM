function [ Movie ] = txt2mat( txtFile, xdim, ydim, varargin )
% TXT2MAT - Read movie stored in txt file and convert it to a mat file
%
% Syntax:
% -------
% [ Movie ] = txt2mat( txtFile, xdim, ydim, Nfr )
% [ Movie ] = txt2mat( txtFile, xdim, ydim, Nst, Nfr ) - read frames Nst:Nfr
%
% Inputs:
% -------
% txtFile - Name of input txt file
% xdim    - x dimension
% ydim    - y dimension
% Nst     - Initial frame number
% Nfr     - Number of frames to read (end frame number)
%
% Outputs:
% --------
% Movie   - [xdim, ydim, Nfr] Movie
%
% Ver 1. Written by Oren Solomon,  Technion I.I.T 08-11-2015
% Ver 2. Written by Oren Solomon,  Technion I.I.T 17-01-2016: Added ability to
% start from a specific frame
%

if nargin == 5
    Nst = varargin{1};
    Nfr = varargin{2};
else
    Nst = 1;
    Nfr = varargin{1};
end

% Initialize Movie
Movie = zeros(xdim, ydim, Nfr - Nst + 1);

% Open file in READ mode
fid = fopen(txtFile, 'r');

kk = 1;
for ii = Nst:Nfr
    % Set file indicator at beginning of line (relative to beginning of file)
    fseek(fid, xdim*ydim*(ii - 1), 'bof');
    
    % Read line, convert it to a double precision [xdim, ydim] matrix
    Movie(:, :, kk) = double(reshape(fread(fid,xdim*ydim,'uint8=>uint8'),xdim,ydim));
    
    % Advance frame number
    kk = kk + 1;
end

% Close file
fclose(fid);







% % % % Initialize Movie
% % % Movie = zeros(xdim, ydim, Nfr);
% % % 
% % % % Open file in READ mode
% % % fid = fopen(txtFile, 'r');
% % % 
% % % for ii = 1:Nfr
% % %     % Set file indicator at beginning of line (relative to beginning of file)
% % %     fseek(fid, xdim*ydim*(ii - 1), 'bof');
% % %     
% % %     % Read line, convert it to a double precision [xdim, ydim] matrix
% % %     Movie(:, :, ii) = double(reshape(fread(fid,xdim*ydim,'uint8=>uint8'),xdim,ydim));
% % % end
% % % 
% % % % Close file
% % % fclose(fid);