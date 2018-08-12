function [ tStr ] = SearchLine(FileName, match_str, SwapValue)
% SEARCHLINE Searches for a line in "fid" which begins with "match_str" and reads the
% relevant parameter
%
% Syntax:
% -------
% [ tStr ] = SearchLine(FileName, match_str, SwapValue)
%
% Inputs:
% -------
% FileName  - Input file name to read / write
% match_str - Parameter name to seek in the file
% SwapValue - (optional) A new parameter value for match_str
%
% Output:
% -------
%  tStr     - Output value of the desired parameter ([] if the user only wants to swap a line)
%
% Ver 1. Written by Oren Solomon, Technion I.I.T. 16-10-2016
%

%% ---------------------- Initializations ---------------------------------
% Check number of inputs - determines type of function operation
if nargin > 2
    SwapValueFlag = 1;
else 
    SwapValueFlag = 0;
end

% Open file to read
fid = fopen(FileName, 'r');

% Initialize tStr
tStr = [];

% Text file buffer 
FileBuffer = {};
kk = 1;

%% ------------------- Read each line in file -----------------------------
% Set cursor to beginning of file
fseek(fid, 0, 'bof');

while ~feof(fid) 
    % Read a line from the file
    tline = fgetl(fid);
    
    % File buffer
    FileBuffer{kk} = tline;
    kk = kk + 1;
    
    % Determine type of line
    if ~isempty(tline)                                                          % Not an empty line
        if tline(1) ~= '%'                                                      % Line is not a comment line
            % Line beginning with '%' is a comment and should be disregarded
            str_ind  = strfind(tline, match_str);
            
            % We found the string that we wanted
            if ~isempty(str_ind)
                ind1 = strfind(tline, '=');
                ind2 = strfind(tline, ';');                                    % Every line must end with ';'
                
                if ~SwapValueFlag
                    % Read the value
%                     tStr = tline(ind1(1) + 1:ind2 - 1);
                    tStr = strtrim(tline(ind1(1) + 1:ind2 - 1));
                    
                    % Clean ' if there are any in the string
                    k    = strfind(tStr, '''');
                    if ~isempty(k)
                        tStr = tStr(k(1) + 1:k(2) - 1);
                    end

                    % Exit the loop
                    break;
                else
                    % Rewrite the value
                    tline = [tline(1:ind1(1)) ' ' num2str(SwapValue) ';'];
                    
                    % Replace line
                    FileBuffer{kk - 1} = tline;
                end
            end
        end
    end
end

% Close file
fclose(fid);

%% ------------------ Rewrite file when needed ----------------------------
% Rewrite the file, if needed
if SwapValueFlag
    % Open file for writing
    fid = fopen(FileName, 'w');
    
    % Write new line, one at a time
    for kk = 1:numel(FileBuffer) - 1
        fprintf(fid, '%s\n', FileBuffer{kk});
    end
    
    % Write last line
    fprintf(fid, '%s', FileBuffer{end}); 
end

% Close file
fclose(fid);