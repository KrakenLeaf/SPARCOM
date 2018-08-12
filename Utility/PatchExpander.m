function [ PatchOut ] = PatchExpander( PatchIn, Pix, method )
%PATCHEXPANDER - expands a patch in its borders according to the specified method.
%
% Syntax:
% -------
% [ PatchOut ] = PatchExpander( PatchIn, Pix, method )
%
% Inputs:
% -------
% PatchIn  - Input patch to manipulate 
% Pix      - Border length for expansion
% Method   - Expansion method:
%            'same'    - Do nothing
%            'neumann' - Copy patch borders (gradient is zero)
%            'shrink'  - Cut a center from the patch
%            'zero'    - Pad with zeros
%
% Output:
% -------
% PatchOut - Patch after expansion
%
% Written by Oren Solomon, Technion I.I.T, Ver 1.
%

%% Initializations
[M, N, L] = size(PatchIn);

switch lower(method)
    % Same - do nothing
    case 'same'
        TmpImage = PatchIn;
    % Neumann - duplicate the edges by a factor perc
    case 'neumann'
        TmpImage = [PatchIn(:, 1:Pix, :) PatchIn PatchIn(:, end - Pix + 1:end, :)];
        TmpImage = [TmpImage(1:Pix, :, :); TmpImage; TmpImage(end - Pix + 1:end, :, :) ];
    % Shrink - restore patch to its original size
    case 'shrink' 
        TmpImage = PatchIn(Pix + 1: end - Pix, Pix + 1: end - Pix, :);
    % Zero - put zeros in the boundaries
    case 'zero'
        TmpImage = zeros(M, N, L);
        TmpImage(Pix + 1:end - Pix, Pix + 1:end - Pix, :) = PatchIn(Pix + 1:end - Pix, Pix + 1:end - Pix, :);
    otherwise
        error('PatchExpander: Method not supported.');
end

% Output
PatchOut = TmpImage;