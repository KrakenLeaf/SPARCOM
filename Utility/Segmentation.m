function [ E_out ] = Segmentation( E_in, P )
%SEGMENTATION - performs binary segmentation of the image, keeps only the
%largest P*ValMsx entries.
%
% Format:
% -------
% [ E_out ] = Segmentation( E_in, P )
%
% INPUTS:
% -------
% E_in          - Input vector (or image).
% P             - A calue between [0,1]. If P = 0, E_out = E_in. If P = 1, E_out
%                 is a sparse vector with non-zero entries which are the largest value of E_in
%
% OUTPUTS:
% --------
% E_out         - Output result, same dimentions as E_in.
%

if (P < 0 | P > 1)
   error('SEGMENTATION: Percentage must be in [0,1].'); 
end

% Take size of the input
[M, N] = size(E_in);

% Transform to vector
E_in_temp = E_in(:);

% Find largest value
MaxVal = max(E_in_temp);

% Place zero in every entry which is lower than P*MaxVal
Ind = find(E_in_temp < P*MaxVal);

E_in_temp(Ind) = 0;

% Output
E_out = reshape(E_in_temp,M,N);

