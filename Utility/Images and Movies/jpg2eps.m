function jpg2eps( jpeg_file )
%JPG2EPS Transforms a jpeg image to an eps file, for LaTeX.
%
% Syntax:
% -------
% [ flag ] = jpg2eps( jpeg_file )
%
% Input:
% ------
% jpeg_file - jpeg file to convert to eps
%
% Ver 1: WRitten by Oren SOlomon, Technion I.I.T.
%


k = strfind(jpeg_file, '.jpg');

if (~isempty(k))
    rgb = imread(jpeg_file);
    image(rgb);
    axis equal;
    set(gca, 'YTick', []);set(gca, 'XTick', []);
    saveas(gcf,[jpeg_file(1:end-4) '.eps'],'epsc');
else
    error('File is not JPEG, or does not include .jpg suffix.');
end



