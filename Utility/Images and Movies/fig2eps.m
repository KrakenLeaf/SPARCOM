function fig2eps(FigName)
%FIG2EPS saves the current Fig file to an eps file
%
% SYntax:
% -------
% fig2eps( FigName )
%
% Input:
% ------
% FigName - Desired name for the eps file
%
% Ver 1. Written by Oren Solomon, Technion I.I.T.
%

saveas(gca, [FigName '.eps'], 'epsc');

