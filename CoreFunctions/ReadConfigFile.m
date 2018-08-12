function [ GenParams, MovieParams, PSFParams, WienerParams, AlgParams ] = ReadConfigFile( InputConfigFile, varargin )
% READCONFIGFILE - This function read a TXT file which containts all of the
% configuration definitions for the SRFM script.
%
% Syntax:
% -------
% [ FullMovieLength, MovieBlockLength, GenParams, MovieParams, PSFParams, WienerParams, AlgParams ] = ReadConfigFile( fid )
%
% Inputs:
% -------
% InputConfigFile - Name of configuration TXT file
% inFolder        - Directory in which the configuration file is located (optional)
%
% Outputs:
% --------
% GenParams       - General parameters
% MovieParams     - Movie parameters
% PSFParams       - PSF parameters
% WienerParams    - Filtering parameters
% AlgParams       - Solver parameters
%
% Ver 1. Written by Oren Solomon, Technion I.I>T. 14-12-2015
%

%% Initializations
% ------------------------------------------------------------------------------------------------------------------
% InputConfigFile directory
if nargin > 1
    inFolder = varargin{1};
else
    inFolder = '.';
end

% Open file
fid = fopen(fullfile(inFolder, InputConfigFile), 'r');

%% Parse each line in InputConfigFile
% ------------------------------------------------------------------------------------------------------------------
% Global parameters
MovieParams.FullMovieLength           = str2num(SearchLine(fid, 'MovieParams.FullMovieLength'));           % Total number of frames in the movie       
MovieParams.MovieBlockLength          = str2num(SearchLine(fid, 'MovieParams.MovieBlockLength'));          % Divide the movie into blocks with this length

% Post processing parameters
GenParams.PostProc.Method             = SearchLine(fid, 'GenParams.PostProc.Method');                      % Post processing parameters of the SR image      
GenParams.PostProc.Sigma              = str2num(SearchLine(fid, 'GenParams.PostProc.Sigma'));
GenParams.PostProc.ColorMap           = SearchLine(fid, 'GenParams.PostProc.ColorMap');
GenParams.PostProc.ShowImage          = str2num(SearchLine(fid, 'GenParams.PostProc.ShowImage'));

% Patch related parameters
GenParams.Patches.DoPatches           = str2num(SearchLine(fid, 'GenParams.Patches.DoPatches'));
GenParams.Patches.PatchDimensions.X   = str2num(SearchLine(fid, 'GenParams.Patches.PatchDimensions.X'));   % XY size of a single patch
GenParams.Patches.PatchDimensions.Y   = str2num(SearchLine(fid, 'GenParams.Patches.PatchDimensions.Y'));
GenParams.Patches.Overlap.X           = str2num(SearchLine(fid, 'GenParams.Patches.Overlap.X'));           % Overlap [in pixels] between patches
GenParams.Patches.Overlap.Y           = str2num(SearchLine(fid, 'GenParams.Patches.Overlap.Y'));

% Single patch parameters
MovieParams.Xstart                    = str2num(SearchLine(fid, 'MovieParams.Xstart'));                    % Take a block from the original movie
MovieParams.Ystart                    = str2num(SearchLine(fid, 'MovieParams.Ystart'));
MovieParams.BlockSize                 = str2num(SearchLine(fid, 'MovieParams.BlockSize'));

% Denoising parameters
MovieParams.DoDenoising               = str2num(SearchLine(fid, 'MovieParams.DoDenoising'));               % 1 - Do denoising
MovieParams.DenoisingPer              = str2num(SearchLine(fid, 'MovieParams.DenoisingPer'));              % Discard P percent of the smallest eigenvalues

% Movie sequence parameters
MovieParams.SaveFolder                = SearchLine(fid, 'MovieParams.SaveFolder');                         % Location in which to save the output
MovieParams.InputDataType             = SearchLine(fid, 'MovieParams.InputDataType');                      % Input data type: 'mat', 'txt'
MovieParams.Name                      = SearchLine(fid, 'MovieParams.Name'); 
MovieParams.Path                      = SearchLine(fid, 'MovieParams.Path');
MovieParams.StartFrame                = str2num(SearchLine(fid, 'MovieParams.StartFrame'));
MovieParams.EndFrame                  = str2num(SearchLine(fid, 'MovieParams.EndFrame'));
MovieParams.SliceSize                 = str2num(SearchLine(fid, 'MovieParams.SliceSize'));                 % Extrude a slice of SliceSize X SliceSize pixels from the center of the F.T. of the image (most of the information)

% PSF parameters
PSFParams.PSF_Name                    = SearchLine(fid, 'PSFParams.PSF_Name');                             % Name of a mat file containing a 2D PSF
try
    PSFParams.PSF_Path                = eval(SearchLine(fid, 'PSFParams.PSF_Path'));                       % Location of the PSF
catch
    PSFParams.PSF_Path                = SearchLine(fid, 'PSFParams.PSF_Path');                             % Location of the PSF
end

% Wiener prefilter parameters
WienerParams.NSR                      = str2num(SearchLine(fid, 'WienerParams.NSR'));                      % Ratio beteen noise and PSF magnitudes. Larger values for noisy measurements

% Restoration algorithm parameters
AlgParams.AlgType                     = SearchLine(fid, 'AlgParams.AlgType');                              % 'storm' - perform storm like CS based localization per frame. 'corr' - perform correlations
AlgParams.CSFactor                    = str2num(SearchLine(fid, 'AlgParams.CSFactor'));                    % Emitter's grid resolution compared to the number pixels in each frame

AlgParams.MMVNumVecs                  = str2num(SearchLine(fid, 'AlgParams.MMVNumVecs'));

% single-molecule solvers: MatOMP, MatFISTA
% Correlations solvers: pFISTA variants 
AlgParams.Solver                      = SearchLine(fid, 'AlgParams.Solver');                               % Maximum time-lag for auto-correlation. 1 - variance estimation only

% MatFISTA parameters 
% <><><><><><><><><><>
AlgParams.MatFISTA.IterMax            = str2num(SearchLine(fid, 'AlgParams.MatFISTA.IterMax'));            % Maximum number of iterations
AlgParams.MatFISTA.L0                 = str2num(SearchLine(fid, 'AlgParams.MatFISTA.L0'));                 % Lipschitz constant of the gradient
AlgParams.MatFISTA.Lambda             = str2num(SearchLine(fid, 'AlgParams.MatFISTA.Lambda'));             % Regularization parameter
AlgParams.MatFISTA.NonNegOrth         = str2num(SearchLine(fid, 'AlgParams.MatFISTA.NonNegOrth'));         % Orthogonal projection onto the non-negative orthant. 0 - false, 1 - true
AlgParams.MatFISTA.Beta               = str2num(SearchLine(fid, 'AlgParams.MatFISTA.Beta'));               % In order to avoide numerical problems
AlgParams.MatFISTA.LambdaBar          = str2num(SearchLine(fid, 'AlgParams.MatFISTA.LambdaBar'));          % In order to avoide numerical problems

% MatOMP parameters 
% <><><><><><><><><>
AlgParams.MatOMP.IterMax              = str2num(SearchLine(fid, 'AlgParams.MatOMP.IterMax'));              % Maximum number of iterations
AlgParams.MatOMP.Tol                  = str2num(SearchLine(fid, 'AlgParams.MatOMP.Tol'));                  % Stopping criterion tolerance
AlgParams.MatOMP.NonNegOrth           = str2num(SearchLine(fid, 'AlgParams.MatOMP.NonNegOrth'));           % Orthogonal projection onto the non-negative orthant. 0 - false, 1 - true

% pFISTA_diag parameters
% <><><><><><><><><><><>
AlgParams.pFISTA_diag.Beta            = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag.Beta'));
AlgParams.pFISTA_diag.Lambda          = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag.Lambda'));          % l1 regulariazation parameters - Must be before LambdaBar !!
AlgParams.pFISTA_diag.LambdaBar       = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag.LambdaBar'));
AlgParams.pFISTA_diag.L0              = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag.L0'));              % Lipschitz constant - keep empty
AlgParams.pFISTA_diag.N               = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag.N')); 
AlgParams.pFISTA_diag.IterMax         = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag.IterMax'));
AlgParams.pFISTA_diag.NonNegOrth      = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag.NonNegOrth'));      % Projection onto the non-negative orthant (non-negative x)
AlgParams.pFISTA_diag.LargeScale      = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag.LargeScale'));

% pFISTA_Iterative parameters
% <><><><><><><><><><><><><><>
AlgParams.pFISTA_diag.NumOfSteps      = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag.NumOfSteps')); % Total number of "enveloping" iterations
AlgParams.pFISTA_diag.eps             = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag.eps'));

% pFISTA_diag_TV parameters
% <><><><><><><><><><><><><>
AlgParams.pFISTA_diag_TV.Beta         = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_TV.Beta'));
AlgParams.pFISTA_diag_TV.Lambda       = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_TV.Lambda'));          % l1 regulariazation parameters - Must be before LambdaBar !!
AlgParams.pFISTA_diag_TV.LambdaBar    = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_TV.LambdaBar'));
AlgParams.pFISTA_diag_TV.L0           = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_TV.L0'));              % Lipschitz constant - keep empty
AlgParams.pFISTA_diag_TV.N            = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_TV.N')); 
AlgParams.pFISTA_diag_TV.IterMax      = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_TV.IterMax'));
AlgParams.pFISTA_diag_TV.NonNegOrth   = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_TV.NonNegOrth'));      % Projection onto the non-negative orthant (non-negative x)
AlgParams.pFISTA_diag_TV.LargeScale   = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_TV.LargeScale'));
AlgParams.pFISTA_diag_TV.DenoiseIter  = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_TV.DenoiseIter'));     % Inner loop, number of denoisig iterations
AlgParams.pFISTA_diag_TV.TV           = SearchLine(fid, 'AlgParams.pFISTA_diag_TV.TV');                       % 'iso': Isotropic TV norm, 'l1': Anisotropic TV norm
AlgParams.pFISTA_diag_TV.mon          = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_TV.mon')); 		      % 0 - no monotonicity (regular), 1 - monotone version

% pFISTA_diag_analysis parametes
% <><><><><><><><><><><><><><><>
AlgParams.pFISTA_diag_analysis.Beta         = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_analysis.Beta'));
AlgParams.pFISTA_diag_analysis.Lambda       = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_analysis.Lambda'));          % l1 regulariazation parameters - Must be before LambdaBar !!
AlgParams.pFISTA_diag_analysis.LambdaBar    = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_analysis.LambdaBar'));
AlgParams.pFISTA_diag_analysis.L0           = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_analysis.L0'));              % Lipschitz constant - keep empty
AlgParams.pFISTA_diag_analysis.N            = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_analysis.N')); 
AlgParams.pFISTA_diag_analysis.IterMax      = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_analysis.IterMax'));
AlgParams.pFISTA_diag_analysis.NonNegOrth   = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_analysis.NonNegOrth'));      % Projection onto the non-negative orthant (non-negative x)
AlgParams.pFISTA_diag_analysis.LargeScale   = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_analysis.LargeScale'));
AlgParams.pFISTA_diag_analysis.mu		    = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_analysis.mu'));
AlgParams.pFISTA_diag_analysis.AnalysisType = SearchLine(fid, 'AlgParams.pFISTA_diag_analysis.AnalysisType');
AlgParams.pFISTA_diag_analysis.mon          = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_analysis.mon')); 		      % 0 - no monotonicity (regular), 1 - monotone version

% pFISTA_diag_BM3D parametes
% <><><><><><><><><><><><><>
AlgParams.pFISTA_diag_BM3D.Beta            = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_BM3D.Beta'));
AlgParams.pFISTA_diag_BM3D.Lambda          = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_BM3D.Lambda'));          % l1 regulariazation parameters - Must be before LambdaBar !!
AlgParams.pFISTA_diag_BM3D.LambdaBar       = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_BM3D.LambdaBar'));
AlgParams.pFISTA_diag_BM3D.L0              = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_BM3D.L0'));              % Lipschitz constant - keep empty
AlgParams.pFISTA_diag_BM3D.N               = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_BM3D.N')); 
AlgParams.pFISTA_diag_BM3D.IterMax         = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_BM3D.IterMax'));
AlgParams.pFISTA_diag_BM3D.NonNegOrth      = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_BM3D.NonNegOrth'));      % Projection onto the non-negative orthant (non-negative x)
AlgParams.pFISTA_diag_BM3D.LargeScale      = str2num(SearchLine(fid, 'AlgParams.pFISTA_diag_BM3D.LargeScale'));

% Close file
fclose(fid);

%% Auxiliary functions
% ------------------------------------------------------------------------------------------------------------------
function [ tStr ] = SearchLine(fid, match_str)
% Search for a line in "fid" which begins with "match_str" and read the
% relevant parameter

% Initialize tStr
tStr = [];

% Set cursor to beginning of file
fseek(fid, 0, 'bof');

while ~feof(fid) 
    % Read a line from the file
    tline = fgetl(fid);
    
    % Determine type of line
    if ~isempty(tline)                                                          % Not an empty line
        if tline(1) ~= '%'                                                      % Line is not a comment line
            % Line beginning with '%' is a comment and should be disregarded
            str_ind  = strfind(tline, match_str);
            
            % We found the string that we wanted
            if ~isempty(str_ind)
                ind1 = strfind(tline, '=');
                ind2 = strfind(tline, ';');                                    % Every line must end with ';'

                % Read the value
                tStr = strtrim(tline(ind1(1) + 1:ind2 - 1));
                
                % Clean ' if there are any in the string
                k    = strfind(tStr, '''');
                if ~isempty(k)
                    tStr = tStr(k(1) + 1:k(2) - 1);
                end
                
                % Exit the loop
                break;
            end
        end
    end
end

























