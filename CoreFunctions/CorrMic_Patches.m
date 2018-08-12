function [ MovieData ] = CorrMic_Patches( GenParams, MovieParams, PSFParams, WienerParams, AlgParams, StartFrame, EndFrame )
%CORRMIC_PATCHES - Performs correlation based Super-Resolution Microscopy.
%                  This function performs a patch based analysis.
%
%
%
%
% Written by Oren Solomon, V.1 12-08-2015
%
% -------------------------------------------------------------------------------------------------------------------------------------------------------
%% Description:
% -------------------------------------------------------------------------------------------------------------------------------------------------------
% This script performs Super-Resolution processing on microscopy
% data, using time correlations.

global VERBOSE;

disp('Correlations based Super-Resolution imaging of fluorescence microsopy.');
disp('Ver 2 - 15/07/2015: Patch based analysis');
disp(' ');

% -------------------------------------------------------------------------------------------------------------------------------------------------------
%% Initial processing of the input movie
% -------------------------------------------------------------------------------------------------------------------------------------------------------
[ MovieData ] = InitialMovieProcessing( GenParams, MovieParams, AlgParams, StartFrame, EndFrame );

% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%%                                                                      Patches related code                                                                               %
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
if VERBOSE; tic;fprintf('Extract patches for each frame...');disp(' ');end;

% If we work on the entire image
if ~GenParams.Patches.DoPatches
    GenParams.Patches.PatchDimensions.X = MovieData.NumOfRows;
    GenParams.Patches.PatchDimensions.Y = MovieData.NumOfCols;
end

% Segment the movie into patches
% <><><><><><><><><><><><><><><>
% Prepare patch parameters
FullImageDimensions.X = MovieData.NumOfRows;                        % Full image dimensions
FullImageDimensions.Y = MovieData.NumOfCols;
PatchParams.Type = 'extract';

% Extract patches for each frame in the movie
for ii = 1:MovieData.NumOfFrames
    % NOTE: Every column in MovieData.Patches.Movie for different frame indices is the same patch.
    % Each column is a vectorized patch of the original image.
    [ MovieData.Patches.Movie(:,:,ii), ~ ] = PatchManipulator_V2( MovieData.Movie(:,:,ii), GenParams.Patches.PatchDimensions, ...
        GenParams.Patches.Overlap, FullImageDimensions, PatchParams );
end

% Determine the size of the patches matrix (same for each frame)
[MovieData.Patches.Rows, MovieData.Patches.Cols] = size(MovieData.Patches.Movie(:,:,1));
if VERBOSE; fprintf('Done.\n');toc;disp(' ');end;
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%%                                                                    End of patches related code                                                                          %
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% Perform Preparations
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
MovieData.NumOfRows = GenParams.Patches.PatchDimensions.X;
MovieData.NumOfCols = GenParams.Patches.PatchDimensions.Y;
[ PSFModel, FourierMat, W_RepMat, H, MidRange ] = FreqPrep( GenParams, MovieParams, PSFParams, WienerParams, AlgParams, MovieData );

% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                                                                    Work on each patch separately                                                                         %
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
for tt = 1:MovieData.Patches.Cols
    % Additional preparations
    % ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    if VERBOSE;disp(' ');
    disp('------------------------------------------');
    disp(['       Working on patch #: ' num2str(tt) '/' num2str(MovieData.Patches.Cols)]);
    disp('------------------------------------------');
    disp(' ');end;
    
    % Update the current movie to the relevant patch
    CurrentPatch = reshape(MovieData.Patches.Movie(:,tt,:), GenParams.Patches.PatchDimensions.X, GenParams.Patches.PatchDimensions.Y, MovieData.NumOfFrames);
    
    % Perform FFT on each frame
    % <><><><><><><><><><><><><>
    tic;if VERBOSE;fprintf('Performing FFT analysis on each frame...');end;
    MovieData.SpatFreqMovie = (1/GenParams.Patches.PatchDimensions.X/GenParams.Patches.PatchDimensions.Y)*fft2(CurrentPatch);
    MovieData.SpatFreqMovie = fftshift(MovieData.SpatFreqMovie);
    if VERBOSE; fprintf('Done.\n');toc;disp(' ');end;
    
    % For each frame, perform Wiener filtering - elementwise multiplication
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    if VERBOSE; tic;fprintf('Performing Wiener filtering on each frame...');end;
    if WienerParams.NSR >= 0
        % Apply Wiener
        TempStack = MovieData.SpatFreqMovie.*W_RepMat;
    else
        % Do not apply Wiener - relevant only for the correlations with vectorization
        TempStack = MovieData.SpatFreqMovie;
    end
    if VERBOSE; fprintf('Done.\n');toc;disp(' ');end;
    
    % Update stack after applying the Wiener prefilter
    MovieData.SpatFreqMovieAfterWiener = TempStack;
    
    % -------------------------------------------------------------------------------------------------------------------------------------------------------
    %% Choose type of restoration algorithm
    % -------------------------------------------------------------------------------------------------------------------------------------------------------
    switch lower(AlgParams.AlgType)
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        % ----                                   Perform STORM type reconstruction                                     ---- %
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        case 'storm'
            % Initialize empty stack
            M = GenParams.Patches.PatchDimensions.X*AlgParams.CSFactor;
            N = GenParams.Patches.PatchDimensions.Y*AlgParams.CSFactor;
            TmpSTORM = zeros(M, N, MovieData.NumOfFrames);
            
            % Choose type of solver
            % <><><><><><><><><<><>
            if VERBOSE; tic;fprintf(['Runnig solver: ' AlgParams.Solver '\r\n']);end;
            switch lower(AlgParams.Solver)
                % pFISTA_diag
                % <><><><><><>
                case 'pfista_diag'
                    for ii = 1:MovieData.NumOfFrames      % Perform on each frame separately
                        AlgParams.pFISTA_diag.N = (GenParams.Patches.PatchDimensions.X*AlgParams.CSFactor)^2;
                        TmpSTORM(:, :, ii) = reshape(pFISTA_diag_US(TempStack(:, :, ii), H, AlgParams.pFISTA_diag), sqrt(AlgParams.pFISTA_diag.N), sqrt(AlgParams.pFISTA_diag.N));
                        TmpSTORM(:, :, ii) = Segmentation( TmpSTORM(:, :, ii), AlgParams.SegThreshold );  % 0.75
                        
                        if (GenParams.DebugMode.ShowCompMovie)
                            imshowpair(TmpSTORM(:, :, ii),imresize(MovieData.Movie(:,:,ii), AlgParams.CSFactor),'montage');
                            title('Super resolved image and blurred image');
                            pause(1/GenParams.DebugMode.FrameRate);
                        end
                    end
                % Matrix FISTA
                % <><><><><><>
                case 'matfista'
                    for ii = 1:MovieData.NumOfFrames      % Perform on each frame separately
                        TmpSTORM(:, :, ii) = FISTA(TempStack(:, :, ii), FourierMat.Fx, FourierMat.Fy.',N, N, AlgParams.MatFISTA);
                    end
                % Matrix OMP
                % <><><><><>
                case 'matomp'
                    for ii = 1:MovieData.NumOfFrames      % Perform on each frame separately
                        [TmpSTORM(:, :, ii), Supp] = MatOMP(TempStack(:, :, ii), FourierMat.Fx, conj(FourierMat.Fy), AlgParams.MatOMP);
                    end
                otherwise
                    error('Solver not supported for STORM type method.');
            end
            
            MovieData.SR_Stack = TmpSTORM;
            MovieData.SR_image = sum(TmpSTORM, 3);  % Sum each pixel over the frames dimension
            toc;fprintf('Done.\n');
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
            % ----                               Perform CORRELATIONS type restoration                                     ---- %
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        case 'corr'
            % Vectorize the measurements and the partial Fourier matrices
            % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            if VERBOSE; tic;fprintf('Vectorizing measurements...');end;
            % Vectorize measurements (after Wiener) - vectorization is performed per frame
            TmpMat = zeros(MovieParams.SliceSize*MovieParams.SliceSize, MovieData.NumOfFrames);
            for ii = 1:MovieData.NumOfFrames
                % Vectorized frame X Time index
                TmpMat(:,ii) = reshape(TempStack(:,:,ii), MovieParams.SliceSize*MovieParams.SliceSize, 1);
            end
            MovieData.VecAfterWiener = TmpMat;
            if VERBOSE; fprintf('Done.\n');toc;disp(' ');end;
            
            % -------------------------------------------------------------------------------------------------------------------------------------------------------
            %% Time correlations calculation
            % -------------------------------------------------------------------------------------------------------------------------------------------------------
            % Estimate correlations - NOTE: See if we can make this computation more efficient
            % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>   
            if VERBOSE; tic;fprintf('Performing correlations analysis...');end;
            MovieData.CTF_Frame = CalcCorr_Mat_V3( MovieData.VecAfterWiener, [(MovieParams.SliceSize)^2, MovieData.NumOfFrames], AlgParams.MMVNumVecs );

            if VERBOSE; fprintf('Done.\n');toc;disp(' ');end;
            
            % Solve MMV problem: Perform support estimation
            % <><><><><><><><><><><><><><><><><><><><><><><>
            if VERBOSE; tic;disp(['Running solver: ' AlgParams.Solver]);
                disp('------------------------------------'); end;
            
            % The patch is considered to contain data
            switch lower(AlgParams.Solver)
                case 'pfista_diag'
                    AlgParams.pFISTA_diag.N = (GenParams.Patches.PatchDimensions.X*AlgParams.CSFactor)^2;
                    [ X ]    = pFISTA_diag( MovieData.CTF_Frame, H, AlgParams.pFISTA_diag );
                case 'pfista_iterative'
                    AlgParams.pFISTA_diag.N = (GenParams.Patches.PatchDimensions.X*AlgParams.CSFactor)^2;
                    [ X, ~ ] = pFISTA_Iterative( MovieData.CTF_Frame, H, AlgParams.pFISTA_diag );
                case 'pfista_diag_tv'
                    AlgParams.pFISTA_diag_TV.N = (GenParams.Patches.PatchDimensions.X*AlgParams.CSFactor)^2;
                    [ X ] = pFISTA_diag_TV( MovieData.CTF_Frame, H, AlgParams.pFISTA_diag_TV );
                case 'pfista_diag_analysis'
                    AlgParams.pFISTA_diag_analysis.N = (MovieData.NumOfCols*AlgParams.CSFactor)^2;
                    [ X ] = pFISTA_diag_analysis( MovieData.CTF_Frame, H, AlgParams.pFISTA_diag_analysis );
                case 'pfista_diag_bm3d'
                    AlgParams.pFISTA_diag_BM3D.N = (MovieData.NumOfCols*AlgParams.CSFactor)^2;
                    [ X ] = pFISTA_diag_BM3D( MovieData.CTF_Frame, H, AlgParams.pFISTA_diag_BM3D );
                otherwise
                    error('Solver not supported for CORRELATIONS type method.');
            end
            
            % Return to un-normalized solution
            if ~(strcmpi(AlgParams.Solver, 'fpg_iterative') || strcmpi(AlgParams.Solver, 'pfista_diag') ||...
                    strcmpi(AlgParams.Solver, 'pfista_iterative') || strcmpi(AlgParams.Solver, 'pfista_diag_TV') || strcmpi(AlgParams.Solver, 'pfista_diag_analysis')  || strcmpi(AlgParams.Solver, 'pfista_diag_bm3d'))
                if VERBOSE; tic;disp(' ');disp('------------------------------------');fprintf('Un-normalizing solution...');end;
                U = FourierMat.NrmMat*(diag(sum(X, 2))*FourierMat.NrmMat');
                
                % Super-resolution image - U is a diagonal matrix
                MovieData.SR_image = reshape(diag(U), MovieParams.BlockSize*AlgParams.CSFactor, MovieParams.BlockSize*AlgParams.CSFactor);
                if VERBOSE; fprintf('Done.\n');toc;disp(' ');end;
            else
                % sum(X, 2) is a vector
                MovieData.SR_image = reshape(sum(X, 2), GenParams.Patches.PatchDimensions.X*AlgParams.CSFactor, GenParams.Patches.PatchDimensions.X*AlgParams.CSFactor);
            end
                        
            if VERBOSE; fprintf('Done.\n');toc;disp(' ');end;
        otherwise
            error('Reconstruction method not supported.');
    end
    % Store HR patches in matrix
    % <><><><><><><><><><><><><>
    MovieData.Patches.HR_Patch(:,tt) = MovieData.SR_image(:);
end

% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                                              Combine all High-Resolution patches to a single HR image                                                                    %
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
% Combine parameters
PatchParams.Type = 'combine';
PatchDimensions.X = GenParams.Patches.PatchDimensions.X*AlgParams.CSFactor;
PatchDimensions.Y = GenParams.Patches.PatchDimensions.Y*AlgParams.CSFactor;
Dimensions.X = FullImageDimensions.X*AlgParams.CSFactor;
Dimensions.Y = FullImageDimensions.Y*AlgParams.CSFactor;
Overlap.X = GenParams.Patches.Overlap.X*AlgParams.CSFactor;
Overlap.Y = GenParams.Patches.Overlap.Y*AlgParams.CSFactor;
PatchParams.FFT = 0;

% Combine
[ MovieData.SR_image, ~ ] = PatchManipulator_V2( MovieData.Patches.HR_Patch, PatchDimensions, Overlap, Dimensions, PatchParams );

% -------------------------------------------------------------------------------------------------------------------------------------------------------
%% Vizualization
% -------------------------------------------------------------------------------------------------------------------------------------------------------
MovieData.SR_image_after_ps = ImagePostProcessing( MovieData.SR_image, GenParams.PostProc.Method, GenParams.PostProc );

% Convolve with the low resolution PSF
MovieData.SR_image_after_ps = imfilter(MovieData.SR_image_after_ps, fftshift(PSFModel.PSF).^1);




















