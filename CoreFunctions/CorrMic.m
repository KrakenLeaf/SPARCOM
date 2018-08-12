function [ MovieData ] = CorrMic( GenParams, MovieParams, PSFParams, WienerParams, AlgParams, StartFrame, EndFrame )
%CORRMIC - Performs correlation based Super-Resolution Microscopy.
%          This function performs the analysis on a SINGLE patch.
%
%
%
%
% Written by Oren Solomon,
% V.1 12-08-2015: Basic formulation
% V.2 29-10-2015: Matrix formulation
%
% -------------------------------------------------------------------------------------------------------------------------------------------------------
%% Description:
% -------------------------------------------------------------------------------------------------------------------------------------------------------
% This script performs Super-Resolution processing on microscopy
% data, using time correlations.

global VERBOSE;

disp('Correlations based Super-Resolution imaging of fluorescence microsopy.');
disp('**********************************************************************');
disp(' ');
disp('Ver 1 - 15/07/2015: Initial formulation');
disp('Ver 2 - 01/10/2015: Matrix formulation of the correlations problem');
disp(' ');
disp('**********************************************************************');
disp(' ');

% -------------------------------------------------------------------------------------------------------------------------------------------------------
%% Initial processing of the input movie
% -------------------------------------------------------------------------------------------------------------------------------------------------------
[ MovieData ] = InitialMovieProcessing( GenParams, MovieParams, AlgParams, StartFrame, EndFrame );

if strcmpi(AlgParams.AlgType, 'corr') || strcmpi(AlgParams.AlgType, 'storm')
    % -------------------------------------------------------------------------------------------------------------------------------------------------------
    %% Spatial frequency calculations
    % -------------------------------------------------------------------------------------------------------------------------------------------------------
    % Perform frequency preparations
    % <><><><><><><><><><><><><><><>
    [ PSFModel, FourierMat, W_RepMat, H, MidRange ] = FreqPrep( GenParams, MovieParams, PSFParams, WienerParams, AlgParams, MovieData );
    
    % Perform FFT on each frame
    % <><><><><><><><><><><><><>
    if VERBOSE; tic;fprintf('Performing FFT analysis on each frame...');end;
    MovieData.SpatFreqMovie = (1/MovieData.NumOfRows/MovieData.NumOfCols)*fft2(MovieData.Movie);
    MovieData.SpatFreqMovie = fftshift(MovieData.SpatFreqMovie);
%     MovieData.SpatFreqMovie = fftshift(MovieData.Movie); % A test modification - if sampling directly in Fourier space
    if VERBOSE; fprintf('Done.\n');toc;disp(' ');end;
    
    % For each frame, perform Wiener filtering - elementwise multiplication
    if WienerParams.NSR >= 0
        % Apply Wiener
        if VERBOSE; tic;fprintf('Performing Wiener filtering on each frame...');end;
        TempStack = MovieData.SpatFreqMovie.*W_RepMat;
        if VERBOSE; fprintf('Done.\n');toc;disp(' ');end;
    else
        % Do not apply Wiener - relevant only for the correlations with vectorization
        TempStack = MovieData.SpatFreqMovie;
    end
    
    % Update stack after applying the Wiener prefilter
    MovieData.SpatFreqMovieAfterWiener = TempStack;
elseif strcmpi(AlgParams.AlgType, 'space_corr')
    % -------------------------------------------------------------------------------------------------------------------------------------------------------
    %% Space domain calculations
    % -------------------------------------------------------------------------------------------------------------------------------------------------------
    TempStack = MovieData.Movie;
    
    % Arrange the parameters - NEEDS TO BE CHANGED IN THE FUTURE
    sPSFname = 'SOFI';
    
    sPSFParams{1} = 1;                                        % m
    sPSFParams{2} = 0;                                        % theta      [rad]
    sPSFParams{3} = 16;                                       % PixelSize  [nm]
    sPSFParams{4} = 625;                                      % lamda      [nm]
    sPSFParams{5} = 1.4;                                      % NA
    sPSFParams{6} = AlgParams.CSFactor*MovieParams.SliceSize; % PSFsize
    sPSFParams{7} = AlgParams.CSFactor;                       % BinSize
    
    % Create the dictionary
    % <><><><><><><><><><><>
    if VERBOSE; tic;fprintf('Building the dictionary...');end;
    [ Udic ] = GenPSFDic(sPSFname, sPSFParams );
    
    % Normalize columns of Udic
    if ~strcmpi(AlgParams.Solver, 'fpg_iterative')
        [ Udic, UdicNrmMat ] = NormMat( Udic, 'l2' );
    end
    if VERBOSE; fprintf('Done.\n');toc;disp(' ');end;
end

% -------------------------------------------------------------------------------------------------------------------------------------------------------
%% Choose type of restoration algorithm
% -------------------------------------------------------------------------------------------------------------------------------------------------------
switch lower(AlgParams.AlgType)
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % ----                                   Perform STORM type restoration                                        ---- %
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    case 'storm'
        % Initialize empty stack
        M = MovieData.NumOfRows*AlgParams.CSFactor;
        N = MovieData.NumOfCols*AlgParams.CSFactor;
        TmpSTORM = zeros(M, N, MovieData.NumOfFrames);
        
        if (GenParams.DebugMode.ShowCompMovie)
            figure;
        end
        % Choose type of solver
        % <><><><><><><><><<><>
        if VERBOSE; tic;fprintf(['Running solver: ' AlgParams.Solver '\r\n']); end;
        switch lower(AlgParams.Solver)
            % pFISTA_diag
            % <><><><><><>
            case 'pfista_diag'
                for ii = 1:MovieData.NumOfFrames      % Perform on each frame separately
                    if VERBOSE; disp(['Frame #: ' num2str(ii) ' / ' num2str(MovieData.NumOfFrames)]); end;
                    if VERBOSE; fprintf(['~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  \r\n']); end;
                    
                    AlgParams.pFISTA_diag.N = (MovieParams.BlockSize*AlgParams.CSFactor)^2;
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
                    [ TmpSTORM(:, :, ii) ] = Segmentation( TmpSTORM(:, :, ii), AlgParams.SegThreshold );  % 0.75
                    
                    if (GenParams.DebugMode.ShowCompMovie)
                        imshowpair(TmpSTORM(:, :, ii),imresize(MovieData.Movie(:,:,ii), AlgParams.CSFactor),'montage');
                        title('Super resolved image and blurred image');
                        pause(1/GenParams.DebugMode.FrameRate);
                    end
                end
                % Matrix OMP
                % <><><><><>
            case 'matomp'
                for ii = 1:MovieData.NumOfFrames      % Perform on each frame separately
                    [TmpSTORM(:, :, ii), Supp] = MatOMP(TempStack(:, :, ii), FourierMat.Fx, conj(FourierMat.Fy), AlgParams.MatOMP);
                    [ TmpSTORM(:, :, ii) ] = Segmentation( TmpSTORM(:, :, ii), AlgParams.SegThreshold );  % 0.75
                    
                    if (GenParams.DebugMode.ShowCompMovie)
                        imshowpair(TmpSTORM(:, :, ii),imresize(MovieData.Movie(:,:,ii), AlgParams.CSFactor),'montage');
                        title('Super resolved image and blurred image');
                        pause(1/GenParams.DebugMode.FrameRate);
                    end
                end
            otherwise
                error('Solver not supported for STORM type method.');
        end
        
        MovieData.SR_Stack = TmpSTORM;
        MovieData.SR_image = sum(TmpSTORM, 3);  % Sum each pixel over the frames dimension
        if VERBOSE; toc;fprintf('Done.\n');disp(' ');end;
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
        %         MovieData.CTF_Frame = CalcCorr_Mat_V1( MovieData.VecAfterWiener, [(MovieParams.SliceSize)^2, MovieData.NumOfFrames], AlgParams.MMVNumVecs, GenParams.parFlag );
        
        % DEBUG
        if (GenParams.DebugMode.ShowImages)
            CorrIndX = 2;
            CorrIndY = 4;
            figure;
            subplot(2,1,1);
            plot(1:length(MovieData.Corr(1, 1, :)), squeeze(real(MovieData.Corr(CorrIndX, CorrIndY, :))), '-*');grid on;
            ylabel('Real');title(['Real and imaginary parts of correlation of pixel [' num2str(CorrIndX) ', ' num2str(CorrIndY) ']']);
            subplot(2,1,2);
            plot(1:length(MovieData.Corr(1, 1, :)), squeeze(imag(MovieData.Corr(CorrIndX, CorrIndY, :))), '-*');grid on;
            xlabel('Time-lag');ylabel('Imaginary');
            
            figure;
            subplot(2,1,1);
            plot(1:length(MovieData.Corr(1, 1, :)), squeeze(abs(MovieData.Corr(CorrIndX, CorrIndY, :))), '-*');grid on;
            ylabel('Magnitude');title(['Magnitude and phase of correlation of pixel [' num2str(CorrIndX) ', ' num2str(CorrIndY) ']']);
            subplot(2,1,2);
            plot(1:length(MovieData.Corr(1, 1, :)), squeeze(angle(MovieData.Corr(CorrIndX, CorrIndY, :))), '-*');grid on;
            xlabel('Time-lag');ylabel('Phase');
        end
        
        if VERBOSE; fprintf('Done.\n');toc;disp(' ');end;
        
        % Solve MMV problem: Perform support estimation
        % <><><><><><><><><><><><><><><><><><><><><><><>
        if VERBOSE; tic;disp(['Running solver: ' AlgParams.Solver]);
            disp('------------------------------------'); end;
        
        % Run your favorite solver
        switch lower(AlgParams.Solver)
            case 'pfista_diag'
                AlgParams.pFISTA_diag.N = (MovieData.NumOfCols*AlgParams.CSFactor)^2;
                [ X ]    = pFISTA_diag( MovieData.CTF_Frame, H, AlgParams.pFISTA_diag );
            case 'pfista_iterative'
                AlgParams.pFISTA_diag.N = (MovieData.NumOfCols*AlgParams.CSFactor)^2;
                [ X, ~ ] = pFISTA_Iterative( MovieData.CTF_Frame, H, AlgParams.pFISTA_diag );
            case 'pfista_diag_tv'
                AlgParams.pFISTA_diag_TV.N = (MovieData.NumOfCols*AlgParams.CSFactor)^2;
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
            MovieData.SR_image = reshape(sum(X, 2), MovieParams.BlockSize*AlgParams.CSFactor, MovieParams.BlockSize*AlgParams.CSFactor);
        end    
        
        if VERBOSE; fprintf('Done.\n');toc;disp(' ');end;
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
        % ----                       Perform CORRELATIONS type restoration in the space domain                         ---- %
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    case 'space_corr'
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
        
        % Run your favorite solver
        switch lower(AlgParams.Solver)
            case 'fpg_affine'
                [ X, L ] = FPG_Affine_AB( MovieData.CTF_Frame, Udic, Udic, AlgParams.WeightsVec, AlgParams.FPG_Affine );
            case 'fpg_iterative'
                [ X, L ] = FPG_Iterative( MovieData.CTF_Frame, Udic, Udic, AlgParams.WeightsVec, AlgParams.FPG_Affine );
            otherwise
                error('Solver not supported for CORRELATIONS type method.');
        end
        
        % Return to un-normalized solution
        if ~(strcmpi(AlgParams.Solver, 'fpg_iterative'))
            if VERBOSE; tic;disp(' ');disp('------------------------------------');fprintf('Un-normalizing solution...');end;
            U = UdicNrmMat*(diag(sum(X, 2))*UdicNrmMat');
            
            % Super-resolution image - U is a diagonal matrix
            MovieData.SR_image = reshape(diag(U), MovieParams.BlockSize*AlgParams.CSFactor, MovieParams.BlockSize*AlgParams.CSFactor);
            if VERBOSE; fprintf('Done.\n');toc;disp(' ');end;
        else
            % sum(X, 2) is a vector
            MovieData.SR_image = reshape(sum(X, 2), MovieParams.BlockSize*AlgParams.CSFactor, MovieParams.BlockSize*AlgParams.CSFactor);
        end
        
        
        if VERBOSE; fprintf('Done.\n');toc;disp(' ');end;
    otherwise
        error('Reconstruction method not supported.');
end

% -------------------------------------------------------------------------------------------------------------------------------------------------------
%% Vizualization and post-processing
% -------------------------------------------------------------------------------------------------------------------------------------------------------
MovieData.SR_image_after_ps = ImagePostProcessing( MovieData.SR_image, GenParams.PostProc.Method, GenParams.PostProc );

% Convolve with the low resolution PSF
MovieData.SR_image_after_ps = imfilter(MovieData.SR_image_after_ps, fftshift(PSFModel.PSF).^1);

























