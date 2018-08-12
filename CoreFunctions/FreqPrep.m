function [ PSFModel, FourierMat, W_RepMat, H, MidRange ] = FreqPrep( GenParams, MovieParams, PSFParams, WienerParams, AlgParams, MovieData )
%FREQPREP - Prepare matrices needed for the SRFM process

global VERBOSE;

% -------------------------------------------------------------------------------------------------------------------------------------------------------
%% Load PSF model to memory
% -------------------------------------------------------------------------------------------------------------------------------------------------------
if VERBOSE; tic;fprintf('Loading PSF to memory...');end;
load(fullfile(PSFParams.PSF_Path,PSFParams.PSF_Name));
eval(['PSFModel.PSF = ifftshift(' PSFParams.PSF_Name ');']);          % PSF with ifftshift - this is actually the centered PSF in MATLAB coordinates
PSFModel.PSF = PSFModel.PSF/max(max(PSFModel.PSF)); % Normalize PSF to values between [0, 1]

% Take only a fragmet of the PSF - "patch"
% <><><><><><><><><><><><><><><><><><><><>
MidRangePSF = floor(MovieData.NumOfRows/2);
Blk1_PSF = PSFModel.PSF(1:MidRangePSF + mod(MovieData.NumOfRows,2), 1:MidRangePSF + mod(MovieData.NumOfRows,2),:);
Blk2_PSF = PSFModel.PSF(1:MidRangePSF + mod(MovieData.NumOfRows,2), end - MidRangePSF + 1:end,:);
Blk3_PSF = PSFModel.PSF(end - MidRangePSF + 1:end, 1:MidRangePSF + mod(MovieData.NumOfRows,2),:);
Blk4_PSF = PSFModel.PSF(end - MidRangePSF + 1:end, end - MidRangePSF + 1:end,:);
PSFModel.PSF = [Blk1_PSF Blk2_PSF;Blk3_PSF Blk4_PSF];

% -------------------------------------------------------------------------------------------------------------------------------------------------------
%% Spatial frequency calculations
% -------------------------------------------------------------------------------------------------------------------------------------------------------
% Perform FFT on the PSF - Note that the PSF is not centralized
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
PSFModel.PSF_FFT       = (1/MovieData.NumOfRows/MovieData.NumOfCols)*fft2(PSFModel.PSF); % Correlations mode
PSFModel.PSF_FFT       = fftshift(PSFModel.PSF_FFT);
PSFModel.PSF_Magnitude = abs(PSFModel.PSF_FFT);                                          % Magnitude of the PSF
PSFModel.PSF_Phase     = angle(PSFModel.PSF_FFT);                                        % Phase of the PSF

if VERBOSE; fprintf('Done.\n');toc;disp(' ');end;

% Create the Wiener prefilter to the measurements in the spatial frequency domain
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
if VERBOSE; tic;fprintf('Preparing the Wiener filter or H matrix...');end;
MidRange = floor(MovieParams.SliceSize/2);

% Do not do Wiener filtering - if not required
% --------------------------------------------
if WienerParams.NSR < 0
    Blk1 = PSFModel.PSF_FFT(1:MidRange + mod(MovieParams.SliceSize,2), 1:MidRange + mod(MovieParams.SliceSize,2),:);
    Blk2 = PSFModel.PSF_FFT(1:MidRange + mod(MovieParams.SliceSize,2), end - MidRange + 1:end,:);
    Blk3 = PSFModel.PSF_FFT(end - MidRange + 1:end, 1:MidRange + mod(MovieParams.SliceSize,2),:);
    Blk4 = PSFModel.PSF_FFT(end - MidRange + 1:end, end - MidRange + 1:end,:);
    PSFModel.PSF_FFT = [Blk1 Blk2;Blk3 Blk4];
    
    % Perform column normalization on A
    Nfactor = norm(PSFModel.PSF_FFT(:), 2);
    
    % Prepare the PSF_FFT matrix - if no Wiener filter was applied
    H = spdiags(PSFModel.PSF_FFT(:), 0, sparse(size(PSFModel.PSF_FFT, 1)^2, size(PSFModel.PSF_FFT, 2)^2));              % No Nfactor 
%     H = spdiags(PSFModel.PSF_FFT(:)/Nfactor, 0, sparse(size(PSFModel.PSF_FFT, 1)^2, size(PSFModel.PSF_FFT, 2)^2));
    
    % No Wiener filter
    W_RepMat = [];
else
% Prepare the Wiener filter
% -------------------------
    W = conj(PSFModel.PSF_FFT)./((PSFModel.PSF_Magnitude).^2 + WienerParams.NSR);           % The filter is the same for every frame in the movie. This is a matrix
    W_RepMat = repmat(W,1,1,MovieData.NumOfFrames);                                         % No need for a FOR loop    
    
    H = [];
end
if VERBOSE; fprintf('Done.\n');toc;disp(' ');end;

% Construct partial Fourier matrices - if needed
% <><><><><><><><><><><><><><><><><><><><><><><>
if ~(strcmpi(AlgParams.Solver, 'pFISTA_diag') || strcmpi(AlgParams.Solver, 'pFISTA_Iterative') ||...
        strcmpi(AlgParams.Solver, 'pFISTA_diag_TV') ||  strcmpi(AlgParams.Solver, 'pFISTA_diag_analysis') ||  strcmpi(AlgParams.Solver, 'pFISTA_diag_BM3D'))
    if VERBOSE; tic;fprintf('Constructing partial Fourier matrices...');end;
    
    % Determine Fourier matrices size
    if ~GenParams.Patches.DoPatches
        % Work on entire image
        Xdim = MovieData.NumOfRows;
        Ydim = MovieData.NumOfCols;
    else
        % Work on patches
        Xdim = GenParams.Patches.PatchDimensions.X;
        Ydim = GenParams.Patches.PatchDimensions.Y;
    end
    
    % Full Fourier matrices
    Fx = dftmtx(Xdim*AlgParams.CSFactor);
    Fy = dftmtx(Ydim*AlgParams.CSFactor);
    
    if (MovieParams.SliceSize < 0 || isempty(MovieParams.SliceSize))
        MovieParams.SliceSize = GenParams.Patches.PatchDimensions.X;
        %     MovieParams.SliceSize = GenParams.Patches.PatchDimensions.Y;                      % Just as a reminder
    end
    
    FourierMat.Fx = fftshift([Fx(1:MidRange + mod(MovieParams.SliceSize,2),:); Fx(end - MidRange + 1:end,:)], 1);
    FourierMat.Fy = fftshift([Fy(1:MidRange + mod(MovieParams.SliceSize,2),:); Fy(end - MidRange + 1:end,:)], 1);
    if VERBOSE; fprintf('Done.\n');toc;disp(' ');end;
    
    % Cosntruct the measurement matrix if we work in correlation mode
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    if strcmpi(AlgParams.AlgType,'corr')
        % Construct the measurement matrix
        % <><><><><><><><><><><><><><><><>
        if VERBOSE; tic;fprintf('Performing Kronecker product...');end;
        FourierMat.Kron = kron(FourierMat.Fy, FourierMat.Fx);                                   % Perform Kronecker product calculation on the partial Fourier matrices
        
        % Multiply by H if no Wiener was applied
        if WienerParams.NSR < 0
            FourierMat.Kron = H*FourierMat.Kron;
        end
        
        % Normalize the Kronecker matrix
        if ~strcmpi(AlgParams.Solver, 'fpg_iterative')
            [ FourierMat.Kron, FourierMat.NrmMat ] = NormMat( FourierMat.Kron, 'l2' );
        end
        if VERBOSE; fprintf('Done.\n');toc;disp(' ');end;
    end
else
    FourierMat = [];
end
