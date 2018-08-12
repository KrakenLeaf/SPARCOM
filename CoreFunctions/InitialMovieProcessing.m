function [ MovieData ] = InitialMovieProcessing( GenParams, MovieParams, AlgParams, StartFrame, EndFrame )
%INITIALMOVIEPROCESSING Summary of this function goes here
%   Detailed explanation goes here
%
%

global VERBOSE;

% -------------------------------------------------------------------------------------------------------------------------------------------------------
%% Input checks
% -------------------------------------------------------------------------------------------------------------------------------------------------------
if (AlgParams.CSFactor < 1 && ~isinteger(AlgParams.CSFactor))
    error('SRFM: AlgParams.CSFactor must be an integer equal or greater than 1.');
end

if (MovieParams.SliceSize > MovieParams.BlockSize)
   error('SRFM: MovieParams.BlockSize should be less than or equal MovieParams.BlockSize.'); 
end

if VERBOSE
    disp('------------------------------------');
    disp(['Reconstruction method: ' AlgParams.AlgType]);
    disp(' ');
    disp(['Patch starts at: ' num2str(MovieParams.Xstart) ' ,' num2str(MovieParams.Ystart)]);
    disp('------------------------------------');
    disp(' ');
end

% -------------------------------------------------------------------------------------------------------------------------------------------------------
%% Load data
% -------------------------------------------------------------------------------------------------------------------------------------------------------
if VERBOSE; tic;fprintf('Loading movie data...');end;
% Determine type of input file
switch lower(MovieParams.InputDataType)
    case 'mat'
        MovieDataTmp = load(fullfile(MovieParams.Path, [MovieParams.Name '.mat']));
        MovieData.Movie = MovieDataTmp.Movie(:,:,StartFrame:EndFrame);
    case 'txt'
        MovieData.Movie = txt2mat(fullfile(MovieParams.Path, MovieParams.Name), MovieParams.BlockSize, MovieParams.BlockSize, StartFrame, EndFrame);
    case 'tif'
        MovieData.Movie = tif2mat(fullfile(MovieParams.Path, MovieParams.Name));
        MovieData.Movie = MovieData.Movie(:, :, StartFrame:EndFrame);
    otherwise
        error('SRFM: Input data type not supported');
end

if VERBOSE; fprintf('Done.\n');toc;disp(' ');end;

% Normalize movie intensity to be between [0, 256]
MovieData.Movie = MovieData.Movie/max(max(max(MovieData.Movie)))*256;

% perform denoising if needed
% <><><><><><><><><><><><><><>
% SVD filtering
if MovieParams.DoDenoising
    if VERBOSE; tic;fprintf('Beginning SVD based denoising...');end;
    TmpSVD          = SVDfilt(MovieData.Movie, 'ind', MovieParams.DenoisingPer);
    MovieData.Movie = (TmpSVD{2}); % abs
    if VERBOSE; fprintf('Done.\n');toc;disp(' ');end;
end

% Remove temporal median of the movie - to remove constant background
MovieData.Movie  = MovieData.Movie - repmat(median(MovieData.Movie, 3), [1 1 size(MovieData.Movie, 3)]);

MovieData.MaxVal = max(max(mean(abs(MovieData.Movie), 3)));

% Take only a fragmet of the full image - "patch"
% <><><><><><><><><><><><><><><><><><><><><><><><>
MovieData.Movie = MovieData.Movie(MovieParams.Xstart:MovieParams.Xstart + MovieParams.BlockSize - 1, MovieParams.Ystart:MovieParams.Ystart + MovieParams.BlockSize - 1,:);

% Expand the patches
% [ Tmp ] = PatchExpander( MovieData.Movie, AlgParams.ExpandPix, AlgParams.ExpandMethod );

% Determine dimentions of the movie
% <><><><><><><><><><><><><><><><><>
MovieData.NumOfRows   = length(MovieData.Movie(:,1,1));             % Determine number of rows of the movie
MovieData.NumOfCols   = length(MovieData.Movie(1,:,1));             % Determine number of columns of the movie
MovieData.NumOfFrames = length(MovieData.Movie(1,1,:));             % Determine the length of the movie

% perform denoising if needed
% <><><><><><><><><><><><><><>
if MovieParams.DoDenoising == 1
    if VERBOSE; tic;fprintf('Beginning SVD based denoising...');end;
    MovieData.Movie = SVD_denoise( MovieData.Movie, MovieParams.DenoisingPer );
    if VERBOSE; fprintf('Done.\n');toc;disp(' ');end;
end
