% *************************************************************************************************************************************
% Super-Resolution Fluorescence Microscopy (SPARCOM)
%
% Description:
% ------------
% Two methods are considered: 1. single-molecule like algorithm 
%                             2. Correlations based algorithm
%
% This is an envelope script which runs the main functions:
%                             1. CorrMic         - single patch analysis
%                             2. CorrMic_Patches - patch based analysis
%
% Written by Oren Solomon, Technion , I.I.T.
% Ver 1
% *************************************************************************************************************************************
clc;
clear;
close all;

% Add relevant folders - can be performed only once
AddPathList;

% Measure total time
TotalTime = tic;

% -------------------------------------------------------------------------------------------------------------------------------------------------------
%% Parameters
% -------------------------------------------------------------------------------------------------------------------------------------------------------
global VERBOSE SAVE_FOLDER_PREFIX; 

VERBOSE = 1;

InternalSaveFlag = 1;   % 1 - save, 0 don't save

% Read configuration TXT file
InputConfigFile = 'SRFM_Config_1.txt';
[ GenParams, MovieParams, PSFParams, WienerParams, AlgParams ] = ReadConfigFile( InputConfigFile );

% -------------------------------------------------------------------------------------------------------------------------------------------------------
%% Create save folder
% -------------------------------------------------------------------------------------------------------------------------------------------------------
if InternalSaveFlag
    % Save Output - Create a new directory inside SaveFolder
    try 
        DestFolder = [SAVE_FOLDER_PREFIX date];
    catch
        DestFolder = ['SRFM_Results_' date];
    end
    if ~exist(fullfile(MovieParams.SaveFolder, DestFolder), 'dir')
        mkdir(MovieParams.SaveFolder, DestFolder);
    end
    
    % Current time stamp - Time stamp corresponds to the beginning of execution
    TimeStamp = datestr(clock);
    TimeStamp(ismember(TimeStamp, ' -:')) = ['_'];
    
    % Copy the configuration file to the specified folder
    copyfile(InputConfigFile, fullfile(MovieParams.SaveFolder, DestFolder, [InputConfigFile(1:end-4) '_' TimeStamp '.txt']));
end

% -------------------------------------------------------------------------------------------------------------------------------------------------------
%% Run the algorithm
% -------------------------------------------------------------------------------------------------------------------------------------------------------
if VERBOSE
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
    disp('SPARCOM: Sparsity-based Super-Resolution Correlations Microscopy');
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
    disp(' ');
    disp(['Data folder: ' MovieParams.Path]);
    disp(['SRF        = ' num2str(AlgParams.CSFactor)]);
    disp(['Block size = ' num2str(MovieParams.BlockSize)]);
    disp(' ');
    disp(' ');
end

% Divide the movie into blocks
DivBlocks = [1:MovieParams.MovieBlockLength:MovieParams.FullMovieLength MovieParams.FullMovieLength + 1];

% Work on each movie block separately
for ii = 1:length(DivBlocks) - 1
    % Display
    % <><><><>
    disp('*** ------------------- *** ------------------- ***');
    disp(['          Working on frames: ' num2str(DivBlocks(ii)) ' - ' num2str(DivBlocks(ii + 1) - 1)]);
    disp('*** ------------------- *** ------------------- ***'); disp(' ');
    
    % Perform patch based analysis
    % <><><><><><><><><><><><><><>
    if GenParams.Patches.DoPatches
        % Apply algorithm
        [ MovieData ] = CorrMic_Patches( GenParams, MovieParams, PSFParams, WienerParams, AlgParams, DivBlocks(ii),  DivBlocks(ii + 1) - 1 );
        
        % Accumulate the HR patches for later use
        HR_Patches{ii} = MovieData.Patches.HR_Patch;
        
    % Perform single patch analysis
    % <><><><><><><><><><><><><><><>
    else
        % Apply algorithm
        [ MovieData ] = CorrMic( GenParams, MovieParams, PSFParams, WienerParams, AlgParams, DivBlocks(ii),  DivBlocks(ii + 1) - 1 );  
    end
    
    % Accumulate the SR image from the current block
    SR_image_stack(:, :, ii) = MovieData.SR_image_after_ps;
end

% Construct super-resolution image
SR_image = sum(SR_image_stack, 3);

% -------------------------------------------------------------------------------------------------------------------------------------------------------
%% Construct SR image and save in desired location
% -------------------------------------------------------------------------------------------------------------------------------------------------------
if InternalSaveFlag
    if VERBOSE; disp(['Saving results in folder: ' MovieParams.SaveFolder]); end;
    
    % Save SR results
    save([fullfile(MovieParams.SaveFolder, DestFolder) '/SR_image_stack_' TimeStamp '.mat'], 'SR_image_stack');
    save([fullfile(MovieParams.SaveFolder, DestFolder) '/SR_image_' TimeStamp '.mat'], 'SR_image');
end

% Display total time
disp(['Total processing time: ' num2str(toc(TotalTime)) ' seconds.']);

%% Comparison - not part of the processing
% -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load(fullfile(MovieParams.Path, PSFParams.PSF_Name));

Ref = tif2mat(fullfile(MovieParams.Path, MovieParams.Name));
Ref = Ref(MovieParams.Xstart:MovieParams.Xstart + MovieParams.BlockSize - 1, MovieParams.Ystart:MovieParams.Ystart + MovieParams.BlockSize - 1, :);

figure;
h(1) = subplot(1,2,1);
imagesc(imresize(sum(Ref, 3), AlgParams.CSFactor)); title('Diffraction Limited'); axis image;colormap(gray.^0.4);
h(2) = subplot(1,2,2);
imagesc(SR_image); title('SR reconstruction'); axis image;colormap pink;
linkaxes(h, 'xy');

disp('Done processing.');


















