function [ SucFlag ] = SaveLogger( Data, FigHandle, SaveDir, GenParams, MovieParams, PSFParams, WienerParams, AlgParams)
% SAVELOGGER - This function saves a mat file or a figure and keeps a log file of the data (one for Data and one for FigHandle).
% 
% Syntax:
% -------
% [ SucFlag ] = SaveLogger( Data, FigHandle, SaveDir, GenParams, MovieParams, PSFParams, WienerParams, AlgParams)
%
% Inputs:
% -------
% Data         - Data to save. Use [] if you do not wish to save anything.
% FigHandle    - Figure handle to save. Use [] if you do not wish to save anything.
% SaveDir      - Desired save directory.
% GenParams    - GenParams to write in the log files.
% MovieParams  - MovieParams to write in the log files.
% PSFParams    - PSFParams to write in the log files.
% WienerParams - WienerParams to write in the log files.
% AlgParams    - AlgParams to write in the log files.
%
% Outputs:
% --------
% SucFlag      - 0 if writing to file succeeded, -1 if not.
%
% Written by Oren Solomon,
%                           Ver 1, 16-08-2015

%% Initial flags
% --------------------------------------------------------
DataFlag     = 0;
FigHandleFlag = 0;

% Determine if we want to save a mat file, a figure or both
TL = 0;
if ~isempty(Data)
    DataFlag = 1;              % Save mat file
    TL = TL + 1;
end
if ~isempty(FigHandle)
    FigHandleFlag = 1;          % Save figure
    TL = TL + 1;
end

% If there is no data to save
if FigHandleFlag == 0 && DataFlag == 0
    disp('No data to save.');
    SucFlag = 0;
    return
end

% Determine if SaveDir exists. If not, create it
if ~exist(SaveDir, 'dir')
    mkdir(SaveDir);
end

%% Open log file to determine name of last saved item
% --------------------------------------------------------
fid = fopen([SaveDir '\SRFM_figure_log.txt'], 'a+');          % Open or create new file for reading and writing. Append data to the end of the file 

% Determine if file opening succeeded
if fid == -1, error('SaveLogger: Cannot open file.'), end;

% Move cursor to beginning of the text file
frewind(fid);

% Initial name buffer and cursor in the cell array
ii = 1;
NameBuffer = {};
while ~feof(fid)
   % Read line
   tline = fgetl(fid); 
   
   % Find if this line contains a name
   if strfind(tline, 'NAME')
       NameBuffer{ii} = tline;
       ii = ii + 1;
   end
end
SucFlag = fclose(fid);

% Determine last name used
% --------------------------------------------------------
% Name formats: for a figure   - 'FIG_XXX'
%               for a mat file - 'MAT_XXX'
if ii == 1      % Empty or non-existing log file
    LastIndex = 0;
else
    LastName  = NameBuffer{length(NameBuffer)}(7:end);
    LastIndex = str2num(LastName(5:end));
end

%% Save figure
if FigHandleFlag
    saveas(FigHandle, [SaveDir '\FIG_' num2str(LastIndex + 1) '.fig'], 'fig');
end
if DataFlag
    eval(['save([SaveDir ''\MAT_'' num2str(LastIndex + 1)], ''Data'');']);
end

%% Open log file and write log data
% --------------------------------------------------------
% Write twice if necessary (mat file / figure)
for kk = 1:TL
    % Write name of figure / mat file
    if FigHandleFlag
        fid = fopen([SaveDir '\SRFM_figure_log.txt'], 'a+');          % Open or create new file for reading and writing. Append data to the end of the file
        
        % Determine if file opening succeeded
        if fid == -1 
            SucFlag = fid;
            error('SaveLogger: Cannot open file.')
        end;
        
        fprintf(fid, 'NAME: FIG_%06d\r\n', LastIndex + 1);
        FigHandleFlag = 0;
    else
        fid = fopen([SaveDir '\SRFM_mat_log.txt'], 'a+');          % Open or create new file for reading and writing. Append data to the end of the file
        
        % Determine if file opening succeeded
        if fid == -1
            SucFlag = fid;
            error('SaveLogger: Cannot open file.')
        end;
        
        fprintf(fid, 'NAME: MAT_%06d\r\n', LastIndex + 1);
    end
    fprintf(fid, '\r\n');
    
    % Write parameters of figure / mat file
    % --------------------------------------------------------
    % General movie parameters
    fprintf (fid, 'General movie parameters\r\n');
    fprintf(fid, '************************************************\r\n');
    fprintf(fid, 'MovieParams.Name                    - %s\r\n', MovieParams.Name);
    fprintf(fid, 'MovieParams.Path                    - %s\r\n', MovieParams.Path);
    fprintf(fid, 'MovieParams.StartFrame              - %s\r\n', num2str(MovieParams.StartFrame));
    fprintf(fid, 'MovieParams.EndFrame                - %s\r\n', num2str(MovieParams.EndFrame));
    fprintf(fid, '\r\n');
    
    % PSF parameters
    fprintf (fid, 'PSF parameters\r\n');
    fprintf(fid, '************************************************\r\n');
    fprintf(fid, 'PSFParams.PSF_Name                  - %s\r\n', PSFParams.PSF_Name);
    fprintf(fid, 'PSFParams.PSF_Path                  - %s\r\n', PSFParams.PSF_Path);
    fprintf(fid, 'PSFParams.Xres                      - %s\r\n', num2str(PSFParams.Xres));
    fprintf(fid, 'PSFParams.Yres                      - %s\r\n', num2str(PSFParams.Yres));
    fprintf(fid, '\r\n');
    
    % Algorithm type parameters
    fprintf (fid, 'Algorithm parameters\r\n');
    fprintf(fid, '************************************************\r\n');
    fprintf(fid, 'AlgParams.AlgType                   - %s\r\n', AlgParams.AlgType);
    fprintf(fid, 'WienerParams.NSR                    - %s\r\n', num2str(WienerParams.NSR));
    fprintf(fid, 'AlgParams.CSFactor                  - %s\r\n', num2str(AlgParams.CSFactor));
    fprintf(fid, 'AlgParams.SegThreshold              - %s\r\n', num2str(AlgParams.SegThreshold));
    fprintf(fid, '\r\n');
    
    % Patches parameters
    fprintf (fid, 'Patches parameters\r\n');
    fprintf(fid, '************************************************\r\n');
    try
        DoPatchesFlag = GenParams.Patches.DoPatches;
    catch
        DoPatchesFlag = 0;
    end
    if DoPatchesFlag
        fprintf(fid, 'GenParams.Patches.DoPatches         - %s\r\n', num2str(GenParams.Patches.DoPatches));
        fprintf(fid, 'GenParams.Patches.PatchDimensions.X - %s\r\n', num2str(GenParams.Patches.PatchDimensions.X));
        fprintf(fid, 'GenParams.Patches.PatchDimensions.Y - %s\r\n', num2str(GenParams.Patches.PatchDimensions.Y));
        fprintf(fid, 'GenParams.Patches.Overlap.X         - %s\r\n', num2str(GenParams.Patches.Overlap.X));
        fprintf(fid, 'GenParams.Patches.Overlap.Y         - %s\r\n', num2str(GenParams.Patches.Overlap.Y));
    else
        fprintf(fid, 'GenParams.Patches.DoPatches         - %s\r\n', num2str(0));
        fprintf(fid, 'MovieParams.Xstart                  - %s\r\n', num2str(MovieParams.Xstart));
        fprintf(fid, 'MovieParams.Ystart                  - %s\r\n', num2str(MovieParams.Ystart));
        fprintf(fid, 'MovieParams.BlockSize               - %s\r\n', num2str(MovieParams.BlockSize));
    end
    fprintf(fid, 'MovieParams.SliceSize               - %s\r\n', num2str(MovieParams.SliceSize));
    fprintf(fid, 'GenParams.NoisePatchThr             - %s\r\n', num2str(GenParams.NoisePatchThr));
    fprintf(fid, '\r\n');
    
    % Optimization solver parameters
    fprintf (fid, 'Optimization solver parameters\r\n');
    fprintf(fid, '************************************************\r\n');
    fprintf(fid, 'AlgParams.Solver                    - %s\r\n', AlgParams.Solver);
    fprintf(fid, 'AlgParams.EigThreshold              - %s\r\n', num2str(AlgParams.EigThreshold));
    fprintf(fid, 'AlgParams.AvOverFrames              - %s\r\n', num2str(AlgParams.AvOverFrames));
    switch lower(AlgParams.Solver)
        case 'matfista'
            fprintf(fid, 'AlgParams.MatFISTA.IterMax         - %s\r\n', num2str(AlgParams.MatFISTA.IterMax));
            fprintf(fid, 'AlgParams.MatFISTA.L0              - %s\r\n', num2str(AlgParams.MatFISTA.L0));
            fprintf(fid, 'AlgParams.MatFISTA.Lambda          - %s\r\n', num2str(AlgParams.MatFISTA.Lambda));
            fprintf(fid, 'AlgParams.MatFISTA.NonNegOrth      - %s\r\n', num2str(AlgParams.MatFISTA.NonNegOrth));
            fprintf(fid, 'AlgParams.MatFISTA.Beta            - %s\r\n', num2str(AlgParams.MatFISTA.Beta));
            fprintf(fid, 'AlgParams.MatFISTA.LambdaBar       - %s\r\n', num2str(AlgParams.MatFISTA.LambdaBar));
        case 'matomp'
            fprintf(fid, 'AlgParams.MatOMP.IterMax            - %s\r\n', num2str(AlgParams.MatOMP.IterMax));
            fprintf(fid, 'AlgParams.MatOMP.Tol                - %s\r\n', num2str(AlgParams.MatOMP.Tol));
            fprintf(fid, 'AlgParams.MatOMP.NonNegOrth         - %s\r\n', num2str(AlgParams.MatOMP.NonNegOrth));
        case 'cvx'
            fprintf(fid, 'AlgParams.CVX_Lambda                - %s\r\n', num2str(AlgParams.CVX_Lambda));
        case 'fpg_mmv'
            fprintf(fid, 'AlgParams.FPG_MMV.IterMax           - %s\r\n', num2str(AlgParams.FPG_MMV.IterMax));
            fprintf(fid, 'AlgParams.FPG_MMV.L                 - %s\r\n', num2str(AlgParams.FPG_MMV.L));
            fprintf(fid, 'AlgParams.FPG_MMV.Lambda            - %s\r\n', num2str(AlgParams.FPG_MMV.Lambda));
        case 'omp_mmv'
            fprintf(fid, 'AlgParams.OMP_MMV.IterMax           - %s\r\n', num2str(AlgParams.OMP_MMV.IterMax));
            fprintf(fid, 'AlgParams.OMP_MMV.ResThreshold      - %s\r\n', num2str(AlgParams.OMP_MMV.ResThreshold));
            fprintf(fid, 'AlgParams.OMP_MMV.ResvsSolThreshold - %s\r\n', num2str(AlgParams.OMP_MMV.ResvsSolThreshold));
    end
    fprintf(fid, '\r\n');
    
    % Post processing parameters
    fprintf (fid, 'Post processing parameters\r\n');
    fprintf(fid, '************************************************\r\n');
    fprintf(fid, 'GenParams.PostProc.Method           - %s\r\n', GenParams.PostProc.Method);
    if strcmpi(GenParams.PostProc.Method, 'gaussian')
        try
            fprintf(fid, 'GenParams.PostProc.Sigma     - %s\r\n', num2str(GenParams.PostProc.Sigma));
        end
        try
            fprintf(fid, 'GenParams.PostProc.ColorMap  - %s\r\n', GenParams.PostProc.ColorMap);
        end
    end
    fprintf(fid, '\r\n');
    
    fprintf(fid, '* --------------------------------------------------------------------------------------- *\r\n');
    SucFlag = fclose(fid);
end

















