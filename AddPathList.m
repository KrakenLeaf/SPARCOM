% This scipt adds the relevant folders for SPARCOM.
%
% Oren Solomon, Technion I.I.T. 
%

disp('Adding folders to the MATLAB path...');
disp('------------------------------------ ');
disp(' ');

% Add relative folders to the MATLAB path
disp('Adding relative folder: Utility');
addpath(genpath(fullfile(pwd, 'Utility')));
disp('Adding relative folder: Math');
addpath(genpath(fullfile(pwd, 'Math')));
disp('Adding relative folder: OptEngine');
addpath(genpath(fullfile(pwd, 'OptEngine')));
disp('Adding relative folder: CoreFunctions');
addpath(genpath(fullfile(pwd, 'CoreFunctions')));
disp(' ');

%% Third party toolboxes
% ------------------------------------------------------------------
ToolBoxDirectory = fullfile(pwd, 'External');
ToolBoxName = {
               'HNO',...                % HNO     
               'rwt-master',...         % RWT - Rice Wavelet Toolbox V 3.0
               'Bioformats 5.2.0',...   % Bioformats 5.2.0
               'FISTA_TV',...           % Total variation
              };

for ii = 1:length(ToolBoxName)
    PackFolder = fullfile(ToolBoxDirectory, ToolBoxName{ii});
    disp([ToolBoxName{ii} ' : adding folder ' PackFolder]);
    addpath(genpath(PackFolder));
end

%% ------------------------------------------------------------------------
% Save path permanently
disp(' ');
disp('Saving MATLAB path permanently...');
savepath;

disp(' ');
disp('Done.');
