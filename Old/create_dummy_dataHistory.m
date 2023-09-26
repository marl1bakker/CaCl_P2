%% Create dummy dataHistory and add to processed files in SaveFolder

% Be sure that the .JSON pipeline file that you will use contains
% only the functions that were used to create the .dat file. So, one for
% the fluo and another for the HbO and HbR.

% Load pipeline file:
% myPipeline = jsondecode(fileread('F:\PATH_TO_PIPELINE_FILE\myPipeline.json'));
% myPipeline = jsondecode(fileread('/media/mbakker/GDrive/P2/GCaMP/PipeLineConfigFiles/GCaMP_umIT_Pipeline_2-5-23.json'));
myPipeline = jsondecode(fileread('/media/mbakker/GDrive/P2/GCaMP/PipeLineConfigFiles/Hemocompute.json'));
% Create dummy dataHistory
for ii = 1:length(myPipeline)
    fcnInfo = dir(which(myPipeline(ii).name));% Find the function's .m file info.
    [folder,name,~] = fileparts(fcnInfo.name); % Get name (without ext)
    % Create dummy data history
    dummy_dh(ii) = genDataHistory(struct('name',name,'folder',folder, 'datenum',fcnInfo.datenum),...
        'out = DummyFcnStr',myPipeline(ii).opts,{'none'},myPipeline(ii).inputFileName);
end
%% Add dummy dataHistory to all processed files with name "fileName".
targetFileName = 'fluo_567.mat';
% saveDir = 'F:\PATH_TO_PROJECT_SAVEFOLDER';
saveDir = '/media/mbakker/GDrive/P2/GCaMP';

% Insert dummy dataHistory in .mat files:
fileList = dir(fullfile(saveDir,'**',targetFileName));
for ii = 1:length(fileList)    
    md = matfile(fullfile(fileList(ii).folder, fileList(ii).name),'Writable',true);
    md.dataHistory = dummy_dh;
end
clear md




