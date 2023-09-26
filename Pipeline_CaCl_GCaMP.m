% Preprocessing_Pipeline_CaCl_GCaMP

% spo2 compute?
% gen. seeds/av of regions and timecourses

%% After umIT
load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat')
Recordings = [RecordingOverview.A1; RecordingOverview.A2; RecordingOverview.A3];
Mice = [RecordingOverview.Mouse; RecordingOverview.Mouse; RecordingOverview.Mouse];

SaveDir = '/media/mbakker/GDrive/P2/GCaMP';
Overwrite = 0; %don't overwrite stuff

%% single mouse
for ind = 1:size(Recordings,1)
    Mouse = Mice{ind};
    DataFolder = Recordings{ind};
    
    if( ~strcmp(DataFolder(end), filesep) )
        DataFolder = [DataFolder filesep];
    end

    SaveFolder = [SaveDir filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];      
%     anaReg = matfile([SaveFolder 'anaReg.mat'] ,'Writable',true);
    disp(Mouse)
    disp(DataFolder(end-9:end-1))
    
    %% Get Timecourses
    GetTimecourses(SaveFolder, 'average', Overwrite, 'hemoCorr_fluo') 
    GetTimecourses(SaveFolder, 'average', Overwrite, 'HbO') 

    %% Make corr Matrix single subject
%     SingleSubjectCorrMatrix(SaveFolder, Overwrite)

    %% Make SPCM single subject
%     SingleSubjectSPCM(SaveFolder, 'hemoCorr_fluo')
    
end

%% Mice combined

CombinedCorrMatrix('hemoCorr_fluo', 'A1')
CombinedCorrMatrix('hemoCorr_fluo', 'A2')
CombinedCorrMatrix('hemoCorr_fluo', 'A3')

CombinedCorrMatrix('HbO', 'A1')
CombinedCorrMatrix('HbO', 'A2')
CombinedCorrMatrix('HbO', 'A3')

CombinedSPCM('hemoCorr_fluo', 'A1')
CombinedSPCM('hemoCorr_fluo', 'A2')
CombinedSPCM('hemoCorr_fluo', 'A3')

CombinedSPCM('HbO', 'A1')
CombinedSPCM('HbO', 'A2')
CombinedSPCM('HbO', 'A3')


