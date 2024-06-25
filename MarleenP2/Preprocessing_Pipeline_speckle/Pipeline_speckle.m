% Pipeline_speckle
% make sure you did Preprocessing_Pipeline_speckle

%% Per mouse
load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview_Speckle.mat')
RecordingOverview = sortrows(RecordingOverview,9);
Recordings = [RecordingOverview.A1; RecordingOverview.A2; RecordingOverview.A3];
Mice = [RecordingOverview.Mouse; RecordingOverview.Mouse; RecordingOverview.Mouse];
SaveDirs = [RecordingOverview.SaveDirectory; RecordingOverview.SaveDirectory;RecordingOverview.SaveDirectory;];

for ind = 1:size(Recordings,1)
    Mouse = Mice{ind};
    DataFolder = Recordings{ind};
    if matches(DataFolder, 'empty')
        continue
    end
    
    fs = strfind(DataFolder, '-');
    Acq = DataFolder(fs(end-1)+1:fs(end)-1);
    clear fs
    
    if( ~strcmp(DataFolder(end), filesep) )
        DataFolder = [DataFolder filesep];
    end
    
    SaveFolder = [SaveDirs{ind} filesep Mouse filesep DataFolder(end-6:end)]; 

    anaReg = matfile([SaveFolder 'anaReg.mat'] ,'Writable',true);
    fprintf(['.......... \n \n' Mouse '\n']);
    disp(DataFolder(end-9:end-1))

    %% check if preprocessing pipeline is done
    Preproc_functions = {'ImagesClassification', 'MarkMovement', 'FlipLR', ...
        'BFI', 'ROI', 'Coreg_acquisitions'};

    for indfcndone = 1:size(Preproc_functions, 2)
        eval(['varlist = who(anaReg, ''' Preproc_functions{indfcndone} ''');'])
        if ( isempty(varlist) || ...
                eval(['~isfield(anaReg.' Preproc_functions{indfcndone} ', ''ended'')']) )
            disp('Preprocessing_Pipeline_CaCl_GCaMP not done yet. Do that first.')
            cont = 1;
            break 
        end
    end
    
    if exist('cont', 'var')
        continue
    end

    %% Table for Right vs left sides
    disp('Make ratio table...');
    varlist = who(anaReg,'Ratio');
    if( isempty(varlist) || ~isfield([anaReg.Ratio], 'ended') )
        anaReg.Ratio = [];
        ratiotable.started = datestr(now);
        
        try
            RightvsLeft_MakeTable(SaveFolder, 1)
            ratiotable.ended = datestr(now);
            anaReg.Ratio = ratiotable;
            
        catch
            ratiotable.error = datestr(now);
            anaReg.Ratio = ratiotable;
            
            disp(['Ratio table error ' DataFolder]);
            return;
        end
        clear ratiotable*
        disp('Ratio table done');
    else
        disp('Ratio table already done');
    end

end


%% Together
