% DataFolder = '/media/mbakker/Microstroke-II/Marleen/GCaMP/';
% DataFolder = '/media/mbakker/PJM - HDD - 2/Marleen/GCaMP/';

function Pipeline_CaCl_GCaMP

%% Pick right recording
% CompareMovementRecordings('/media/mbakker/Microstroke-II/Marleen/GCaMP/');
% CompareMovementRecordings('/media/mbakker/PJM - HDD - 2/Marleen/GCaMP/');

load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat')
Recordings = [RecordingOverview.A1; RecordingOverview.A2; RecordingOverview.A3];

SaveDir = '/media/mbakker/SSD-2TB/P2/GCaMP';


%% Start going per recording
for ind = 1:size(Recordings,1)
    DataFolder = Recordings{ind};
    if( ~strcmp(DataFolder(end), filesep) )
        DataFolder = [DataFolder filesep];
    end
    
    %The place that you get the data from is different than the one you
    %save it in. 
    seps = strfind(DataFolder, filesep);    
    SaveFolder = [SaveDir DataFolder(seps(end-2):end)];
    if ~exist(SaveFolder, 'dir')
        mkdir([SaveFolder]);
    end
    
    anaReg = matfile([SaveFolder 'anaReg.mat'] ,'Writable',true);
    
    %% ImagesClassification
    disp(['ImagesClassification ', DataFolder]);
    varlist = who(anaReg,'ImagesClassification');
    if( isempty(varlist) || ~isfield(anaReg.ImagesClassification, 'ended') ) %|| is "or" maar kijkt eerst of de eerste klopt voordat het naar de tweede kijkt
        %   disp('Running...')
        anaReg.ImagesClassification = [];
        ImgClass.started = datestr(now);
        try
            ImagesClassification(DataFolder, SaveFolder, 1, 1, 1, 0, 'Internal-main');
            ImgClass.ended = datestr(now);
        catch
            ImgClass.error = datestr(now);
            anaReg.ImagesClassification = ImgClass;
            disp(['ImagescClassification error' DataFolder])
            return;
        end
        anaReg.ImagesClassification = ImgClass;
        clear ImgClass;
    else
        disp('already done.');
    end
    
    disp('ImageClassification Done')
    

    %% Correct For Dichroic Mirror shift
disp(['Dichroic correction ', SaveFolder]);
    varlist = who(anaReg,'CorrectDichroic');
    
       if( isempty(varlist) || ~isfield(anaReg.CorrectDichroic, 'ended') ) 
        anaReg.CorrectDichroic = [];
        CorrDi.started = datestr(now);
        try
            CorrectDichroic(SaveFolder)
            CorrDi.ended = datestr(now);
        catch
            CorrDi.error = datestr(now);
            anaReg.CorrectDichroic = CorrDi;
            disp(['CorrectDichroic error' DataFolder])
            return;
        end
        anaReg.CorrectDichroic = CorrDi;
        clear ImgClass;
       end
    
        %% Remove fluo 475 files
    % *** BE CAREFUL WITH THIS ***
    if contains(SaveFolder, 'SSD-2TB') %% second statement to make sure you dont delete from other places
        answer = questdlg(['Do you want to delete ' SaveFolder ' fluo files?'],...
            'WARNING - DELETING FILES',...
            'Yes, delete', 'No, cancel', 'No, cancel');
        switch answer
            case 'Yes, delete'
                
                eval(['delete ' SaveFolder 'fluo_475.dat'])
                eval(['delete ' SaveFolder 'fluo_475.mat'])
                disp('fluo 475 files deleted')
            case ' No, cancel'
                disp('files not deleted')
        end
    end
    
    %% step 3:
%     disp(' ')
    
%     %% Intra-Coregistration
%     disp('IntraCoregCalc');
%     varlist = who(anaReg,'IntraCoregCalc');
%     if( isempty(varlist) || ~isfield(anaReg.IntraCoregCalc, 'ended') )
%         disp('Running...')
%         anaReg.IntraCoregCalc = [];
%         IntraCoregCalc.started = datestr(now);
%         try
%             IntraCoRegistration(DataFolder);
%             IntraCoregCalc.ended = datestr(now);
%         catch
%             IntraCoregCalc.error = datestr(now);
%             anaReg.IntraCoregApp = IntraCoregCalc;
%             disp(['IntraCorecCalc error' DataFolder])
%             return;
%         end
%         anaReg.IntraCoregCalc = IntraCoregCalc;
%         clear IntraCoregCalc;
%     else
%         disp('already done.');
%     end
%     
%     disp('IntraCoregApp');
%     varlist = who(anaReg,'IntraCoregApp');
%     if( isempty(varlist) || ~isfield(anaReg.IntraCoregApp, 'ended') )
%         disp('Running...')
%         anaReg.IntraCoregApp = [];
%         IntraCoregApp.started = datestr(now);
%         try
%             ApplyIntraCoreg(DataFolder);
%             IntraCoregApp.ended = datestr(now);
%         catch
%             IntraCoregApp.error = datestr(now);
%             anaReg.IntraCoregApp = IntraCoregApp;
%             disp(['IntraCoregApplication error' DataFolder])
%             return;
%         end
%         anaReg.IntraCoregApp = IntraCoregApp;
%         clear IntraCoregApp;
%     else
%         disp('already done.');
%     end
%     
    
end
end