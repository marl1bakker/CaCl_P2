%% Manually search for the best recordings:
% DataFolder = '/media/mbakker/PJM - HDD - 2/Marleen/GCaMP/';
% CompareMovementRecordings('/media/mbakker/PJM - HDD - 2/Marleen/GCaMP/');
% DataFolder = '/media/mbakker/Microstroke-II/Marleen/GCaMP/';
% CompareMovementRecordings('/media/mbakker/Microstroke-II/Marleen/GCaMP/');

%% Transfer recordings
% % Then make sure you have all the right recordings in the same place:
% load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat')
% Recordings = [RecordingOverview.A1; RecordingOverview.A2; RecordingOverview.A3];
% 
% for ind = 1:size(RecordingOverview, 1) % go per mouse
%     MouseFolders = dir([RecordingOverview.Folder{ind} ...
%         filesep RecordingOverview.Mouse{ind}]);
%     MouseFolders = MouseFolders(3:end); % get rid of ' and ''
%     folder = [MouseFolders.isdir];
%     MouseFolders = MouseFolders(folder == 1); %get rid of non-folder files
%     clear folder
%     
%     for idx = 1:size(MouseFolders,1)
%         CurrentFolder = [MouseFolders(idx).folder filesep MouseFolders(idx).name];
%         if sum(contains(Recordings, CurrentFolder)) == 0 %folder is not approved
%             DestinationFolder = [RecordingOverview.Folder{ind} '_Non_Used'];
%             movefile(CurrentFolder, DestinationFolder);
%         end
%     end
% end

%% Pipeline_CaCl_GCaMP
ManualInput = 0;

%% Pick right recording
load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat')
Recordings = [RecordingOverview.A1; RecordingOverview.A2; RecordingOverview.A3];
Mice = [RecordingOverview.Mouse; RecordingOverview.Mouse; RecordingOverview.Mouse];

SaveDir = '/media/mbakker/GDrive/P2/GCaMP';


%% Start going per recording
for ind = 1:size(Recordings,1)
    Mouse = Mice{ind};
    DataFolder = Recordings{ind};
    
    if( ~strcmp(DataFolder(end), filesep) )
        DataFolder = [DataFolder filesep];
    end
    
    %The place that you get the data from is different than the one you save it in.
%     %old
%     seps = strfind(DataFolder, filesep);
%     SaveFolder = [SaveDir DataFolder(seps(end-2):end)];
%     
    % Save to work with UmIToolbox
    SaveFolder = [SaveDir filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];    
    
    if ~exist(SaveFolder, 'dir')
        mkdir([SaveFolder]);
    end
    
    anaReg = matfile([SaveFolder 'anaReg.mat'] ,'Writable',true);
    disp(Mouse)
    disp(DataFolder(end-9:end-1))
    
    %% ImagesClassification
    disp(['ImagesClassification... ']);
    varlist = who(anaReg,'ImagesClassification');
    if( isempty(varlist) || ~isfield(anaReg.ImagesClassification, 'ended') ) %|| is "or" maar kijkt eerst of de eerste klopt voordat het naar de tweede kijkt
        %   disp('Running...')
        anaReg.ImagesClassification = [];
        ImgClass.started = datestr(now);
        try
            ImagesClassification(DataFolder, SaveFolder, 1, 1, 1, 0, 'Internal-main');
            ImgClass.ended = datestr(now);
            disp('ImageClassification Done')
        catch
            ImgClass.error = datestr(now);
            anaReg.ImagesClassification = ImgClass;
            disp(['ImagescClassification error' DataFolder])
            return;
        end
        anaReg.ImagesClassification = ImgClass;
        clear ImgClass;
    else
        disp('ImagesClassification already done');
    end
    
    %% Correct For Dichroic Mirror shift
    disp(['Dichroic correction...']); 
    varlist = who(anaReg,'CorrectDichroic');
    
    if( isempty(varlist) || ~isfield(anaReg.CorrectDichroic, 'ended') )...
            && ManualInput == 1
        anaReg.CorrectDichroic = [];
        CorrDi.started = datestr(now);
        try
            CorrectDichroic(SaveFolder)
            CorrDi.ended = datestr(now);
            disp('Dichroic Correction done');
        catch
            CorrDi.error = datestr(now);
            anaReg.CorrectDichroic = CorrDi;
            disp(['CorrectDichroic error' DataFolder]);
            return;
        end
        anaReg.CorrectDichroic = CorrDi;
        clear Corrdi;
    else
       disp('Dichroic Correction already done');
    end
    
    %% Flip data
    disp('Flip left and right...');
    
    for datfile = {'green', 'red', 'fluo_567'}
        varlist = who(anaReg,['FlipLR_' char(datfile)]);
        
        if( isempty(varlist) || ...
                eval(['~isfield([anaReg.FlipLR_' char(datfile) '], "ended")']) )
            eval(['anaReg.FlipLR_' char(datfile) '= [];']);
            eval(['flip' char(datfile) '.started = datestr(now);']);
            
            try
                fid = fopen([SaveFolder char(datfile) '.dat']);
                dat = fread(fid, inf, '*single');
                dat = reshape(dat, 512, 512, []);
%                 imagesc(dat(:,:,23))
                dat = fliplr(dat);
%                 figure()
%                 imagesc(dat(:,:,23))
                
                disp(['Overwriting ' char(datfile) '.dat file'])
                fclose(fid);
                fid = fopen([SaveFolder char(datfile) '.dat'], 'w');
                fwrite(fid,dat,'single');
                fclose(fid);
                eval(['flip' char(datfile) '.ended = datestr(now);']);
                eval(['anaReg.FlipLR_' char(datfile) '= flip' char(datfile) ';']);
                
            catch
                eval(['flip' char(datfile) '.error = datestr(now);']);
                eval(['anaReg.FlipLR_' char(datfile) '= flip' char(datfile)]);
                
                disp(['Flip left and right error ' DataFolder char(datfile)]);
                return;
            end
            
            clear flip*
            disp(['Flip left and right ' char(datfile) ' done']);
        else
            disp(['Flip left and right ' char(datfile) ' already done']);
        end
    end
    
    
    %% Remove fluo 475 files
    % *** BE CAREFUL WITH THIS ***
    if contains(SaveFolder, 'GDrive') && exist([SaveFolder 'fluo_475.dat']) ...  %to make sure you dont delete from other places
            && ManualInput == 1
        answer = questdlg(['Do you want to delete ' SaveFolder ' fluo 475 files?'],...
            'WARNING - DELETING FILES',...
            'Yes, delete', 'No, cancel', 'No, cancel');
        switch answer
            case 'Yes, delete'
                
                eval(['delete ' SaveFolder 'fluo_475.dat'])
                eval(['delete ' SaveFolder 'fluo_475.mat'])
                disp('Fluo 475 files deleted')
            case ' No, cancel'
                disp('Fluo 475 files not deleted')
        end
    end
    
    
end

% umIToolbox