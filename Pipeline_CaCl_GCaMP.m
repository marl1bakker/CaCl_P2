%% Manually search for the best recordings:
% CompareMovementRecordings('/media/mbakker/PJM - HDD - 2/Marleen/GCaMP/');
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
ManualInput = 1;

%% Pick right recording
load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat')
Recordings = [RecordingOverview.A1; RecordingOverview.A2; RecordingOverview.A3];
Mice = [RecordingOverview.Mouse; RecordingOverview.Mouse; RecordingOverview.Mouse];

SaveDir = '/media/mbakker/GDrive/P2/GCaMP';

load('/media/mbakker/GDrive/P2/GCaMP/LogBook.mat') % to check later if it's already done by toolbox

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
    disp('ImagesClassification... ');
    varlist = who(anaReg,'ImagesClassification');
    if( isempty(varlist) || ~isfield(anaReg.ImagesClassification, 'ended') ) %|| is "or" maar kijkt eerst of de eerste klopt voordat het naar de tweede kijkt
%           disp('Running...')
        anaReg.ImagesClassification = [];
        ImgClass.started = datestr(now);
        try
            ImagesClassification(DataFolder, SaveFolder, 1, 1, 1, 0, 'Internal-main');
            ImgClass.ended = datestr(now);
            disp('ImageClassification Done')
        catch
            ImgClass.error = datestr(now);
%             anaReg.ImagesClassification = ImgClass;
            disp(['ImagescClassification error ' DataFolder])
            return;
        end
        anaReg.ImagesClassification = ImgClass;
        clear ImgClass;
    else
        disp('ImagesClassification already done');
    end
    
    %% MarkMovedFrames
    disp('Mark frames with movement...');
    varlist = who(anaReg,'MarkMovement');
    if( isempty(varlist) || ~isfield(anaReg.MarkMovement, 'ended') ) %|| is "or" maar kijkt eerst of de eerste klopt voordat het naar de tweede kijkt
        anaReg.MarkMovement = [];
        MarkMov.started = datestr(now);
        try
            MarkMovedFrames(DataFolder, SaveFolder);
            MarkMov.ended = datestr(now);
            disp('Movement Marking Done')
        catch
            MarkMov.error = datestr(now);
            disp(['Movement Marking error ' DataFolder])
%             return;
        end
        anaReg.MarkMovement = MarkMov;
        clear MarkMov;
    else
        disp('Movement Marking already done');
    end
    
    %% Correct For Dichroic Mirror shift
    disp('Dichroic correction...'); 
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
%             anaReg.CorrectDichroic = CorrDi;
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
    
        
    %% ROI creating + mask
    if ManualInput == 1
        if exist([SaveDir filesep Mouse filesep 'ImagingReferenceFrame.mat'], 'file') && ...
                exist([SaveDir filesep Mouse filesep 'ROImasks_data.mat'], 'file')
            disp('Mask and ROI already created')
        else
            fid = fopen([SaveFolder 'fluo_567.dat']);
            fluo_im = fread(fid, 512*512, '*single' );
            fluo_im = reshape(fluo_im, 512, 512);
            disp(Mouse)
            fclose(fid);
            
            fid = fopen([SaveFolder 'green.dat']);
            green_im = fread(fid, 512*512, '*single' );
            green_im = reshape(green_im, 512, 512);
            imagesc(green_im)
            fclose(fid); 
            
            ROImanager(fluo_im)
            f = msgbox(["In ROImanager, do the following steps:"; ...
                "-  Image – Set origin – New – drag to bregma, right mouse-click to set"; ...
                "-  Image – Set origin – Align image to origin – drag to lambda, right mouse-click to set and correct for tilt in frame"; ...
                "-  Image – Set pixel size – 50 pixels per mm"; ...
                "-  Image – Mask – Draw new – Draw the mask by clicking points, you can still drag the points after setting them, double click inside the mask to confirm"; ...
                "-  Image – Image reference file… - Export – save as ImagingReferenceFrame.mat in the mouse folder"; ...
                "-  Create – Mouse Allen Brain Atlas – Areas – Select areas – close window to select them all – drag atlas to what seems to be fitting – double click inside atlas to confirm"; ...
                "-  File – Save as… - ROImasks_data.mat in mouse folder"]);            
            input('Make reference and ROI')

            SeedGenerator(SaveFolder) %new 15-6-23, have to test
        end
    elseif ~exist([SaveDir filesep Mouse filesep 'ImagingReferenceFrame.mat'], 'file')
        disp(['To do: ROI & Mask for ' Mouse ])
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
    
    %% umIT
% umIToolbox

% hemocorrection
% hemocompute
% align frames
% normalize LPF
    

% check what has been done by toolbox already
LogBookMouse = LogBook(find(matches(LogBook.Subject, Mouse)),:); %get everything for mouse

%get acquisition
Acquisition = DataFolder(end-5:end-1);
LogBookMouse = LogBookMouse(find(matches(LogBookMouse.Acquisition, Acquisition)),:);
    
    %% HemoCorrection
    
    disp('Hemodynamic correction...'); 
    varlist = who(anaReg,'HemoCorrection');
    
    if( isempty(varlist) || ~isfield(anaReg.HemoCorrection, 'ended') )... %if not done in pipeline before    
            && ManualInput == 1
        
        
            % and if not done in toolbox before

        
        anaReg.HemoCorrection = [];
        hemcorr.started = datestr(now);
        try
            [outData, metaData] = run_HemoCorrection(SaveFolder, varargin);
            
            hemcorr.ended = datestr(now);
            disp('Hemodynamic Correction done');
        catch
            hemcorr.error = datestr(now);
            anaReg.HemoCorrection = hemcorr;
            disp(['Hemodynamic Correction error ' DataFolder]);
            return;
        end
        anaReg.HemoCorrection = hemcorr;
        clear hemcorr;
    else
       disp('Hemodynamic Correction already done');
    end
    
    
    
end





%% After umIT
load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat')
% RecordingOverview = RecordingOverview(1:4,:); % subset M13-M16 A1
% Recordings = RecordingOverview.A1;
% Mice = RecordingOverview.Mouse;
Recordings = [RecordingOverview.A1; RecordingOverview.A2; RecordingOverview.A3];
Mice = [RecordingOverview.Mouse; RecordingOverview.Mouse; RecordingOverview.Mouse];

SaveDir = '/media/mbakker/GDrive/P2/GCaMP';

for ind = 1:size(Recordings,1)
    Mouse = Mice{ind};
    DataFolder = Recordings{ind};
    
    if( ~strcmp(DataFolder(end), filesep) )
        DataFolder = [DataFolder filesep];
    end

    SaveFolder = [SaveDir filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];      
    anaReg = matfile([SaveFolder 'anaReg.mat'] ,'Writable',true);
    disp(Mouse)
    disp(DataFolder(end-9:end-1))
    
        %% Get Timecourses
    disp('Get Timecourses for ROI, fluo & HbO...');
    varlist = who(anaReg,'Timecourses');
    if( isempty(varlist) || ~isfield(anaReg.Timecourses, 'ended') ) %|| is "or" maar kijkt eerst of de eerste klopt voordat het naar de tweede kijkt
        anaReg.Timecourses = [];
        Timcor.started = datestr(now);
        try
            GetTimecourses(SaveFolder) 
            Timcor.ended = datestr(now);
            disp('Timecourse calculation Done')
        catch
            Timcor.error = datestr(now);
            disp(['Timecourse calculation error ' SaveFolder])
%             return;
        end
        anaReg.Timecourses = Timcor;
        clear Timcor;
    else
        disp('Timecourse Calculation already done');
    end

    %% Make corr Matrix single subject
    SingleSubjectCorrMatrix(SaveFolder)
    
end
