%% preprocessing GCaMP pipeline CaCl project
% Start here for the pipeline. After this step, do
% Umit_Pipeline_CaCl_GCaMP. After that, do Pipeline_CaCl_GCaMP. 

%% Manually search for the best recordings, make RecordingOverview:
% CompareMovementRecordings('/media/mbakker/PJM - HDD - 2/Marleen/GCaMP/');
% CompareMovementRecordings('/media/mbakker/Microstroke-II/Marleen/GCaMP/');
% CompareMovementRecordings('/media/mbakker/SSD-2TB/GCaMP/');

%% Transfer recordings
% % Then make sure you have all the right recordings in the same place:
% % This is commented because there's no need to keep doing it, but it
% % would not give a problem if it was uncommented either.
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
SaveDirs = [RecordingOverview.SaveDirectory; RecordingOverview.SaveDirectory;RecordingOverview.SaveDirectory;];

% Recordings = [RecordingOverview.A3];
% Mice = [RecordingOverview.Mouse];
% SaveDirs = [RecordingOverview.SaveDirectory];

% Recordings = [RecordingOverview.A2; RecordingOverview.A3];
% Mice = [RecordingOverview.Mouse; RecordingOverview.Mouse];
% SaveDirs = [RecordingOverview.SaveDirectory; RecordingOverview.SaveDirectory];

% SaveDir = '/media/mbakker/GDrive/P2/GCaMP';

%% left-right
% DataFolder = '/media/mbakker/Microstroke-II/Marleen/LR-GCAMP/Rawdat';
% SaveFolder = '/media/mbakker/Microstroke-II/Marleen/LR-GCAMP/M99/A1-R1/CtxImg';

%% Start going per recording
for ind = 1:size(Recordings,1)
    Mouse = Mice{ind};
    DataFolder = Recordings{ind};
    fs = strfind(DataFolder, '-');
    Acq = DataFolder(fs(end-1)+1:fs(end)-1);
    clear fs

    if( ~strcmp(DataFolder(end), filesep) )
        DataFolder = [DataFolder filesep];
    end
    
    % Save to work with UmIToolbox
    SaveFolder = [SaveDirs{ind} filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];    
    
    if ~exist(SaveFolder, 'dir')
        mkdir(SaveFolder);
    end
    
    anaReg = matfile([SaveFolder 'anaReg.mat'] ,'Writable',true);
    anaReg.tostart = 1;
    fprintf(['.......... \n \n' Mouse '\n']);
    disp(DataFolder(end-9:end-1))
    
    %% check if this pipeline is already done (saves printing space)
    Preproc_functions = {'ImagesClassification', 'MarkMovement', 'CorrectDichroic',...
        'FlipLR_fluo_567', 'FlipLR_green', 'FlipLR_red', 'Coreg_acquisitions'};
    teller = 0;

    for indfcndone = 1:size(Preproc_functions, 2)
        % eval(['varlist = who(anaReg, ''' Preproc_functions{indfcndone} ''');'])
        varlist = who(anaReg, Preproc_functions{indfcndone});
        if ( isempty(varlist) || ...
                ~isfield(anaReg.(Preproc_functions{indfcndone}), 'ended') )

            teller = teller +1;
            break
        end
    end

    if teller == 0
        disp('Preprocessing_Pipeline_CaCl_GCaMP already done. Mouse skipped.')
        continue % if all functions here are done, go to next acquisition.
    end
    
    %% ImagesClassification
    disp('ImagesClassification... ');
    varlist = who(anaReg,'ImagesClassification');
    if( isempty(varlist) || ~isfield(anaReg.ImagesClassification, 'ended') ) %|| is "or" maar kijkt eerst of de eerste klopt voordat het naar de tweede kijkt
        anaReg.ImagesClassification = [];
        ImgClass.started = datestr(now);
        try
            ImagesClassification(DataFolder, SaveFolder, 1, 1, 0);
            ImgClass.ended = datestr(now);
            disp('ImageClassification Done')
        catch
            ImgClass.error = datestr(now);
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
        end
        anaReg.MarkMovement = MarkMov;
        clear MarkMov;
    else
        disp('Movement Marking already done');
    end

    %% Correct For Dichroic Mirror shift
    disp('Dichroic Correction...'); 
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
            disp(['CorrectDichroic error' DataFolder]);
            return;
        end
        anaReg.CorrectDichroic = CorrDi;
        clear CorrDi;
    elseif datetime(getfield(anaReg.CorrectDichroic, 'ended')) < datetime(getfield(anaReg.ImagesClassification, 'ended')) ...
            && ManualInput == 1
        disp('Redid Imagesclassification so have to redo CorrectDichroic...')
                
        anaReg.CorrectDichroic = [];
        CorrDi.started = datestr(now);
        try
            CorrectDichroic(SaveFolder)
            CorrDi.ended = datestr(now);
            disp('Dichroic Correction done');
        catch
            CorrDi.error = datestr(now);
            disp(['CorrectDichroic error' DataFolder]);
            return;
        end
        anaReg.CorrectDichroic = CorrDi;
        clear CorrDi;
    else
       disp('Dichroic Correction already done');
    end

    %% Flip data
    disp('Flip left and right...');

    for datfile = {'green', 'red', 'fluo_567'}
        varlist = who(anaReg,['FlipLR_' char(datfile)]);

        if( isempty(varlist) || ...
                eval(['~isfield([anaReg.FlipLR_' char(datfile) '], "ended")']) ) || ...
                eval(['datetime(getfield(anaReg.FlipLR_' char(datfile) ', ''ended'')) < datetime(getfield(anaReg.ImagesClassification, ''ended''))'])

            eval(['anaReg.FlipLR_' char(datfile) '= [];']);
            eval(['flip' char(datfile) '.started = datestr(now);']);

            try
                fid = fopen([SaveFolder char(datfile) '.dat']);
                dat = fread(fid, inf, '*single');
                dat = reshape(dat, 512, 512, []);
                dat = fliplr(dat);

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
    % temp after correction of coregistration dichroic. Check if atlasmask
    % fits on brain. 

    % klad

    disp('Mask and ROI creation...')

    if ManualInput == 1
        if exist([SaveDirs{ind} filesep Mouse filesep 'ImagingReferenceFrame.mat'], 'file') && ...
                exist([SaveDirs{ind} filesep Mouse filesep 'ROImasks_data.mat'], 'file') && ...
                exist([SaveDirs{ind} filesep Mouse filesep 'BigROI.mat'], 'file')
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

            green_im = green_im./mean(green_im(:));
            green_im = (green_im-min(green_im(:)))./(max(green_im(:))-min(green_im(:)));
            green_im(green_im < 0) = 0;
            green_im(green_im > 1) = 1;
            green_im = adapthisteq(green_im);
            figure
            imagesc(green_im)

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

            % SeedGenerator only necessary if you do centroids, not whole
            % regions.
%             SeedGenerator(SaveFolder) %new 15-6-23, have to test
            ClusterRois(SaveFolder, 1)
            Correct_for_rotation_ROI(SaveFolder)

        end
    elseif ~exist([SaveDirs{ind} filesep Mouse filesep 'ImagingReferenceFrame.mat'], 'file')
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

    %% Coregistration A1 - A2 - A3    % Save to work with UmIToolbox
%     SaveFolder = [SaveDirs{ind} filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];    

    %It should work, finished it 23-8-23. Haven't tested
    %it. Already did the green.dat file for M13, A2-R1. None of the other
    %files though. Is the correction that we do for A1 saved in the
    %fluo_567 AND the green and red files during the ROI making? If yes,
    %there is no problem. If no, we have to make sure that all the images
    %are upright first. 
    %Have to recalculate the HbO and HbR files for A2 and A3, as well as
    %any fluo files. 
    disp('Coregistration A1-A2-A3...')
    varlist = who(anaReg,'Coreg_acquisitions');

    if isequal(Acq, 'A1') 
        disp('No coregistration, Acquisition is A1')
        anaReg.Coreg_acquisitions = [];
        Coreg.started = datestr(now);
        Coreg.ended = datestr(now);
        anaReg.Coreg_acquisitions = Coreg;
        clear Coreg

    elseif ManualInput == 0
        disp('Coregistration skipped')

    elseif ( isempty(varlist) || ~isfield(anaReg.Coreg_acquisitions, 'ended') || ...
            datetime(getfield(anaReg.Coreg_acquisitions, 'ended')) < datetime(getfield(anaReg.ImagesClassification, 'ended')) )
        anaReg.Coreg_acquisitions = [];
        Coreg.started = datestr(now);

        eval(['allinfomouse = RecordingOverview(ismember(RecordingOverview.' Acq ', ''' DataFolder(1:end-1) '''),:);']);
        pathFixed = allinfomouse.A1{:}; %get folder of A1 for this mouse
        pathFixed = [SaveDirs{ind} filesep Mouse filesep pathFixed(end-4:end) filesep 'CtxImg' filesep];

        try
            Coregistration(pathFixed, SaveFolder);
            Coreg.ended = datestr(now);
            disp('Coregistration done.')
        catch

            % disp('Coregistration try 2')
            % try
            %     CoregistrationManual(pathFixed, SaveFolder)
            %     Coreg.ended = datestr(now);
            %     disp('Coregistration done.')
            % catch
                Coreg.error = datestr(now);
                disp(['Coregistration error ' DataFolder])
                %             return;
            % end
        end
        anaReg.Coreg_acquisitions = Coreg;
        clear Coreg

    else
        disp('Coregistration already done');

    end
    
    
end


