%Preprocessing_Pipeline_speckle

%% Manually search for the best recordings, make RecordingOverview:
% CompareMovementRecordings('/media/mbakker/PJM - HDD - 2/Marleen/Speckle/', 'Speckle');
% CompareMovementRecordings('/media/mbakker/Microstroke-II/Marleen/Speckle/', 'Speckle');

%% Transfer recordings
% % Then make sure you have all the right recordings in the same place:
% % This is commented because there's no need to keep doing it, but it
% % would not give a problem if it was uncommented either.
% load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview_Speckle.mat')
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
%             if ~isfolder([RecordingOverview.Folder{ind} '_Non_Used'])
%                 mkdir([RecordingOverview.Folder{ind} '_Non_Used'])
%             end
%             DestinationFolder = [RecordingOverview.Folder{ind} '_Non_Used'];
%             movefile(CurrentFolder, DestinationFolder);
%         end
%     end
% end


%% Pipeline_CaCl_GCaMP
ManualInput = 1;

%% Pick right recording
% To check if you need to flip left and right:
% load('/media/mbakker/GDrive/P2/TEST_Speckle/RecordingOverview_Speckle_LRtest.mat')
% Need to flip.

load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview_Speckle.mat')
Recordings = [RecordingOverview.A1; RecordingOverview.A2; RecordingOverview.A3];
Mice = [RecordingOverview.Mouse; RecordingOverview.Mouse; RecordingOverview.Mouse];
SaveDirs = [RecordingOverview.SaveDirectory; RecordingOverview.SaveDirectory;RecordingOverview.SaveDirectory;];

%% Start going per recording
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
    
    % Save to work with UmIToolbox
%     SaveFolder = [SaveDirs{ind} filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];    
    SaveFolder = [SaveDirs{ind} filesep Mouse filesep DataFolder(end-6:end)];    
    
    if ~exist(SaveFolder, 'dir')
        mkdir(SaveFolder);
    end
    
    anaReg = matfile([SaveFolder 'anaReg.mat'] ,'Writable',true);
    anaReg.tostart = 1;
    fprintf(['.......... \n \n' Mouse '\n']);
    disp(DataFolder(end-9:end-1))
    
    %% check if this pipeline is already done (saves printing space)
    Preproc_functions = {'ImagesClassification', 'MarkMovement', 'FlipLR'};
    teller = 0;
    
    for indfcndone = 1:size(Preproc_functions, 2)
        eval(['varlist = who(anaReg, ''' Preproc_functions{indfcndone} ''');'])
        if ( isempty(varlist) || ...
                eval(['~isfield(anaReg.' Preproc_functions{indfcndone} ', ''ended'')']) )
                
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
    
    %% Flip data
    disp('Flip left and right...');
    varlist = who(anaReg,'FlipLR');
    if( isempty(varlist) || ~isfield([anaReg.FlipLR], 'ended') )
        anaReg.FlipLR = [];
        flip.started = datestr(now);
        
        try
            fid = fopen([SaveFolder 'speckle.dat']);
            dat = fread(fid, inf, '*single');
            dat = reshape(dat, 1024, 1024, []);
            dat = fliplr(dat);
            
            fclose(fid);
            fid = fopen([SaveFolder 'speckle.dat'], 'w');
            disp('Overwriting speckle.dat file')
            fwrite(fid,dat,'single');
            flip.ended = datestr(now);
            anaReg.FlipLR = flip;
            fclose(fid);
            
        catch
            flip.error = datestr(now);
            anaReg.FlipLR = flip;
            
            disp(['Flip left and right error ' DataFolder]);
            return;
        end
        
        clear flip*
        disp('Flip left and right done');
    else
        disp('Flip left and right already done');
    end
    
    
    
    
end




fid = fopen('speckle.dat');
dat = fread(fid, inf, '*single');
dat = reshape(dat, 1024, 1024, []);
fclose(fid)

sdat = movstd(dat, 10, 0, 3);
mdat = movmean(dat, 10,3);
K = sdat./mdat;
clear sdat mdat
BFI = 1./(K.^2);
clear K

for ind = 1:size(BFI,3)
    BFI(:,:,ind) = medfilt2(BFI(:,:,ind), [5 5]);
end

medianBFI = zeros(1,size(BFI,3));
for ind = 1:size(BFI,3)
    medianBFI(ind) = median(BFI(500:850,200:750,ind), 'all');
end

for ind = 4000:4100
    imagesc(BFI(:,:,ind),[1 2000])
    title(num2str(ind))
    pause(0.05)
end






