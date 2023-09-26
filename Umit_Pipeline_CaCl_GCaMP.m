%% Umit_Pipeline_CaCl_GCaMP
% This part was first done by the Umit toolbox. It is now done function by
% function in this pipeline.

% umIToolbox

% hemocorrection
% align frames (coreg)
% normalize LPF
% hemocompute
% align frames hbo
% align frames hbr


% spo2 compute?
% cluster ROI
% gen. seeds and timecourses

% % check what has been done by toolbox already
% LogBookMouse = LogBook(find(matches(LogBook.Subject, Mouse)),:); %get everything for mouse


%% Pipeline_CaCl_GCaMP
ManualInput = 1;

%% Pick right recording
load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat')
Recordings = [RecordingOverview.A1; RecordingOverview.A2; RecordingOverview.A3];
% Recordings = RecordingOverview{1,4};
Mice = [RecordingOverview.Mouse; RecordingOverview.Mouse; RecordingOverview.Mouse];
% Mice = {'M13'};

SaveDir = '/media/mbakker/GDrive/P2/GCaMP';
% load('/media/mbakker/GDrive/P2/GCaMP/LogBook.mat') % to check later if it's already done by toolbox

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
    SaveFolder = [SaveDir filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];
    
    if ~exist(SaveFolder, 'dir')
        mkdir(SaveFolder);
    end
    
    disp(Mouse)
    disp(DataFolder(end-9:end-1))
    
    %% check if previous pipeline (Preprocessing_Pipeline_CaCl_GCaMP) is completely done
    anaReg = matfile([SaveFolder 'anaReg.mat'] ,'Writable',true);
    Preproc_functions = {'ImagesClassification', 'MarkMovement', 'CorrectDichroic',...
        'FlipLR_fluo_567', 'FlipLR_green', 'FlipLR_red', 'Coreg_acquisitions'};
    teller = 0;
    
    for indfcndone = 1:size(Preproc_functions, 2)
        if eval(['~isfield(anaReg.' Preproc_functions{indfcndone} ', ''ended'')'])
            disp(['Exited pipeline because ' Preproc_functions{indfcndone} ...
                ' is not done for ' Mouse])
            teller = teller +1;
            break
        end
    end
    
    if teller > 0
        continue % if a function of the preprocessing was not done, go to next acquisition.
    end
    
    %% HemoCorrection
    disp('Hemodynamic correction...');
    varlist = who(anaReg,'HemoCorrection');
    
    if( isempty(varlist) || ~isfield(anaReg.HemoCorrection, 'ended') )... %if not done in pipeline before
            && ManualInput == 1
        
        anaReg.HemoCorrection = [];
        hemcorr.started = datestr(now);
        try
            % need fluo.dat for this
            fid = fopen([SaveFolder 'fluo_567.dat']);
            dat = fread(fid, inf, '*single');
            metadata = load([SaveFolder 'fluo_567.mat']);
            hemoCorr_fluo = HemoCorrection(SaveFolder, dat, metadata, {'Red','Green'});
            
            %change nans to 0, to be able to do the normalization later
            hemoCorr_fluo(isnan(hemoCorr_fluo)) = 0;
            
            fid = fopen([SaveFolder 'hemoCorr_fluo.dat'],'w');
            fwrite(fid,hemoCorr_fluo,'*single');
            fclose(fid);
            
            Make_HemoCorr_Matfile(SaveFolder);
            
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
    
    %% Normalize and filter
    %     continue
    disp('Normalization...');
    
    varlist = who(anaReg,'Normalization');
    if( isempty(varlist) || (~isfield(anaReg.Normalization, 'ended')) )
        anaReg.Normalization = [];
        Normalization.started = datestr(now);
        try
            
%             %only temporary for the ones where you didnt change to zero at hemcorr
%             %yet
%             fid = fopen([SaveFolder 'hemoCorr_fluo.dat']);
%             hemoCorr_fluo = fread(fid, inf, '*single');
%             hemoCorr_fluo(isnan(hemoCorr_fluo)) = 0;
%             fclose(fid);
%             fid = fopen([SaveFolder 'hemoCorr_fluo.dat'], 'w');
%             fwrite(fid,hemoCorr_fluo,'*single');
%             fclose(fid);
            
%             %only needed for now, can delete when old ones are fixed.
%             if ~exist([SaveFolder 'hemoCorr_fluo.mat'], 'file')
%                 Make_HemoCorr_Matfile(SaveFolder);
%             end
            
            NormalisationFiltering(SaveFolder, 'hemoCorr_fluo', 0.3, 3, 1, 0);
            
            Normalization.ended = datestr(now);
        catch e
            Normalization.error = datestr(now);
            Normalization.def = e;
            anaReg.Normalization = Normalization;
            disp(['normalization error' DataFolder])
            continue;
        end
        anaReg.Normalization = Normalization;
        clear Normalization fid
    else
        disp('Normalization already done')
    end
    
%     continue
    %% HbO HbR
    disp('HbO/HbR calculation...')
    varlist = who(anaReg,'HbOHbR');
    
    if( isempty(varlist) || (~isfield(anaReg.HbOHbR, 'ended')) )
        anaReg.HbOHbR = [];
        HbOHbRvar.started = datestr(now);
        try
            %         HbOHbRCalculation(SaveFolder, 1);  %0 is negeer niet als files al bestaan, 1 is negeer het wel en overwrite
            [~, ~] = HemoCompute(SaveFolder, SaveFolder, 'gcamp', {'red', 'green'}, 1);
            
            HbOHbRvar.ended = datestr(now);
        catch e
            disp(['HbOHbR error' DataFolder])
            HbOHbRvar.error = datestr(now);
            HbOHbRvar.def = e;
            anaReg.HbOHbR = HbOHbRvar;
            throw(e)
        end
        anaReg.HbOHbR = HbOHbRvar;
        clear HbOHbRvar
    else
        disp('HbO/HbR calculation already done')
    end
    
    %% SpO2 calculation
    % disp('spO2 calculation')
    % varlist = who(anaReg,'spO2');
    %
    % if( isempty(varlist) || (~isfield(anaReg.spO2, 'ended')) )
    %     anaReg.spO2 = [];
    %     spO2var.started = datestr(now);
    %     try
    %         spO2Calculation(DataFolder, 1);  %0 is negeer niet als files al bestaan, 1 is negeer het wel en overwrite
    %         spO2var.ended = datestr(now);
    %     catch e
    %         disp(['spO2 error' DataFolder])
    %         spO2var.error = datestr(now);
    %         spO2var.def = e;
    %         anaReg.spO2 = spO2var;
    %         throw(e)
    %     end
    %     anaReg.spO2 = spO2var;
    %     clear spO2var
    % else
    %     disp('already done.')
    % end
    
    
    
end
