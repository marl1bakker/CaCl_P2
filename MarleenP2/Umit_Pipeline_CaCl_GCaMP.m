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
Overwrite = 0;

%% Pick right recording
load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat')
Recordings = [RecordingOverview.A1; RecordingOverview.A2; RecordingOverview.A3];
Mice = [RecordingOverview.Mouse; RecordingOverview.Mouse; RecordingOverview.Mouse];
SaveDirs = [RecordingOverview.SaveDirectory; RecordingOverview.SaveDirectory;RecordingOverview.SaveDirectory;];

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
    
%     disp(Mouse)
    fprintf(['.......... \n \n' Mouse '\n']);
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
    
    clear Preproc_functions teller infcndone 
    
    %% check if this pipeline is already done (saves printing space)
	Umit_functions = {'HemoCorrection', 'Normalization', 'HbOHbR', ...
        'TimecourseFluo', 'TimecourseHbOHbR', 'AlignAllBrains'};
    teller = 0;
    
    for indfcndone = 1:size(Umit_functions, 2)
        eval(['varlist = who(anaReg, ''' Umit_functions{indfcndone} ''');'])
        if ( isempty(varlist) || ...
                eval(['~isfield(anaReg.' Umit_functions{indfcndone} ', ''ended'')']) )
            teller = teller +1;
            break
        end
    end
    
    if teller == 0
        disp('Umit_Pipeline_CaCl_GCaMP already done. Mouse skipped.')
        continue % if all functions here are done, go to next acquisition.
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

            hemoCorr_fluo = NormalisationFiltering(SaveFolder, 'hemoCorr_fluo', 0.3, 3, 1, 0);
            
            fid = fopen([SaveFolder 'hemoCorr_fluo.dat'],'w');
            fwrite(fid,hemoCorr_fluo,'*single');
            fclose(fid);   
            
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
    

    %% Timecourses Fluo
    disp('Timecourse calculation Fluo...')
    varlist = who(anaReg,'TimecourseFluo');
    
    if( isempty(varlist) || (~isfield(anaReg.TimecourseFluo, 'ended')) )
        anaReg.TimecourseFluo = [];
        Timecoursesvar.started = datestr(now);
        
        % check if you have all the necessary things
        if isfield(anaReg.HemoCorrection, 'ended') && isfield(anaReg.Normalization, 'ended')
 
            try %if yes, try to calculate timecourses
            GetTimecourses(SaveFolder, 'average', Overwrite, 'hemoCorr_fluo')
            Timecoursesvar.ended = datestr(now);
            
            catch e
            disp(['Timcoureses Fluo error' DataFolder])
            Timecoursesvar.error = datestr(now);
            Timecoursesvar.def = e;
%             anaReg.TimecoursesFluo = Timecoursesvar;
            throw(e)
            end
        
        else % if functions before are not done
            disp('Timecourses fluo not calculated, HemoCorr or Norm not done')
            Timecoursesvar.error = datestr(now);
            
        end
        
        anaReg.TimecourseFluo = Timecoursesvar;
        clear Timecoursesvar
    else
        disp('Timecourses Fluo calculation already done')
    end
    
    %% Timecourses HbO/HbR
    disp('Timecourse calculation HbO/HbR...')
    varlist = who(anaReg,'TimecourseHbOHbR');
    
    if( isempty(varlist) || (~isfield(anaReg.TimecourseHbOHbR, 'ended')) )
        anaReg.TimecourseHbOHbR = [];
        Timecoursesvar.started = datestr(now);
        
        % check if you have all the necessary things
        if isfield(anaReg.HbOHbR, 'ended')
 
            try %if yes, try to calculate timecourses
            GetTimecourses(SaveFolder, 'average', Overwrite, 'HbO')
            GetTimecourses(SaveFolder, 'average', Overwrite, 'HbR')
            Timecoursesvar.ended = datestr(now);
            
            catch e
            disp(['Timcoureses Fluo error' DataFolder])
            Timecoursesvar.error = datestr(now);
            Timecoursesvar.def = e;
%             anaReg.TimecoursesHbOHbR = Timecoursesvar;
            throw(e)
            end
        
        else % if functions before are not done
            disp('Timecourses HbO/HbR not calculated, .dat files dont exist')
            Timecoursesvar.error = datestr(now);
            
        end
        
        anaReg.TimecourseHbOHbR = Timecoursesvar;
        clear Timecoursesvar
    else
        disp('Timecourses HbOHbR calculation already done')
    end        
            

%     %% calculate the tform to align the mice/brains
%     disp('Align to all brains...')
%     varlist = who(anaReg,'AlignAllBrains');
%     
%     if( isempty(varlist) || (~isfield(anaReg.AlignAllBrains, 'ended')) )
%         anaReg.AlignAllBrains = [];
%         Align.started = datestr(now);
%         
%         try
%             AlignAllBrains(SaveFolder)
%             
%         catch e
%             disp(['Align all brains error' DataFolder])
%             Align.error = datestr(now);
%             Align.def = e;
%             anaReg.AlignAllBrains = Align;
% %             throw(e)
%         end
%         anaReg.AlignAllBrains = Align;
%         clear Align
%     else
%         disp('Alignment all brains already done')
%     end
%     
end
