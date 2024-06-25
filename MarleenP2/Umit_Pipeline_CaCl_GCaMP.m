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

% Recordings = [RecordingOverview.A3];
% Mice = [RecordingOverview.Mouse];
% SaveDirs = [RecordingOverview.SaveDirectory];
% Recordings = [RecordingOverview.A2; RecordingOverview.A3];
% Mice = [RecordingOverview.Mouse; RecordingOverview.Mouse];
% SaveDirs = [RecordingOverview.SaveDirectory; RecordingOverview.SaveDirectory];

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

    % temporary, to get rid of stupid things
    anaReg = matfile([SaveFolder 'anaReg.mat'] ,'Writable',true);

    if ~isempty(who(anaReg, 'TimecourseFluo'))
        rmmatvar([SaveFolder 'anaReg.mat'], 'TimecourseFluo')
    end
    if ~isempty(who(anaReg, 'TimecourseHbOHbR'))
        rmmatvar([SaveFolder 'anaReg.mat'], 'TimecourseHbOHbR')
    end    
    if ~isempty(who(anaReg, 'Timecourses'))
        rmmatvar([SaveFolder 'anaReg.mat'], 'Timecourses')
    end


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

    clear Preproc_functions teller indfcndone teller

    RightvsLeft_MakeTable(SaveFolder, 'HbO', 0, 1)
    RightvsLeft_MakeTable(SaveFolder, 'HbR', 0, 1)

    %% check if this pipeline is already done (saves printing space)
    Umit_functions = {'HemoCorrection', 'NoFilt', 'Normalization', 'HbOHbR', ...
        'OutlierMask', 'TimecoursesFluo', 'TimecoursesHbOHbR'};
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

    clear Umit_functions teller indfcndone teller
    %% HemoCorrection
    disp('Hemodynamic correction...');
    varlist = who(anaReg,'HemoCorrection');

    if( isempty(varlist) || ~isfield(anaReg.HemoCorrection, 'ended') ||... %if not done in pipeline before
            datetime(getfield(anaReg.HemoCorrection, 'ended')) < datetime(getfield(anaReg.ImagesClassification, 'ended'))) 

        anaReg.HemoCorrection = [];
        hemcorr.started = datestr(now);
        try
            % need fluo.dat for this
            fid = fopen([SaveFolder 'fluo_567.dat']);
            dat = fread(fid, inf, '*single');
            metadata = load([SaveFolder 'fluo_567.mat']);
            hemoCorr_fluo = HemoCorrection(SaveFolder, dat, metadata, {'Red','Green'});
            fclose(fid);

            hemoCorr_fluo(isnan(hemoCorr_fluo)) = 0;

            fid = fopen([SaveFolder 'hemoCorr_fluo.dat'],'w');
            fwrite(fid, hemoCorr_fluo, '*single');
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

    %% no filt data
    disp('Save unfiltered fluo...');
    varlist = who(anaReg,'NoFilt');

    if( isempty(varlist) || ~isfield(anaReg.NoFilt, 'ended') ||... 
            datetime(getfield(anaReg.NoFilt, 'ended')) < datetime(getfield(anaReg.ImagesClassification, 'ended')))

        anaReg.NoFilt = [];
        nofilt.started = datestr(now);
        try
            fid = fopen([SaveFolder 'hemoCorr_fluo.dat']);
            dat = fread(fid, inf, '*single');
            dat(dat == 0) = NaN;
            dat = reshape(dat, 512, 512, []);
            fclose(fid);
            % normalize
            datnofilt = dat./mean(dat,3,'omitnan');

            fid = fopen([SaveFolder 'fluo_nofilt.dat'],'w');
            fwrite(fid,datnofilt,'*single');
            fclose(fid);

            nofilt.ended = datestr(now);
            disp('Unfiltered fluo saved');
        catch
            nofilt.error = datestr(now);
            anaReg.NoFilt = nofilt;
            disp(['No filt error ' DataFolder]);
            return;
        end
        anaReg.NoFilt = nofilt;
        clear nofilt;
    else
        disp('Unfiltered fluo already saved');
    end

    %% Normalize and filter
    disp('Normalization...');

    varlist = who(anaReg,'Normalization');
    if( isempty(varlist) || (~isfield(anaReg.Normalization, 'ended')) ||...
            datetime(getfield(anaReg.Normalization, 'ended')) < datetime(getfield(anaReg.ImagesClassification, 'ended'))) 
        if ~isfield(anaReg.NoFilt, 'ended')
            disp('Normalization/filtering skipped, save nofilt first')
            continue
        end

        anaReg.Normalization = [];
        Normalization.started = datestr(now);
        try
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
            % return;
            continue
        end
        anaReg.Normalization = Normalization;
        clear Normalization fid
    else
        disp('Normalization already done')
    end

    %% HbO HbR
    disp('HbO/HbR calculation...')
    varlist = who(anaReg,'HbOHbR');

    if( isempty(varlist) || (~isfield(anaReg.HbOHbR, 'ended')) ||...
            datetime(getfield(anaReg.HbOHbR, 'ended')) < datetime(getfield(anaReg.ImagesClassification, 'ended'))) 
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

    %% Outlier mask
    disp('Make outlier mask...');
    varlist = who(anaReg,'OutlierMask');

    if( isempty(varlist) || ~isfield(anaReg.OutlierMask, 'ended') ) %if not done in pipeline before

        anaReg.OutlierMask = [];
        outlier.started = datestr(now);
        try
            MakeOutlierMask(SaveFolder, 1);

            outlier.ended = datestr(now);
        catch e
            disp(['Outlier detection error ' DataFolder])
            outlier.error = datestr(now);
            outlier.def = e;
            anaReg.OutlierMask = outlier;
            throw(e)
        end

        anaReg.OutlierMask = outlier;
        clear outlier
    else
        disp('Outlier Mask already made')
    end

    %% Timecourses Fluo
    disp('Timecourses calculation Fluo...')

    varlist = who(anaReg,'TimecoursesFluo');

    if( isempty(varlist) || (~isfield(anaReg.TimecoursesFluo, 'ended')) ||... 
            datetime(getfield(anaReg.TimecoursesFluo, 'ended')) < datetime(getfield(anaReg.HemoCorrection, 'ended'))) 
        anaReg.TimecoursesFluo = [];
        Timecoursesvar.started = datestr(now);

        % check if you have all the necessary things
        if isfield(anaReg.HemoCorrection, 'ended') && isfield(anaReg.Normalization, 'ended')

            try %if yes, try to calculate timecourses
                GetTimecourses(SaveFolder, 'hemoCorr_fluo', ManualInput)
                Timecoursesvar.ended = datestr(now);

            catch e
                % disp(['Timcoureses Fluo error' DataFolder])
                Timecoursesvar.error = datestr(now);
                Timecoursesvar.def = e;
                %             anaReg.TimecoursesFluo = Timecoursesvar;
                %             throw(e)
            end

        else % if functions before are not done
            disp('Timecourses fluo not calculated, HemoCorr or Norm not done')
            Timecoursesvar.error = datestr(now);

        end

        anaReg.TimecoursesFluo = Timecoursesvar;
        clear Timecoursesvar
    else
        disp('Timecourses Fluo calculation already done')
    end

    %% Timecourses HbO/HbR
    disp('Timecourses calculation HbO/HbR...')

    varlist = who(anaReg,'TimecoursesHbOHbR');

    if( isempty(varlist) || (~isfield(anaReg.TimecoursesHbOHbR, 'ended')) ||... 
            datetime(getfield(anaReg.TimecoursesHbOHbR, 'ended')) < datetime(getfield(anaReg.HbOHbR, 'ended'))) 
        anaReg.TimecoursesHbOHbR = [];
        Timecoursesvar.started = datestr(now);

        % check if you have all the necessary things
        if isfield(anaReg.HbOHbR, 'ended')

            try %if yes, try to calculate timecourses
                GetTimecourses(SaveFolder, 'HbO', ManualInput)
                GetTimecourses(SaveFolder, 'HbR', ManualInput)
                Timecoursesvar.ended = datestr(now);

            catch e
                disp(['Timecourses HbO/HbR error' DataFolder])
                Timecoursesvar.error = datestr(now);
                Timecoursesvar.def = e;
                %             anaReg.TimecoursesHbOHbR = Timecoursesvar;
                %             throw(e)
            end

        else % if functions before are not done
            disp('Timecourses HbO/HbR not calculated, .dat files dont exist')
            Timecoursesvar.error = datestr(now);

        end

        anaReg.TimecoursesHbOHbR = Timecoursesvar;
        clear Timecoursesvar
    else
        disp('Timecourses HbOHbR calculation already done')
    end

    clear DataFolder Mouse SaveFolder varlist anaReg Acq

end

% function rmmatvar(matfile, varname)
% tmp = rmfield(load(matfile), varname);
% save(matfile, '-struct', 'tmp');
% end