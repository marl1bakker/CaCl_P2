% choose one pixel in the brain, and make a correlationmap of that.
% ManualInput 0 - dont ignore existing files
% ManualInput 1 - overwrite existing files

%% P2
% Input Seedname can be several, give like {'seed1', 'seed2'} etc.
%{'VisualROI_R'} SensoryROI_R AuditoryROI_R UnknownROI_R MotorROI_R

%goes with BigROI standard, because the small rois will give way too many
%points.

function CombinedSPCM(dataname, Acquisition, Grouping, GSR, Overwrite)

SaveDir = '/media/mbakker/GDrive/P2/GCaMP';

if ~exist('Acquisition', 'var')
    Acquisition = 'A1';
end

disp(['Combined SPCM ' Acquisition])

if ~exist('Overwrite', 'var')
    Overwrite = 1;
end

if ~exist('GSR', 'var')
    GSR = 0;
end

if ~exist('dataname', 'var')
    dataname = 'hemoCorr_fluo.dat';
elseif  length(dataname) < 4 || ( ~strcmp(dataname(end-3:end), '.dat') )
    dataname = [dataname '.dat'];
end

% if it already exists, dont do it
if GSR == 1
    SaveFolderSPCM = '/SPCM/Combined/withGSR/';
elseif GSR == 0
    SaveFolderSPCM = '/SPCM/Combined/withoutGSR/' ;
end

%don't do if already done
if ( exist([SaveDir SaveFolderSPCM dataname '_Difference_' Acquisition '_VisualROI_R.tiff'], 'file') ...
        && Overwrite == 0 )
    disp('Seed pixel correlation map already done, function exited')
    return
elseif( exist([SaveDir SaveFolderSPCM dataname '_Difference_' Acquisition '_VisualROI_R.tiff'], 'file') ...
        && Overwrite == 1 )
    disp('Seed pixel correlation map already done, OVERWRITING FILES')
end


%% Choose grouping
if ~exist('Grouping', 'var')
    Possibilities = {'Group', 'Sex'};
    [indx, ~] = listdlg('ListString', Possibilities);
    Grouping = Possibilities(indx);
end

load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat', 'RecordingOverview');

[groups, RecordingOverview] = GroupVariables(RecordingOverview, Grouping);

if size(Grouping, 2)>1
    Grouping = 'Combi';
end

% load('/media/mbakker/GDrive/P2/GCaMP/BigROI_ReferenceMask.mat', 'regions');

for indgroup = 1:size(groups,2) % Go per group
    group = groups{indgroup};
    eval([group 'N = 0'])

    disp(group);
    eval(['idx = RecordingOverview.' Grouping ' == ''' group ''';'])
    Mousegroup = RecordingOverview(idx,:);
    AllSPCM = NaN(size(Mousegroup,1), 10, 512*512); % mice x roi x spcm(rho)
    
    for indmouse = 1:size(Mousegroup, 1) %go per mouse
        Mouse = Mousegroup.Mouse{indmouse};
        disp(Mouse)
        eval(['DataFolder = [Mousegroup.' Acquisition '{indmouse} filesep];']);
        SaveFolder = [Mousegroup.SaveDirectory{indmouse} filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];

        if GSR == 0 && exist(['/media/mbakker/GDrive/P2/GCaMP/SPCM/Individual/' Mouse '-' Acq '-SPCM.mat'], 'file')
            load(['/media/mbakker/GDrive/P2/GCaMP/SPCM/Individual/' Mouse '-' Acq '-SPCM.mat'], 'Allrho')
        elseif GSR == 1 && exist(['/media/mbakker/GDrive/P2/GCaMP/SPCM/Individual/' Mouse '-' Acq '-SPCM-GSR.mat'], 'file')
            load(['/media/mbakker/GDrive/P2/GCaMP/SPCM/Individual/' Mouse '-' Acq '-SPCM-GSR.mat'], 'Allrho')
        else
            Allrho = SingleSubjectSPCM(SaveFolder, dataname, 0, GSR);
        end      
        
        Allrho = reshape(Allrho, 10, 512*512);
        AllSPCM(indmouse, :, :) = Allrho;
        
        eval([group 'N = ' group 'N +1;'])
    end %end indmouse
    
    %save group
    eval(['save([SaveDir SaveFolderSPCM ''' group '.mat''], ''AllSPCM'');'])
    eval(['SPCM_' group ' = mean(AllSPCM, 2, ''omitnan'');']);
    
end %end indgroup

clear Allrho 

%% Plot spcm's per group
for indgroup = 1:size(groups,2) % Go per group, nacl/cacl
    group = groups{indgroup};
    eval(['ngroup = ' group 'N;'])
    for indregion = 1:size(regions, 2)
        f = figure;
        eval(['imagesc(reshape(SPCM_' group '(indregion,1,:), 512, 512), [-1 1]);'])
        title([group ' ' Acquisition ' ' regions{indregion}], 'interpreter', 'none')
        subtitle(['n = ' num2str(ngroup)])
        axis off
        colorbar
        if GSR == 0
            f.CurrentAxes.CLim = [0 1];
            saveas(gcf, [SaveDir '/SPCM/Combined/withoutGSR/' dataname '_' group '_' Acquisition '_' regions{indregion} '.tiff'])
        else
            f.CurrentAxes.CLim = [-1 1];
            saveas(gcf, [SaveDir '/SPCM/Combined/withGSR/' dataname '_' group '_' Acquisition '_' regions{indregion} '.tiff'])
        end
        %         pause
        close(f)
    end %of regions
    
end %end of saving spcm's per group

%% Save SPCM, difference
for indregion = 1:size(regions,2)
    f = figure;
    SPCM_diff = SPCM_NaCl(indregion, 1, :) - SPCM_CaCl(indregion, 1, :);
    SPCM_diff = reshape(SPCM_diff, 512, 512);
    imagesc(SPCM_diff)
    title(['Difference ' Acquisition ' ' regions{indregion}], 'interpreter', 'none')
    subtitle(['N CaCl = ' num2str(CaClN) ' --- N NaCl = ' num2str(NaClN)])
    axis off
    colorbar
    if GSR == 0
        f.CurrentAxes.CLim = [0 1];
        saveas(gcf, [SaveDir '/SPCM/Combined/withoutGSR/' dataname '_Difference_' Acquisition '_' regions{indregion} '.tiff'])
    else
        f.CurrentAxes.CLim = [0 1];
        saveas(gcf, [SaveDir '/SPCM/Combined/withGSR/' dataname '_Difference_' Acquisition '_' regions{indregion} '.tiff'])
    end
    close(f)
end

end % end function