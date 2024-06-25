% datatype is 'speckle', 'HbO', or 'HbR'

function RightvsLeft_MakeTable(DataFolder, datatype, manualinput, overwrite)

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

if ~exist('datatype', 'var')
    datatype = 'speckle';
end

if ~exist('manualinput', 'var')
    manualinput = 0;
end

% get mouse info
idx = strfind(DataFolder, filesep); 
if matches(datatype, 'speckle')
    pathROI = DataFolder(1:idx(end-1));
    Mouse = DataFolder(idx(end-2)+1:idx(end-1)-1);
    Acq = DataFolder(idx(end-1)+1:idx(end-1)+2);

    ratiotablepath = '/media/mbakker/GDrive/P2/Speckle/LR_ratio/LR_ratio_table.mat';
else
    pathROI = DataFolder(1:idx(end-2));
    Mouse = DataFolder(idx(end-3)+1:idx(end-2)-1);
    Acq = DataFolder(idx(end-2)+1:idx(end-2)+2);

    ratiotablepath = ['/media/mbakker/GDrive/P2/GCaMP/LR_ratio/LR_ratio_table_' datatype '.mat'];
end

% get table for ratios
disp('get table')
if ~exist(ratiotablepath, 'file')
    sz = [1 7];
    vartypes = ["string", "categorical", "categorical", "categorical", "string", "string", "double"];
    varnames = ["Mouse", "Group", "Sex", "Combi", "Acq", "ROI", "Ratio"];
    LR_ratio_table = table('Size', sz, 'VariableTypes', vartypes, 'VariableNames', varnames);
else
    load(ratiotablepath, 'LR_ratio_table')
end

% check if already done
if ~exist('overwrite', 'var')
    overwrite = 0;
end

Mousetable = LR_ratio_table(matches(LR_ratio_table.Mouse, Mouse),:);
if size(Mousetable(matches(Mousetable.Acq, Acq),:),1) > 0 && overwrite == 0
    disp(['RightvsLeft already done for Mouse ' Mouse ' Acquisition ' Acq])
    return
elseif size(Mousetable(matches(Mousetable.Acq, Acq),:),1) > 0 && overwrite == 1
    disp(['OVERWRITING RightvsLeft for Mouse ' Mouse ' Acquisition ' Acq])
    LR_ratio_table(matches(LR_ratio_table.Mouse, Mouse),:) = [];
end
clear Mousetable


%% Make Mousetable
% Get all the right info, make a table for which you only change the ROI
% and ratio values as you go through the script
load('/home/mbakker/P2_scripts/MarleenP2/MiceCodes.mat', 'Mice')
MouseInfo = Mice(matches(Mice.CodeOfMouse, Mouse),:);

Mousetable = table;
Mousetable.Mouse = Mouse;
Mousetable.Group = MouseInfo.CaClSham;
Mousetable.Sex = MouseInfo.MaleFemale;
Mousetable.Combi = strcat(string(MouseInfo.CaClSham), string(MouseInfo.MaleFemale));
Mousetable.Acq = Acq;
clear Mice


%% check if right things are done
% temp:
if exist([DataFolder 'OutlierMask.mat'], 'file') 
    infovar = who('-file', [DataFolder 'OutlierMask.mat']);
else
    infovar = 'none';
end
if ~sum(contains(infovar, 'OutlierPixels')) && manualinput == 0
    disp(['MakeOutlierMask not done yet for ' Mouse])
    return
elseif ~sum(contains(infovar, 'OutlierPixels')) && manualinput == 1
    MakeOutlierMask(DataFolder, manualinput);
end

% temp:
infovar = who('-file', [pathROI 'BigROI.mat']);
if sum(contains(infovar,'Rotation')) == 0 
    disp('do rotation of ROI before making this')
    return
end

%% get data
disp('open data')
tic
load([DataFolder 'MovMask.mat'], 'MovMask');
load([pathROI 'BigROI.mat'], 'BigROI', 'AtlasMask');
load([DataFolder 'OutlierMask.mat'], 'OutlierPixels', 'OutlierFrames')

% shift atlas if it's another acquisition
if ~matches(Acq, 'A1')
    load([DataFolder 'tform_acq.mat'], 'tform')
    % load([DataFolder 'tform.mat','tform'])
    invtform = invert(tform);
    Rois = fieldnames(BigROI);
    for indroi = 1:size(Rois, 1)
        eval(['BigROI.' Rois{indroi} '= imwarp(BigROI.' Rois{indroi} ', invtform, ''OutputView'', imref2d(size(AtlasMask)), ''interp'', ''nearest'');'])
    end
end

if matches(datatype, 'speckle')
    fid = fopen([DataFolder 'BFI.dat']);
    dat = fread(fid, inf, '*single');
    fclose(fid);
    toc
    dat = reshape(dat,1024*1024, []);
    % get rid of 10 first and last frames, is effect of BFI calculation
    dat(:,1:10) = NaN;
    dat(:,end-10:end) = NaN;

    % get masks
    load([pathROI 'Masks.mat'], 'LeftMask', 'RightMask', 'VesselMask');
    LeftMask = LeftMask .* VesselMask;
    RightMask = RightMask .* VesselMask;
else
    fid = fopen([DataFolder datatype '.dat']);
    dat = fread(fid, inf, '*single');
    fclose(fid);
    toc
    dat = reshape(dat,512*512, []);

    % get masks
    LeftMask = BigROI.VisualROI_L + BigROI.AuditoryROI_L + BigROI.MotorROI_L...
        + BigROI.RetrosplenialROI_L + BigROI.SensoryROI_L;
    se = strel('disk', 5);
    LeftMask = imclose(LeftMask, se);

    RightMask = BigROI.VisualROI_R + BigROI.AuditoryROI_R + BigROI.MotorROI_R...
        + BigROI.RetrosplenialROI_R + BigROI.SensoryROI_R;
    se = strel('disk', 5);
    RightMask = imclose(RightMask, se);
    clear se
    
    % to make sure you dont divide by (close to) 0 for the ratio, so add
    % estimated baseline values for hbo/hbr
    if matches(datatype, 'HbO')
        dat = dat+60;
    elseif matches(datatype, 'HbR')
        dat = dat+40;
    end
end

%% Movement and Outlier removal
dat(:, MovMask(:)==0) = NaN;
dat(:, OutlierFrames(:)==0) = NaN;
eval(['dat(OutlierPixels.' datatype '(:)==0,:) = NaN;']);

%% left vs right
disp('masks left right')

avright = mean(dat(RightMask(:)>0,:),"all", 'omitnan');
avleft = mean(dat(LeftMask(:)>0, :), "all", 'omitnan');
Mousetable.ROI = 'Hemisphere';
Mousetable.Ratio = avright/avleft;
LR_ratio_table = [LR_ratio_table; Mousetable];

% if you want over time...
% avright = mean(dat(RightMask(:)>0,:),1, 'omitnan');
% avleft = mean(dat(LeftMask(:)>0, :), 1, 'omitnan');
% ratio = avright./avleft;

clear LeftMask RightMask avright avleft 

%% BigROI
disp('masks roi')
Rois = {'Visual', 'Auditory', 'Sensory', 'Motor', 'Retrosplenial'};

for indroi = 1:size(Rois, 2)
    ROI = Rois{indroi};

    eval(['LeftMask = BigROI.' ROI 'ROI_L;'])
    eval(['RightMask = BigROI.' ROI 'ROI_R;'])

    if matches(datatype, 'speckle')
        LeftMask = LeftMask .* VesselMask;
        RightMask = RightMask .* VesselMask;
    end

    avright = mean(dat(RightMask(:)>0,:),"all", 'omitnan');
    avleft = mean(dat(LeftMask(:)>0, :), "all", 'omitnan');

    Mousetable.ROI = ROI;
    Mousetable.Ratio = avright/avleft;
    LR_ratio_table = [LR_ratio_table; Mousetable];
end

%% save
disp('saving')
save(ratiotablepath, 'LR_ratio_table');

end

