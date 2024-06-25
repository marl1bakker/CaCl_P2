% datatype is 'fluo' or 'speckle'

function Correct_for_rotation_ROI(DataFolder, datatype)

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

if ~exist('datatype', 'var') || matches(datatype, 'fluo') 
    datatype = 'fluo_567';
    nroffolders = 2;
    % GCaMP has an extra CtxImg folder in the acq, whereas speckle does
    % not. Need to account for this in getting the pathROI
elseif matches(datatype, 'speckle')
    nroffolders = 1;
else
    disp('Datatype not recognized, put fluo or speckle')
    return
end

idx = strfind(DataFolder, filesep); %zoek alle plekken van fileseps in naam
pathROI = DataFolder(1:idx(end-nroffolders));

if ~exist([pathROI 'BigROI.mat'], 'file')
    disp('No BigROI to correct.')
    return
end

%% check if already done
infovar = who('-file', [pathROI 'BigROI.mat']);
if sum(contains(infovar,'Rotation'))
    disp('rotation already done, exited function')
    return
end

%% get data and tform
load([pathROI 'ImagingReferenceFrame.mat'], 'reference_frame')

dims = load([DataFolder datatype '.mat'], 'datSize');
dims = dims.datSize;

if matches(datatype, 'fluo_567')
    fid = fopen([DataFolder datatype '.dat']);
    fl = fread(fid, dims(1)*dims(2), '*single');
    fl = reshape(fl, dims(1), dims(2));
    fclose(fid);
elseif matches(datatype, 'speckle')
    load([DataFolder 'meanBFI.mat'], 'meanBFI')
    % P = prctile(meanBFI(:), [0.1 99.9]);
    P = prctile(meanBFI(:),[0.1 99]);
    BrainIm = (meanBFI-P(2))./(P(1)-P(2));
    BrainIm(BrainIm(:)<0) = 0;
    BrainIm(BrainIm(:)>1) = 1;
    load([pathROI 'Masks.mat'], 'VesselMask');
    BrainIm = BrainIm .* VesselMask;
    fl = BrainIm;

    % temp:
    f1 = figure;
    imagesc(fl, [0 1])
    f2 = figure;
    imagesc(reference_frame, [0 1]) % if not the same, do different P
    f1.Position = [20 20 700 600];
    f2.Position = [800 20 700 600];

    diffp = questdlg('Are the two images (slightly) different in values?');
    close(f1, f2)
    if matches(diffp, 'Yes')
        load([DataFolder 'meanBFI.mat'], 'meanBFI')
        P = prctile(meanBFI(:), [0.1 99.9]);
        BrainIm = (meanBFI-P(2))./(P(1)-P(2));
        BrainIm(BrainIm(:)<0) = 0;
        BrainIm(BrainIm(:)>1) = 1;
        load([pathROI 'Masks.mat'], 'VesselMask');
        BrainIm = BrainIm .* VesselMask;
        fl = BrainIm;

        f1 = figure;
        imagesc(fl, [0 1])
        f2 = figure;
        imagesc(reference_frame, [0 1]) % if not the same, do different P
        f1.Position = [20 20 700 600];
        f2.Position = [800 20 700 600];

        diffpcheck = questdlg('Correct like this?');
        if matches(diffpcheck, 'No')
            disp('Something wrong with Reference Image')
            disp(DataFolder)
            return
        end
        close(f1, f2)
    end
    clear BrainIm meanBFI f1 f2 
end

[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumIterations = 1e6;
optimizer.MinimumStepLength = 1e-5;
optimizer.MaximumStepLength = 1e-2;
optimizer.GradientMagnitudeTolerance = 1e-3;

tform = imregtform(reference_frame, fl, 'affine', optimizer, metric,...
    'DisplayOptimization', false, 'PyramidLevels', 3);

%% validate
f1 = figure;
f1.Name = 'old';
f1.Position = [20 20 700 600];
imshowpair(reference_frame, fl);
dat_corr = imwarp(reference_frame,tform,'OutputView',imref2d(size(fl)));
f2 = figure;
f2.Name = 'new';
f2.Position = [800 20 700 600];
imshowpair(fl, dat_corr);

answer = questdlg('Is coregistration correct?');

close(f1, f2)

if ~matches(answer, 'Yes')
    % disp('DO THIS ONE BY HAND')
    tform = AlignFrame_Manual(fl, reference_frame);
end

clear answer dat_corr f1 f2 fid idx metric nroffolders optimizer

%% Load BigROI and translate
load([pathROI 'BigROI.mat'], 'AtlasMask', 'BigROI', 'regions');

% AtlasMask
AtlasMask = imwarp(AtlasMask,tform,'OutputView',imref2d(size(fl)));
roundatlas = round(AtlasMask); % to get rid of numbers with decimals
AtlasMask(AtlasMask~=roundatlas) = 0;

%BigROI
BigROI_old = BigROI;
for indROI = 1:size(regions,2)
    eval(['MaskROI = BigROI.' regions{indROI} ';'])
    MaskROI = imwarp(MaskROI,tform,'OutputView',imref2d(size(fl)));
    roundROI = round(MaskROI);
    MaskROI(MaskROI~=roundROI) = 0;

    eval(['BigROI.' regions{indROI} ' = MaskROI;'])
end

% regions stay the same. Give mark that you rotated
Rotation = 1;

%% Save
save([pathROI 'BigROI.mat'], 'BigROI', 'AtlasMask', 'regions', 'Rotation');

% %temp
% if matches(datatype, 'speckle') && matches(diffp, 'Yes')
%     load([DataFolder 'meanBFI.mat'], 'meanBFI')
%     P = prctile(meanBFI(:),[0.1 99.9]);
%     BrainIm = (meanBFI-P(2))./(P(1)-P(2));
%     BrainIm(BrainIm(:)<0) = 0;
%     BrainIm(BrainIm(:)>1) = 1;
%     load([pathROI 'Masks.mat'], 'VesselMask');
%     BrainIm = BrainIm .* VesselMask;
% 
%     reference_frame = imwarp(BrainIm,tform,'OutputView',imref2d(size(fl)));
%     save([pathROI 'ImagingReferenceFrame.mat'], 'reference_frame', '-append')
% end
end

function TOut = AlignFrame_Manual(ImFixed, ImVar)
h = figure('Position', [100 100 750 550], 'CloseRequestFcn', @CloseFig);

%ImFixed = double(imread(ImFixedPath));
%ImVar = double(imread(ImVarPath));
% 
% if( contains(ImVarPath, 'og','IgnoreCase',true) )
%     ImVar = fliplr(ImVar);
% end
TOut = [1 0 0; 0 1 0; 0 0 1];
Pf = prctile(ImFixed(:), [10 90]);
Pv = prctile(ImVar(:), [10 90]);

ImFixed = (ImFixed - Pf(1))/(Pf(2) - Pf(1));
ImVar = (ImVar - Pv(1))/(Pv(2) - Pv(1));
ImFixed(ImFixed(:)<0) = 0;
ImFixed(ImFixed(:)>1) = 1;
ImVar(ImVar(:)<0) = 0;
ImVar(ImVar(:)>1) = 1;

ImFixed = adapthisteq(ImFixed);
ImVar = adapthisteq(ImVar);
NewIm = ImVar;

ax = axes('Parent', h, 'Position', [0.2 0.2 0.75 0.75]);
imshowpair(ImFixed, ImVar, 'Parent', ax)
% Alteration Marleen: added Sliderstep, adjusted min and max so that one
% arrow click is one pixel
% Y axis
hsV = uicontrol('Parent', h, 'Position', [5 5 25 510], 'Style', 'slider', 'Min', -size(ImFixed,1), 'Max', size(ImFixed,1), 'SliderStep', [1/size(ImFixed,1) 0.01],'Value', 0, 'Callback', @MoveFrame);
uicontrol('Parent', h, 'Position', [5 510 25 30], 'Style', 'text', 'String', 'Y');
% X axis
hsH = uicontrol('Parent', h, 'Position', [40 5 25 510], 'Style', 'slider', 'Min', -size(ImFixed,2), 'Max', size(ImFixed,2), 'SliderStep', [1/size(ImFixed,2), 0.01],'Value', 0, 'Callback', @MoveFrame);
uicontrol('Parent', h, 'Position', [40 510 25 30], 'Style', 'text', 'String', 'X');
% rotation
hsR = uicontrol('Parent', h, 'Position', [75 5 25 510], 'Style', 'slider', 'Min', -pi/4, 'Max', pi/4, 'SliderStep', [1/size(ImFixed,2), 0.01], 'Value', 0, 'Callback', @MoveFrame);
uicontrol('Parent', h, 'Position', [75 510 25 30], 'Style', 'text', 'String', char(hex2dec('398')));
% scaling
hsS = uicontrol('Parent', h, 'Position', [110 5 25 510], 'Style', 'slider', 'Min', 0, 'Max', 2, 'SliderStep', [1/size(ImFixed,2), 0.01], 'Value', 1, 'Callback', @MoveFrame);
uicontrol('Parent', h, 'Position', [110 510 25 30], 'Style', 'text', 'String', 'S');

hDisp = uicontrol('Parent', h, 'Position', [200 20 75 50], 'Style', 'popupmenu', 'String', 'Tout|Fixe|Mobile|Alternance', 'Callback', @MoveFrame);
uicontrol('Parent', h, 'Position', [120 20 80 50], 'Style', 'text', 'String', 'Affichage:');
uicontrol('Parent', h, 'Position', [300 20 75 50], 'Style', 'pushbutton', 'String', 'Sauvegarde', 'Callback', @Save);


TimerObj = [];
bImShowed = 0;

MoveFrame();
waitfor(h);

    function MoveFrame(~,~,~)
        Rdefaut =  imref2d(size(ImVar));
        tX = mean(Rdefaut.XWorldLimits);
        tY = mean(Rdefaut.YWorldLimits);
        offX = hsH.Value;
        offY = hsV.Value;
        offR = hsR.Value;
        scale = hsS.Value;
        tScale = [scale, 0, 0; 0, scale, 0; 0, 0, 1];
        tTranslationToCenterAtOrigin = [1 0 0; 0 1 0; -tX -tY,1];
        tTranslationBackToOriginalCenter = [1 0 0; 0 1 0; tX tY,1];
        tRotation = [cos(offR) -sin(offR) 0; sin(offR) cos(offR) 0; 0 0 1];
        tTranslation = [1 0 0; 0 1 0; -offX -offY,1];
        tformCenteredRotation = tTranslationToCenterAtOrigin*tRotation*tTranslationBackToOriginalCenter*tTranslation*tScale;
        tformCenteredRotation = affine2d(tformCenteredRotation);
        TOut = tformCenteredRotation;

        NewIm = imwarp(ImVar, tformCenteredRotation, 'OutputView',imref2d(size(ImFixed)));
        
        if( ~isempty(TimerObj) )
            stop(TimerObj);
            delete(TimerObj);
            TimerObj = [];
        end
        
        switch hDisp.Value
            case 1
                imshowpair(ImFixed, NewIm, 'Parent', ax);
            case 2
                imagesc(ax, ImFixed);
                axis image;
                colormap gray;
            case 3
                imagesc(ax, NewIm);
                axis image;
                colormap gray;
            case 4
                TimerObj = timer('Period', 0.75, 'ExecutionMode', 'FixedRate', 'TimerFcn', @TimerUpdate);
                imagesc(ax, ImFixed);
                axis image;
                colormap gray;
                bImShowed = 0;
                start(TimerObj);
        end

        
    end

    function TimerUpdate(~,~,~)
        if( bImShowed )
            imagesc(ax, ImFixed);
            axis image;
                colormap gray;
            bImShowed = 0;
        else
            imagesc(ax, NewIm);
            axis image;
                colormap gray;
            bImShowed = 1;
        end
    end

    function CloseFig(~,~,~)
         Rdefaut =  imref2d(size(ImVar));
        tX = mean(Rdefaut.XWorldLimits);
        tY = mean(Rdefaut.YWorldLimits);
        offX = hsH.Value;
        offY = hsV.Value;
        offR = hsR.Value;
        scale = hsS.Value;
        tScale = [scale, 0, 0; 0, scale, 0; 0, 0, 1];
        tTranslationToCenterAtOrigin = [1 0 0; 0 1 0; -tX -tY,1];
        tTranslationBackToOriginalCenter = [1 0 0; 0 1 0; tX tY,1];
        tRotation = [cos(offR) -sin(offR) 0; sin(offR) cos(offR) 0; 0 0 1];
        tTranslation = [1 0 0; 0 1 0; -offX -offY,1];
        tformCenteredRotation = tTranslationToCenterAtOrigin*tRotation*tTranslationBackToOriginalCenter*tTranslation*tScale;
        tformCenteredRotation = affine2d(tformCenteredRotation);
        TOut = tformCenteredRotation;

        if( ~isempty(TimerObj) )
            stop(TimerObj);
            delete(TimerObj);
        end
        delete(h);

    end
        
    function Save(~,~,~)

    end
end
