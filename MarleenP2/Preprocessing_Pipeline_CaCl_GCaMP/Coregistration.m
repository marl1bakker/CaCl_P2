%P2
% This function is used to coregister different acquisitions of the same
% mouse with each other, so that the brain is at the same place. It needs
% the function CoregistrationManual. 

% changed minor things 7-2-2024, not tested. Added ToBeCorrected input,
% made filesize adaptable, not hardcoded on 512.

% practice:
% path_fixed =  '/media/mbakker/GDrive/P2/GCaMP/M13/A1-R2/CtxImg';
% path_moving = '/media/mbakker/GDrive/P2/GCaMP/M13/A2-R1/CtxImg';

% path_fixed is the folder of the mouse with data to base your fixation on,
% in this case A1. path_moving are the datafolders with the files that you
% want to coregister, in this case A2 and A3. dataname is the file you base
% your coregistration on, in this case green or speckle. ToBeCorrected are
% the files that you want to move. In case of GCaMP, that is green, red,
% and fluo, in case of speckle it is speckle. Give like {'green', 'red', 'fluo_567.dat'}

function Coregistration(path_fixed, path_moving, dataname, ToBeCorrected)

if ~exist('dataname', 'var')
    dataname = 'green.dat';
elseif size(dataname,2)<4 || ~matches(dataname(end-3:end), '.dat')
    dataname = [dataname '.dat'];
end

if( ~strcmp(path_moving(end), filesep) )
    path_moving = [path_moving filesep];
end
if( ~strcmp(path_fixed(end), filesep) )
    path_fixed = [path_fixed filesep];
end

%% Step 1: Load data
load([path_fixed dataname(1:end-4) '.mat'], 'datSize');

fid = fopen([path_fixed dataname]);
dat = fread(fid, datSize(1)*datSize(2), '*single');
dat = reshape(dat, datSize(1),datSize(2));
fclose(fid);
datfixed = dat;

fid = fopen([path_moving dataname]);
dat = fread(fid, datSize(1)*datSize(2), '*single');
dat = reshape(dat, datSize(1),datSize(2));
fclose(fid);
datmov = dat;
clear dat

%% Step 3: Prep, make differences more clear
datfixed = datfixed./mean(datfixed(:));
datfixed = (datfixed-min(datfixed(:)))./(max(datfixed(:))-min(datfixed(:)));
datfixed(datfixed < 0) = 0;
datfixed(datfixed > 1) = 1;
datfixed = adapthisteq(datfixed);

datmov = datmov./mean(datmov(:));
datmov = (datmov-min(datmov(:)))./(max(datmov(:))-min(datmov(:)));
datmov(datmov < 0) = 0;
datmov(datmov > 1) = 1;
datmov = adapthisteq(datmov);

%% Step 4: CoReg
[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumIterations = 1e6;
optimizer.MinimumStepLength = 1e-5;
optimizer.MaximumStepLength = 1e-2;
optimizer.GradientMagnitudeTolerance = 1e-3;
tform = imregtform(datmov, datfixed, 'affine', optimizer, metric,...
    'DisplayOptimization', false, 'PyramidLevels', 3);

%% Step 5:Confirmation
f1 = figure;
f1.Name = 'old';
imshowpair(datfixed, datmov);
dat_corr = imwarp(datmov,tform,'OutputView',imref2d(size(datfixed)));
f2 = figure;
f2.Name = 'new';
imshowpair(datfixed, dat_corr);
% saveas(gcf,[path_moving 'CoregistrationComp.png']);
% 
pause 
answer = questdlg('Does it make sense?', ...
	'Coregistration', ...
	'Yes','No','Yes');
% Handle response
switch answer
    case 'Yes'
        disp('≧◠‿●‿◠≦    ᕙ(^▿^-ᕙ)');
        close(f1, f2)
    case 'No'
        disp('xxx')
%         return
        close(f1, f2)
        tform = AlignFrame_Manual(datfixed, datmov);

        % error('Coregistration not good enough to apply')
end
disp('applying coregistration...')

%% Step 5 B: Save/load
save([path_moving 'tform_acq.mat'],'tform');

% %% Step 6: Apply to all images
% % if( bApply )
% if ~exist('ToBeCorrected', 'var')
%     ToBeCorrected = {'green.dat', 'red.dat', 'fluo_567.dat'};
% end
% 
% for index = 1:size(ToBeCorrected,2)
%     dataname = ToBeCorrected{index};
%     fid = fopen([path_moving dataname]);
%     dat = fread(fid, inf, '*single');
%     dat = reshape(dat, datSize(1),datSize(2),[]);
%     fclose(fid);
% 
%     for ind = 1:size(dat,3)
%         dat(:,:,ind) = imwarp(squeeze(dat(:,:,ind)),tform,'OutputView',imref2d(size(dat(:,:,1))),'interp','nearest');
%     end
% 
%     fid = fopen([path_moving dataname], 'w');
%     fwrite(fid,dat,'*single');
% 	fclose(fid);   
% end


end %of function



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