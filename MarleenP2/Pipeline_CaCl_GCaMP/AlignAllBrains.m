function AlignAllBrains(DataFolder)

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

indfileseps = strfind(DataFolder, filesep); %zoek alle plekken van fileseps in naam
pathROI = DataFolder(1:indfileseps(end-2));

%% check if already done
if exist([DataFolder filesep 'tformAlignAllBrains.mat'], 'file')
    
    % if tform was already calculated, dont redo it unless you redid the ROIs
    TimeROI = dir([pathROI 'BigROI.mat']);
    TimeROI = datetime(TimeROI.date);
    
    Timetform = dir([DataFolder filesep 'tformAlignAllBrains.mat']);
    Timetform = datetime(Timetform.date);
    
    if Timetform > TimeROI
        disp('tform already calculated')
        return
    end
end

clear indfileseps TimeROI Timetform 

%% coregistration
% Make sure the mice brains are in the same location by matching them to
% your BigROIRef
% Get general ROI, fixed image to align to
load('/media/mbakker/GDrive/P2/GCaMP/BigROI_ReferenceMask.mat', 'AtlasMask');
BigROIRef = AtlasMask;
BigROIRefsmall = BigROIRef(135:345, 135:390);
clear AtlasMask;

% Get ROI and mask of the current acquisition
load([pathROI 'BigROI.mat'], 'AtlasMask');
AtlasMasksmall = AtlasMask(135:345, 135:390);
        
%% Try 1: automatic
%compute how much to shift
[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumIterations = 1e10;
optimizer.MinimumStepLength = 1e-4;
optimizer.MaximumStepLength = 1e-2;
optimizer.GradientMagnitudeTolerance = 1e-3;
tform = imregtform(AtlasMasksmall, BigROIRefsmall, 'similarity', optimizer, metric,...
    'DisplayOptimization', false, 'PyramidLevels', 3);

% Confirmation
AtlasMasktemp = imwarp(AtlasMask,tform,'OutputView',imref2d(size(BigROIRef)));
% ssimval = ssim(AtlasMasktemp(135:345, 135:390), BigROIRef(135:345, 135:390));

f1 = figure;
imshowpair(AtlasMasktemp, BigROIRef);
answer = questdlg('Does it make sense?', ...
    'Coregistration', ...
    'Yes','No','Yes');
close(f1)

switch answer
    case 'Yes'
        save([DataFolder filesep 'tformAlignAllBrains.mat'], 'tform');
        
%% Try 2: Manual
    case 'No'
        %% brunos app
        f = msgbox('save as tformAlignAllBrains');
        manualAlignFrames(BigROIRef, AtlasMask);
        
%         %% Sam's app
%         tform = AlignFrame_Manual(BigROIRefsmall, AtlasMasksmall);
%         
%         AtlasMasktemp = imwarp(AtlasMask,tform,'OutputView',imref2d(size(BigROIRef)));
%         f1 = figure;
%         imshowpair(AtlasMasktemp, BigROIRef);
%         
%         answer2 = questdlg('Does it make sense?', ...
%             'Coregistration', ...
%             'Yes','No','Yes');
%         close(f1)
%         
%         switch answer2
%             case 'Yes'
%                 save([DataFolder filesep 'tformAlignAllBrains.mat'], 'tform');
%             case 'No'
%                 disp('tform not saved')
%                 return
%         end

end

if ~exist([DataFolder filesep 'tformAlignAllBrains.mat'], 'file')
    error('Did not save tformAlignAllBrains')
end

end