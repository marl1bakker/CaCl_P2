function Masks_and_ROI(DataFolder, overwrite, ManualInput)

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end
idx = strfind(DataFolder, filesep); 
pathROI = DataFolder(1:idx(end-1));
Mouse = DataFolder(idx(end-2)+1:idx(end-1)-1);

if ~exist('overwrite', 'var')
    overwrite = 0;
end

if exist([pathROI 'Masks.mat'], 'file') && overwrite == 0
    disp('Masks already made, function exited')
    return
elseif exist([pathROI 'Masks.mat'], 'file') && overwrite == 1
    disp('OVERWRITING MASKS')
elseif ~exist([DataFolder 'BFI.dat'], 'file')
    disp('Cannot make masks, run CalculateBFI first.')
    return
elseif exist('ManualInput', 'var') && ManualInput == 0
    disp('Cannot make masks without manual input...')
    return
end

clear idx overwrite ManualInput

%% Get data
% Temp, if you run things again, CalculateBFI will save the meanBFI.
if ~exist([DataFolder 'meanBFI.mat'], 'file')
    fid = fopen([DataFolder 'BFI.dat']);
    BFI = fread(fid, inf, '*single');
    fclose(fid);
    BFI = reshape(BFI, 1024, 1024, []);
    meanBFI = mean(BFI, 3);
    save([DataFolder 'meanBFI.mat'], 'meanBFI');
    clear BFI
else
    load([DataFolder 'meanBFI.mat'], 'meanBFI')
end

%% Get mask for vessels
% If there are damages at the window, there are a lot of values that are
% very high. This would give an unusable BrainIm, if we take 99.9
% percentile because the extremely high values would still be there. If
% this is seen (manually) in a mouse, we add it to the exception cases and
% take only 99 percentile to get the BrainIm to be useful.
% ExceptionCases = {'M27', 'M3', 'M4', 'M5'}; 

% if sum(matches(ExceptionCases, Mouse)) > 0
    P = prctile(meanBFI(:), [0.1 99]);
% else
%     P = prctile(meanBFI(:),[0.1 99.9]);
% end
BrainIm = (meanBFI-P(2))./(P(1)-P(2));
BrainIm(BrainIm(:)<0) = 0;
BrainIm(BrainIm(:)>1) = 1;

options.FrangiBetaTwo = 100; %15 og % 100
options.FrangiScaleRatio = 2; %2
options.FrangiScaleRange = [1 8]; %[1 8]
options.FrangiBetaOne = 0.5; %0.5
options.verbose = false; %true
options.BlackWhite = true;

[outIm, ~, ~] = FrangiFilter2D(BrainIm, options);
VesselMask = outIm>0;
VesselMask = bwareaopen(VesselMask, 20);
VesselMask = ~VesselMask; % make sure that vessels are 0, not 1

f = figure;
imagesc(BrainIm.*VesselMask)
answer = questdlg('Is this a good vessel mask?', 'Vessels', 'Yes', 'No', 'Yes');
close(f)

f2 = figure;
while matches(answer, 'No')
    prompts = {'Give new FrangiBetaTwo (lower nr is more vessels covered)', ...
        'Give new FrangiScaleRatio (lower nr is thicker vessels)',...
        'New FrangiScaleRange Lower', 'NewFrangiScaleRange Higher'};
    dlgtitle = 'Frangi values';
    fieldsize = [1 45; 1 45; 1 45; 1 45];
    definputs = {num2str(options.FrangiBetaTwo), num2str(options.FrangiScaleRatio),...
        num2str(options.FrangiScaleRange(1)), num2str(options.FrangiScaleRange(2))};

    newoptions = inputdlg(prompts, dlgtitle, fieldsize, definputs);
    options.FrangiBetaTwo = str2double(newoptions{1});
    options.FrangiScaleRatio = str2double(newoptions{2});
    options.FrangiScaleRange = [str2double(newoptions{3}), str2double(newoptions{4})];

    [outIm, ~, ~] = FrangiFilter2D(BrainIm, options);
    VesselMask = outIm>0;
    VesselMask = bwareaopen(VesselMask, 20);
    VesselMask = ~VesselMask; % make sure that vessels are 0, not 1

    imagesc(BrainIm.*VesselMask)
    answer = questdlg('Is this a good vessel mask?', 'Vessels', 'Yes', 'No', 'Yes');
end

close(f2)

BrainIm = BrainIm .* VesselMask;

save([pathROI 'Masks.mat'], 'VesselMask', 'options');


%% Make left and right mask for brain, manually
% f = figure;
% imagesc(BrainIm)
% title('Draw right mask')
% poly = drawpolygon;
% Mask.Right = poly.Position;
% Mask.Right = poly2mask(Mask.Right(:,1), Mask.Right(:,2), 1024, 1024);
% close(f)
%
% f = figure;
% imagesc(BrainIm)
% title('Draw left mask')
% poly = drawpolygon;
% Mask.Left = poly.Position;
% Mask.Left = poly2mask(Mask.Left(:,1), Mask.Left(:,2), 1024, 1024);
% close(f)


%% Fit Allen ROI on brain
cd(pathROI)
figure
imagesc(meanBFI, [1 prctile(meanBFI, 99, 'all')])

ROImanager(BrainIm)
f = msgbox(["In ROImanager, do the following steps:"; ...
    "-  Image – Set origin – New – drag to bregma, right mouse-click to set"; ...
    "-  Image – Set origin – Align image to origin – drag to lambda, right mouse-click to set and correct for tilt in frame"; ...
    "-  Image – Set pixel size – 110 pixels per mm"; ...
    "-  Image – Mask – Draw new – Draw the mask by clicking points, you can still drag the points after setting them, double click inside the mask to confirm"; ...
    "-  Image – Image reference file… - Export – save as ImagingReferenceFrame.mat in the mouse folder"; ...
    "-  Create – Mouse Allen Brain Atlas – Areas – Select areas – close window to select them all – drag atlas to what seems to be fitting – double click inside atlas to confirm"; ...
    "-  File – Save as… - ROImasks_data.mat in mouse folder"]);
input('Make reference and ROI')

ClusterRois(DataFolder, 1, 1) %overwrite, and single folder (no ctximg like gcamp)
Correct_for_rotation_ROI(DataFolder, 'speckle')

%% Make left and right side of the brain out of allen roi
load([pathROI 'BigROI.mat'])

LeftMask = BigROI.VisualROI_L + BigROI.AuditoryROI_L + BigROI.MotorROI_L...
    + BigROI.RetrosplenialROI_L + BigROI.SensoryROI_L;
se = strel('disk', 5);
LeftMask = imclose(LeftMask, se);

RightMask = BigROI.VisualROI_R + BigROI.AuditoryROI_R + BigROI.MotorROI_R...
    + BigROI.RetrosplenialROI_R + BigROI.SensoryROI_R;
se = strel('disk', 5);
RightMask = imclose(RightMask, se);
Rotation = 1;

%% Save masks
save([pathROI 'Masks.mat'], 'RightMask', 'LeftMask', 'Rotation', '-append');

end
