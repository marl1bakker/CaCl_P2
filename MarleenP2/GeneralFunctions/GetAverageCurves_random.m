function GetAverageCurves_random(DataFolder, datanamefluocurves, NVCname, datanamehbo, datanamehbr)
%% set up
if ~strcmp(DataFolder(end), filesep) 
    DataFolder = [DataFolder filesep];
end

if ~exist('datanameactivations', 'var')
    datanameactivations = 'hemoCorr_fluo';
end

if ~exist('datanamefluocurves', 'var')
    datanamefluocurves = 'hemoCorr_fluo';
end

if ~exist('NVCname', 'var') && matches(datanamefluocurves, 'hemoCorr_fluo')
    NVCname = 'NVC_ROI_random';
elseif ~exist('NVCname', 'var')
    NVCname = ['NVC_ROI_random_' datanamefluocurves];
end

if ~exist('datanamehbo', 'var')
    datanamehbo = 'HbO';
    datanamehbr = 'HbR';
end

if ~exist([DataFolder datanameactivations '.dat'], 'file')
    disp([datanameactivations ' could not be found, function exited'])
    return
end

if ~exist([DataFolder datanamehbo '.dat'], 'file')
    disp([datanamehbo ' could not be found, function exited'])
    return
end


%% timewindow
imfreq = 15;
timebefore = 5*imfreq; %5 sec before
timeafter = 10*imfreq; % 10 sec after
windowtime = timebefore+timeafter+1;

%% get points of activation
% random points in this case, but same number of points as original data
if exist('/media/mbakker/GDrive/P2/GCaMP/NVC/NrOfActivations.mat', 'file')
    load('/media/mbakker/GDrive/P2/GCaMP/NVC/NrOfActivations.mat', 'ActsTable')
else
    disp('Run NumberOfActivations first!')
end

seps = strfind(DataFolder, filesep);
Mouse = DataFolder(seps(end-3)+1:seps(end-2)-1);
Acq = DataFolder(seps(end-2)+1:seps(end-2)+2);
ActsTable = ActsTable(matches(ActsTable.Mouse, Mouse),:);
ActsTable = ActsTable(matches(ActsTable.Acq, Acq),:);

%ongeveer 100 activations per pixel

%% load data
fid = fopen([DataFolder datanamefluocurves '.dat']);
dF = fread(fid, inf, '*single');
dF = reshape(dF, 512,512,[]);
fclose(fid);

% new 19-7-24 get rid of outliers before detecting activations
% seps = strfind(DataFolder, filesep);
load([DataFolder(1:seps(end-2)) 'ROImasks_data.mat'], 'img_info');
mask = img_info.logical_mask;
load([DataFolder 'OutlierMask.mat'], 'OutlierFrames', 'OutlierPixels');
dF = reshape(dF, 512*512, []);
dF = dF .* OutlierFrames;
dF = reshape(dF, 512, 512, []);
dF = dF .* OutlierPixels.(datanameactivations);
dF = dF .* mask;
dF(dF == 0) = NaN;

fid = fopen([DataFolder datanamehbo '.dat']);
dH = fread(fid, inf, '*single');
dH = reshape(dH,size(dF));
fclose(fid);

fid = fopen([DataFolder datanamehbr '.dat']);
dR = fread(fid, inf, '*single');
dR = reshape(dR,size(dF));
fclose(fid);

dT = dH + dR;

clear datanamehbr datanamefluo fid

%% Get mask for brain and movement and outliers
seps = strfind(DataFolder, filesep);
load([DataFolder(1:seps(end-2)) 'ROImasks_data.mat'], 'img_info');
mask = img_info.logical_mask;

load([DataFolder 'OutlierMask.mat'], 'OutlierFrames', 'OutlierPixels');
OutlierPix = OutlierPixels.HbO + OutlierPixels.HbR + OutlierPixels.hemoCorr_fluo;
OutlierPix(OutlierPix<3) = 0;
OutlierPix(OutlierPix == 3) = 1;

mask = mask.*logical(OutlierPix);
mask = logical(mask);

%take only brain
aF = reshape(aF,[], size(aF,3));
aF = aF(mask(:),:);
dF = reshape(dF,[], size(dF,3));
dat.dF = dF(mask(:),:);
dH = reshape(dH,[], size(dH,3));
dat.dH = dH(mask(:),:);
dR = reshape(dR,[], size(dR,3));
dat.dR = dR(mask(:),:);
dT = reshape(dT,[], size(dT,3));
dat.dT = dT(mask(:),:);
clear dF dH dR dT

%exclude movement and outliers
load([DataFolder 'MovMask.mat'], 'MovMask');
aF = aF.* MovMask;
aF = aF.*OutlierFrames;

clear img_info MovMask seps OutlierMask

datatypes = {'fluo', 'hbo', 'hbr', 'hbt'};
datanames = {'dF', 'dH', 'dR', 'dT'};

%% get average curves per pixel
for inddata = 1:size(datatypes, 2)
    avcurvesperpixel = NaN(size(aF,1), windowtime);

    for pixel = 1:size(aF, 1)
        % get activations of pixel, but exclude beginning and end (where you
        % cant get a full plot)
        indact = find(aF(pixel,:));
        indact(indact<=(timebefore)) = [];
        indact(indact>=(size(aF,2)-timeafter)) = [];

        if size(indact,2)<1 %if theres no activations
            continue
        end

        temp = NaN(size(indact,2), windowtime);

        for ind2 = 1:size(indact,2)
            activation = indact(ind2); %get curves for all activations of this pixel
            temp(ind2,:) = dat.(datanames{inddata})(pixel,activation-timebefore:activation+timeafter);
        end

        avcurvesperpixel(pixel,:) = mean(temp,1,'omitnan');
    end
    % Curves whole brain
    % eval([datatypes{inddata} 'curves.WholeBrain = mean(avcurvesperpixel, 1, ''omitnan'');'])
    curves.(datatypes{inddata}).WholeBrain = mean(avcurvesperpixel, 1, 'omitnan');

    % put all curves per pixel in right format/size (512*512, 216)
    temp = NaN(512*512, windowtime);
    temp(mask,:) = avcurvesperpixel;
    % eval(['avcurvesperpixel_' datatypes{inddata} '= temp;'])
    PixelCurves.(datatypes{inddata}) = temp;
    disp([datatypes{inddata} ' done']);
end
clear activation aF fluotemp hbotemp hbrtemp pixel ind2 indact dF dH dR avcurvesperpixel dat datanames inddata temp

% fluocurves = curves.fluo;
% hbocurves = curves.hbo;
% hbrcurves = curves.hbr;
% hbtcurves = curves.hbt;
% clear curves
% avcurvesperpixel_hbt = avcurvesperpixel_hbo+avcurvesperpixel_hbr;
% hbtcurves.WholeBrain = mean(avcurvesperpixel_hbt, 1, 'omitnan');
%% Whole Brain
% Specs Whole Brain (curves done before)
[tempSpecs] = GetSpecs(curves.fluo, curves.hbo, curves.hbr, 'WholeBrain', timebefore, imfreq, datanameactivations);
if ~exist('Specs', 'var')
    Specs = tempSpecs;
else
    Specs = [Specs; tempSpecs];
end
clear tempSpecs

%% Specific areas (ROI)
% Curves ROI
seps = strfind(DataFolder, filesep);
load([DataFolder(1:seps(end-2)) 'BigROI.mat'], 'regions', 'BigROI', 'AtlasMask');

for indroi = 1:size(regions,2) %go per roi
    % eval(['mapROI = BigROI.' regions{indroi} ';']) %changed, did not test
    mapROI = BigROI.(regions{indroi});
    mapROI = reshape(mapROI, 512*512, []);

    for inddata = 1:size(datatypes, 2) %plot per datatype
        % eval(['ROIdata = avcurvesperpixel_' datatypes{inddata} '(mapROI,:);'])
        % eval([datatypes{inddata} 'curves.' regions{indroi} ' = mean(ROIdata,1, ''omitnan'');'])
        ROIdata = PixelCurves.(datatypes{inddata})(mapROI,:);
        curves.(datatypes{inddata}).(regions{indroi}) = mean(ROIdata,1, 'omitnan');
    end

    % Specs ROI
    [tempSpecs] = GetSpecs(curves.fluo, curves.hbo, curves.hbr, regions{indroi}, timebefore, imfreq, datanameactivations);
    Specs = [Specs; tempSpecs];
    clear tempSpecs
end
clear regions BigROI indroi inddata seps

%% Half brain
% Curves half brain
AtlasMask(AtlasMask == 0) = NaN;

rightmask = zeros(size(AtlasMask));
rightmask(AtlasMask < 6) = 1;
rightmask = bwmorph(rightmask, 'close');
rightmask = imclose(rightmask, strel('disk', 5));

leftmask = zeros(size(AtlasMask));
leftmask(AtlasMask > 5) = 1;
leftmask = bwmorph(leftmask, 'close');
leftmask = imclose(leftmask, strel('disk', 5));

temp = zeros(size(AtlasMask));
temp(rightmask) = 1;
temp(leftmask) = 2;
AtlasMask = temp;

leftmask = leftmask.*OutlierPix;
rightmask = rightmask.*OutlierPix;
leftmask = logical(reshape(leftmask, 512*512,1));
rightmask = logical(reshape(rightmask, 512*512,1));

for inddata = 1:size(datatypes, 2)
    % eval(['left = avcurvesperpixel_' datatypes{inddata} '(leftmask,:);'])
    % eval([datatypes{inddata} 'curves.Left = mean(left,1,''omitnan'');'])
    % 
    % eval(['right = avcurvesperpixel_' datatypes{inddata} '(rightmask,:);'])
    % eval([datatypes{inddata} 'curves.Right = mean(right,1,''omitnan'');'])
    left = PixelCurves.(datatypes{inddata})(leftmask,:);
    curves.(datatypes{inddata}).Left = mean(left,1,'omitnan');
    right = PixelCurves.(datatypes{inddata})(rightmask,:);
    curves.(datatypes{inddata}).Right = mean(right,1,'omitnan');
end
clear temp inddata seps avcurve_fluo avcurve_hbo avcurve_hbr

% Specs half brain
[tempSpecs] = GetSpecs(curves.fluo, curves.hbo, curves.hbr, 'Left', timebefore, imfreq, datanameactivations);
Specs = [Specs; tempSpecs];
[tempSpecs] = GetSpecs(curves.fluo, curves.hbo, curves.hbr, 'Right', timebefore, imfreq, datanameactivations);
Specs = [Specs; tempSpecs];
clear tempSpecs

%% save
fluocurves = curves.fluo;
hbocurves = curves.hbo;
hbrcurves = curves.hbr;
hbtcurves = curves.hbt;

% Specs(matches(Specs.ROI, 'Dummy'),:) = [];
save([DataFolder NVCname '.mat'], 'fluocurves', 'hbocurves', 'hbrcurves', 'hbtcurves', 'Specs');
end



function [Specs] = GetSpecs(fluo, hbo, hbr, roiname, timebefore, imfreq, actname)
VarNames = {'ROI', 'DelaySec', 'DelayFrames', 'ResponseStrength', 'ResponseStrengthRelative'...
    'GCaMPDipBefore', 'GCaMPPeak', 'GCaMPDipAfter',...
    'HbODipBefore', 'HbOPeak', 'HbODipAfter', ...
    'HbRDipBefore', 'HbRPeak', 'HbRDipAfter'};
VarTypes =     {'cell', 'single', 'single', 'single', 'single', ...
    'single', 'single', 'single',...
    'single', 'single', 'single',...
    'single', 'single', 'single'};
Specs = table('Size', size(VarNames), 'VariableNames', VarNames, ...
    'VariableTypes', VarTypes);

Specs.ROI = cellstr(roiname);

% if isempty(maxfluo) %if there's only nan
if anynan(fluo.(roiname)) || anynan(hbo.(roiname))
    Specs.DelaySec = NaN;
    Specs.DelayFrames = NaN;
    Specs.ResponseStrength = NaN;
    return
end

if matches(actname, 'HbO')
    [maxhbo, indhbo] = findpeaks(hbo.(roiname)(timebefore:end));
    indhbo = timebefore + indhbo(1) - 1;
    maxhbo = maxhbo(1);

    [maxfluo, indfluo] = findpeaks(fluo.(roiname)(1:indhbo));
    indfluo = indfluo(end);
    maxfluo = maxfluo(end);
else % if its not based on hbo activations, assume its based on fluo
    [maxfluo, indfluo] = findpeaks(fluo.(roiname)(timebefore:end)); %find first peak after detected activation
    indfluo = timebefore + indfluo(1) - 1;
    maxfluo = maxfluo(1);

    [maxhbo, indhbo] = findpeaks(hbo.(roiname)(indfluo:end));
    indhbo = indfluo + indhbo(1) - 1;
    maxhbo = maxhbo(1);
end


% GCaMPPeak
Specs.GCaMPPeak = maxfluo;
indDips = islocalmin(fluo.(roiname));
indDipBefore = find(indDips(1:indfluo), 1, 'last');
if isempty(indDipBefore)
    indDipBefore = 1;
end
Specs.GCaMPDipBefore = fluo.(roiname)(indDipBefore);
indDipAfter = indDipBefore + find(indDips(indDipBefore+1:end), 1, 'first');
if isempty(indDipAfter)
    indDipAfter = size(indDips, 2);
end
Specs.GCaMPDipAfter = fluo.(roiname)(indDipAfter);

%HbOPeak
Specs.HbOPeak = maxhbo;
indDips = islocalmin(hbo.(roiname));
indDipBefore = find(indDips(1:indhbo), 1, 'last');
if isempty(indDipBefore)
    indDipBefore = 1;
end
Specs.HbODipBefore = hbo.(roiname)(indDipBefore);
indDipAfter = indDipBefore + find(indDips(indDipBefore+1:end), 1, 'first');
if isempty(indDipAfter)
    indDipAfter = size(indDips, 2);
end
Specs.HbODipAfter = hbo.(roiname)(indDipAfter);

%HbRPeak
[maxhbr, indhbr] = findpeaks(hbr.(roiname)(indfluo:end));
maxhbr = maxhbr(1);
indhbr = indhbr(1) + indfluo - 1;
Specs.HbRPeak = maxhbr;
indDips = islocalmin(hbr.(roiname));
indDipBefore = find(indDips(1:indhbr(1)), 1, 'last');
if isempty(indDipBefore)
    indDipBefore = 1;
end
Specs.HbRDipBefore = hbr.(roiname)(indDipBefore);
indDipAfter = indDipBefore + find(indDips(indDipBefore+1:end), 1, 'first');
if isempty(indDipAfter)
    indDipAfter = size(indDips, 2);
end
Specs.HbRDipAfter = hbr.(roiname)(indDipAfter);

%Delays
Specs.DelaySec = (indhbo - indfluo)/imfreq; %in seconds
Specs.DelayFrames = indhbo - indfluo;

% Response Strengths
maxfluo = (maxfluo-1) *100; % in percentage
Specs.ResponseStrength = maxhbo / maxfluo;
increasefluo = Specs.GCaMPPeak - Specs.GCaMPDipBefore;
increasehbo = Specs.HbOPeak - Specs.HbODipBefore;
Specs.ResponseStrengthRelative = increasehbo/increasefluo;

end
