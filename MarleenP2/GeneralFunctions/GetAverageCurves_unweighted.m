% WARNING, OUTLIER REMOVAL NOT IN SCRIPT YET

function GetAverageCurves_unweighted(DataFolder, datanamefluocurves, datanamehbo, datanamehbr, datanamefluoactivations)
%% set up
if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

if ~exist('datanamefluoactivations', 'var')
    datanamefluoactivations = 'hemoCorr_fluo';
end

if ~exist('datanamefluocurves', 'var')
    datanamefluocurves = 'hemoCorr_fluo';
end

if ~exist('datanamehbo', 'var')
    datanamehbo = 'HbO';
    datanamehbr = 'HbR';
end

if ~exist([DataFolder datanamefluoactivations '.dat'], 'file')
    disp([datanamefluoactivations ' could not be found, function exited'])
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
% load data
fid = fopen([DataFolder datanamefluoactivations '.dat']);
dF = fread(fid, inf, '*single');
dF = reshape(dF, 512,512,[]);
fclose(fid);

%Zscore:
zF = (dF - mean(dF, 3))./std(dF,0,3);
%Threshold on Zscore:
aF = zF >= 1.95;
% Removing noise:
for ind = 1:size(aF,3)
    aF(:,:,ind) = bwmorph(bwmorph(aF(:,:,ind), 'close', inf),'open',inf);
end
%Now, we want only the beginning of activations:
aF = aF(:,:,2:end)&~aF(:,:,1:(end-1));
aF = cat(3, false(size(aF,1),size(aF,2)), aF);

clear zF ind dF

% Get mask for brain and movement and outliers
seps = strfind(DataFolder, filesep);
load([DataFolder(1:seps(end-2)) 'ROImasks_data.mat'], 'img_info');
mask = img_info.logical_mask;

% %take only brain
aF = reshape(aF,[], size(aF,3));
aF = aF(mask(:),:);

%exclude movement and outliers
load([DataFolder 'MovMask.mat'], 'MovMask');
aF = aF.* MovMask;
load([DataFolder 'OutlierMask.mat'], 'OutlierMask');
aF = aF.*OutlierMask;

clear img_info MovMask seps OutlierMask

datanames = {datanamefluocurves datanamehbo datanamehbr};
datatypes = {'fluo', 'hbo', 'hbr'};
maxnrofact = max(sum(aF,2)); %what is the most activations per pixel

% needed for going per region/left/right
seps = strfind(DataFolder, filesep);
load([DataFolder(1:seps(end-2)) 'BigROI.mat'], 'regions', 'BigROI', 'AtlasMask');
[rightmask, leftmask] = GetHalfMasks(AtlasMask);
leftmask = double(leftmask(mask(:)));
leftmask(leftmask == 0) = NaN;
rightmask = double(rightmask(mask(:)));
rightmask(rightmask == 0) = NaN;
clear AtlasMask seps

%% get all activations per datatype
for inddata = 1:size(datanames,2)
    fid = fopen([DataFolder datanames{inddata} '.dat']);
    dat = fread(fid, inf, '*single');
    dat = reshape(dat,[], size(aF,2));
    fclose(fid);

    dat = dat(mask(:),:);

    % pixel x activations x windowtime
    AllAct = NaN(size(aF, 1), maxnrofact, windowtime);

    for pixel = 1:size(aF, 1)
        indact = find(aF(pixel,:));
        indact(indact<=(timebefore)) = [];
        indact(indact>=(size(aF,2)-timeafter)) = [];

        if size(indact,2)<1
            continue
        end

        for ind2 = 1:size(indact,2)
            activation = indact(ind2);
            AllAct(pixel, ind2, :) = dat(pixel,activation-timebefore:activation+timeafter);
        end
    end

    %% Curves whole brain
    eval([datatypes{inddata} 'curves.WholeBrain = reshape(mean(mean(AllAct, 1, ''omitnan''), 2, ''omitnan''), 1, windowtime);'])

    %% Curves per ROI
    for indroi = 1:size(regions,2) %go per roi
        eval(['mapROI = BigROI.' regions{indroi} ';'])
        mapROI = reshape(mapROI, 512*512, []);
        mapROI = mapROI(mask(:));
        mapROI = double(mapROI);
        mapROI(mapROI == 0) = NaN;
        ActROI = AllAct .* mapROI;

        eval([datatypes{inddata} 'curves.' regions{indroi} ' = reshape(mean(mean(ActROI, 1, ''omitnan''), 2, ''omitnan''), 1, windowtime);'])

        clear(datanames{inddata});
        clear ActROI pixel dat indact activation ind2
    end

    %% Curves per half brain
    % Curves half brain
    ActLeft = AllAct .* leftmask;
    ActLeft = reshape(mean(mean(ActLeft, 1, 'omitnan'), 2, 'omitnan'), 1, windowtime);
    eval([datatypes{inddata} 'curves.Left = ActLeft;'])

    ActRight = AllAct .* rightmask;
    ActRight = reshape(mean(mean(ActRight, 1, 'omitnan'), 2, 'omitnan'), 1, windowtime);
    eval([datatypes{inddata} 'curves.Right = ActRight;'])    

    clear ActLeft ActRight
    %% Save
    if exist([DataFolder 'NVC_ROI_unweighted.mat'], 'file')
        save([DataFolder 'NVC_ROI_unweighted.mat'], [datatypes{inddata} 'curves'], '-append');
    else
        save([DataFolder 'NVC_ROI_unweighted.mat'], [datatypes{inddata} 'curves']);
    end

    % delete([DataFolder 'NVC_ROI_fluo_nofilt.mat']);

    clear AllAct dat fid indroi mapROI 
    clear([datatypes{inddata} 'curves'])
end

%% Get specs
load([DataFolder 'NVC_ROI_unweighted.mat'], 'fluocurves', 'hbocurves', 'hbrcurves')

Specs = table('Size', [1 7], 'VariableNames', ...
    {'ROI', 'DelaySec', 'DelayFrames', 'ResponseStrength', 'GCaMPPeak', ...
    'HbOPeak', 'HbRPeak'}, 'VariableTypes', ...
    {'cell', 'single', 'single', 'single', 'single', 'single', 'single'});
Specs.ROI(1) = {'Dummy'};

% Whole brain
[tempSpecs] = GetSpecs(fluocurves, hbocurves, hbrcurves, 'WholeBrain', timebefore, imfreq);
Specs = [Specs; tempSpecs];
clear tempSpecs

% Specs ROI
for indroi = 1:size(regions,2) %go per roi
    [tempSpecs] = GetSpecs(fluocurves, hbocurves, hbrcurves, regions{indroi}, timebefore, imfreq);
    Specs = [Specs; tempSpecs];
    clear tempSpecs
end

% Specs half brain
[tempSpecs] = GetSpecs(fluocurves, hbocurves, hbrcurves, 'Left', timebefore, imfreq);
Specs = [Specs; tempSpecs];
[tempSpecs] = GetSpecs(fluocurves, hbocurves, hbrcurves, 'Right', timebefore, imfreq);
Specs = [Specs; tempSpecs];
clear tempSpecs

Specs(matches(Specs.ROI, 'Dummy'),:) = [];

%% save
save([DataFolder 'NVC_ROI_unweighted.mat'], 'Specs', '-append');

end


function [rightmask, leftmask] = GetHalfMasks(AtlasMask)
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

leftmask = reshape(leftmask, 512*512,1);
rightmask = reshape(rightmask, 512*512,1);

end

%%
function [Specs] = GetSpecs(fluo, hbo, hbr, roiname, timebefore, imfreq)
Specs = table('Size', [1 7], 'VariableNames', ...
    {'ROI', 'DelaySec', 'DelayFrames', 'ResponseStrength', 'GCaMPPeak', ...
    'HbOPeak', 'HbRPeak'}, 'VariableTypes', ...
    {'cell', 'single', 'single', 'single', 'single', 'single', 'single'});

Specs.ROI = roiname;

eval(['[maxfluo, indfluo] = findpeaks(fluo.' roiname '(timebefore:end));']) %find first peak after detected activation

if isempty(maxfluo) %if there's only nan
    Specs.DelaySec = NaN;
    Specs.DelayFrames = NaN;
    Specs.ResponseStrength = NaN;
    return
end

indfluo = timebefore + indfluo(1) - 1;
maxfluo = maxfluo(1);
Specs.GCaMPPeak = maxfluo;

eval(['[maxhbo, indhbo] = findpeaks(hbo.' roiname '(indfluo:end));'])
indhbo = indfluo + indhbo(1) - 1;
maxhbo = maxhbo(1);
Specs.HbOPeak = maxhbo;

eval(['[maxhbr, tilde] = findpeaks(hbr.' roiname '(indfluo:end));'])
maxhbr = maxhbr(1);
Specs.HbRPeak = maxhbr;

Specs.DelaySec = (indhbo - indfluo)/imfreq; %in seconds
Specs.DelayFrames = indhbo - indfluo;
maxfluo = (maxfluo-1) *100; % in percentage
Specs.ResponseStrength = maxhbo / maxfluo;

end