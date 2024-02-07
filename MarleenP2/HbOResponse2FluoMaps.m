function HbOResponse2FluoMaps(DataFolder, datanamefluo, datanamehbo)

% DataFolder =  '/media/mbakker/GDrive/P2/GCaMP/M13/A1-R2/CtxImg'
% DataFolder =  '/media/mbakker/GDrive2/P2/GCaMP/M15/A1-R2/CtxImg'
% DataFolder =  '/media/mbakker/GDrive2/P2/GCaMP/M20/A1-R2/CtxImg';

%% set up
if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

if ~exist('datanamefluo', 'var')
    datanamefluo = 'hemoCorr_fluo';
end

if ~exist('datanamehbo', 'var')
    datanamehbo = 'HbO';
end

if ~exist([DataFolder datanamefluo '.dat'], 'file')
    disp([datanamefluo ' could not be found, function exited'])
    return
end

if ~exist([DataFolder datanamehbo '.dat'], 'file')
    disp([datanamehbo ' could not be found, function exited'])
    return
end

%% load data
fid = fopen([DataFolder datanamefluo '.dat']);
dF = fread(fid, inf, '*single');
dF = reshape(dF, 512,512,[]);

fid = fopen([DataFolder datanamehbo '.dat']);
dH = fread(fid, inf, '*single');
dH = reshape(dH,size(dF));

%% get points of activation
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

clear zF;
%% Get mask for brain and movement
seps = strfind(DataFolder, filesep);
load([DataFolder(1:seps(end-2)) 'ROImasks_data.mat'], 'img_info');
mask = img_info.logical_mask;

load([DataFolder 'MovMask.mat'], 'MovMask');

%% Reducing memory usage
%take only brain
dH = reshape(dH,[], size(dH,3));
aF = reshape(aF,[], size(aF,3));
dH = dH(mask(:),:);
aF = aF(mask(:),:);

%exclude movement
dH = dH.* MovMask;
aF = aF.* MovMask;

%% Make map of number of activations per pixel
idx = find(mask);
mapAct = zeros(size(aF,1),1, 'single');
mapAmp = zeros(size(aF,1),1, 'single');
mapTime = zeros(size(aF,1),1, 'single');
mapAcc = zeros(size(idx,1),300, 'single');

aF(:,1:151) = 0;
aF(:,(end-150):end) = 0;
for indI = 1:size(aF,1)
    tmp = find(aF(indI, :));
    act{indI} = tmp;
    mapAct(indI) = length(tmp); %hoe vaak gaat pixel over threshold?
end
map = zeros(512, 512);
map(mask) = mapAct;
mapAct = map;


%%
parfor( indI = 1:size(aF,1), 5) %5 is nr cores
    indexes = idx; %what falls within mask
    Pixels = dH; %to make parfor fast
    Acc = zeros(1,300);
    tAct = act{indI};
    [iY, iX] = ind2sub([512, 512], indexes(indI));
    subsY = iY + (-5:5); %neighbour pixels
    subsY = repelem(subsY,11);
    subsX = iX + (-5:5);
    subsX = repmat(subsX,1, 11);
    indices = sub2ind([512, 512], subsY, subsX);
    indices = find(ismember(indexes, indices));
    subPix = Pixels(indices,:);
    if isempty(tAct) %for if there was no activation above trheshold for that pixel
        mapAcc(indI,:) = 0;
        mapAmp(indI) = 0;
        mapTime(indI) = 0;
    else
        for iAct = 1:length(tAct)
            St = tAct(iAct) - 150; %50 is for 10s at 5 Hz
            En = St + 299;
            tmp = subPix(:, St:En);
            tmp = mean(tmp- mean(tmp(:,120:150),2), 1); %2 sec before act.
            Acc = Acc + tmp;
        end
        Acc = Acc./iAct;
        mapAcc(indI,:) = Acc;
        [A, T] = max(Acc(150:225)); %onset to 5 sec after
        mapAmp(indI) = A;
        mapTime(indI) = T/15;
    end
end

map = zeros(512, 512);
map(mask) = mapTime;
mapTime = map;

map(mask) = mapAmp;
mapAmp = map;



map = zeros(512*512, size(aF, 2));
map(mask(:), :) = aF;
aF = reshape(map, 512, 512, []);
% for ind = 1:500
%     imagesc(aF(:,:,ind), [0 1])
%     pause(0.05)
% end

map = zeros(512*512, size(dH, 2));
map(mask(:), :) = dH;
dH = reshape(map, 512, 512, []);
% for ind = 1:500
%     imagesc(dH(:,:,ind), [-15 15])
%     pause(0.05)
% end



figure
imagesc(mapAct, [0 140])
title('Activity map - m15')

figure
imagesc(mapTime, [0 2])
title('Delay in response - m15')

figure
imagesc(mapAmp, [0 8])
title('Amplitude map - m15')


% dF = reshape(dF, 512, 512, []);
% dH = reshape(dH, 512, 512, []);
pixeldF = reshape(dF(300,200,:), 1, []);
pixeldH = reshape(dH(300,200,:), 1, []);
pixelaF = reshape(aF(300,200,:), 1, []);

figure
tiledlayout(4,1)
nexttile
plot(pixeldF)
title('fluo')
nexttile
plot(pixeldH)
title('hbo')
nexttile
plot(MovMask)
title('movement')
nexttile
plot(pixelaF)
title('activity')

%%
map = zeros(512,512,300);
figure;
for ind=1:300
    map(mask,ind) = mapAcc(:,ind);
    imagesc(map, [-5 5]);
    pause(0.05)
end

map = zeros(512*512,300);
for ind=1:300
    map(mask,ind) = mapAcc(:,ind);
end

map = reshape(map, 512,512,[]);
Time = linspace(-10,10,300);
plot(Time, squeeze(map(400,300,:)))

% hold on

dF = reshape(dF,[], size(dF,3));
dF = dF(mask(:),:);

parfor( indI = 1:size(aF,1), 5) %5 is nr cores
    indexes = idx;
    Pixels = dF; %to make parfor fast
    Acc = zeros(1,300);
    tAct = act{indI};
    [iY, iX] = ind2sub([512, 512], indexes(indI));
    subsY = iY + (-5:5); %neighbour pixels
    subsY = repelem(subsY,11);
    subsX = iX + (-5:5);
    subsX = repmat(subsX,1, 11);
    indices = sub2ind([512, 512], subsY, subsX);
    indices = find(ismember(indexes, indices));
    subPix = Pixels(indices,:);
    for iAct = 1:length(tAct)
        St = tAct(iAct) - 150; %50 is for 10s at 5 Hz
        En = St + 299;
        tmp = subPix(:, St:En);
        tmp = mean(tmp- mean(tmp(:,120:150),2), 1); %2 sec before act.
        Acc = Acc + tmp;
    end
    Acc = Acc./iAct;
    mapAccF(indI,:) = Acc;
    [A, T] = max(Acc(150:225)); %onset to 5 sec after
    mapAmp(indI) = A;
    mapTime(indI) = T/15;
end

mapF = zeros(512*512,300);
for ind=1:300
    mapF(mask,ind) = mapAccF(:,ind);
end

mapF = reshape(mapF, 512,512,[]);
Time = linspace(-10,10,300);
figure
plot(Time, squeeze(mapF(400,300,:)))

figure
tiledlayout(4,1)
nexttile
plot(Time, squeeze(mapF(400,300,:)))
nexttile
plot(Time, squeeze(map(400,300,:)))

fid = fopen('HbRNoFilt.dat');
dR = fread(fid, inf, '*single');
dR = reshape(dR,512,512,[]);
dR = reshape(dR,[], size(dR,3));
dR = dR(mask(:),:);

parfor( indI = 1:size(aF,1), 5) %5 is nr cores
    indexes = idx;
    Pixels = dR; %to make parfor fast
    Acc = zeros(1,300);
    tAct = act{indI};
    [iY, iX] = ind2sub([512, 512], indexes(indI));
    subsY = iY + (-5:5); %neighbour pixels
    subsY = repelem(subsY,11);
    subsX = iX + (-5:5);
    subsX = repmat(subsX,1, 11);
    indices = sub2ind([512, 512], subsY, subsX);
    indices = find(ismember(indexes, indices));
    subPix = Pixels(indices,:);
    for iAct = 1:length(tAct)
        St = tAct(iAct) - 150; %50 is for 10s at 5 Hz
        En = St + 299;
        tmp = subPix(:, St:En);
        tmp = mean(tmp- mean(tmp(:,120:150),2), 1); %2 sec before act.
        Acc = Acc + tmp;
    end
    Acc = Acc./iAct;
    mapAccR(indI,:) = Acc;
    [A, T] = max(Acc(150:225)); %onset to 5 sec after
    mapAmp(indI) = A;
    mapTime(indI) = T/15;
end

mapR = zeros(512*512,300);
for ind=1:300
    mapR(mask,ind) = mapAccR(:,ind);
end
mapR = reshape(mapR, 512, 512, []);

nexttile
plot(Time, squeeze(mapR(400,300,:)))
nexttile
plot(Time, squeeze(mapR(400,300,:))+squeeze(map(400,300,:)))
%%