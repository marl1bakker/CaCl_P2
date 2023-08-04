% cd '/media/mbakker/GDrive/P2/GCaMP/M13/A1-R2/CtxImg'
cd '/media/mbakker/GDrive/P2/GCaMP/M15/A1-R2/CtxImg'
% [~, ~] = HemoComputeNoFilter(pwd, pwd, 'GCaMP', {'red','green'}, 1);
% dF = NormalisationFiltering(pwd,'fluo',0.3, 3, 1,0);
% fid = fopen('mov_aligned.dat');
fid = fopen('hemoCorr_fluo.dat');
dF = fread(fid, inf, '*single');
dF = reshape(dF, 512,512,[]);

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

% fid = fopen('HbO.dat');
fid = fopen([DataFolder 'HbONoFilt.dat']);
dH = fread(fid, inf, '*single');
dH = reshape(dH,size(dF));
clear zF;
%% Segment the brain area
% imagesc(sum(aF,3));
% h = drawpolygon;
% mask = h.createMask;
% load('/media/mbakker/GDrive/P2/GCaMP/M13/ROImasks_data.mat')
load('/media/mbakker/GDrive/P2/GCaMP/M15/ROImasks_data.mat')
mask = img_info.logical_mask;

%% Reducing memory usage
dH = reshape(dH,[], size(dH,3));
aF = reshape(aF,[], size(aF,3));
dH = dH(mask(:),:);
aF = aF(mask(:),:);
%%
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
parfor( indI = 1:size(aF,1), 5) %5 is nr cores
    indexes = idx;
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

map = zeros(512, 512);
map(mask) = mapAct;
mapAct = map;
map(mask) = mapTime;
mapTime = map;
map(mask) = mapAmp;
mapAmp = map;

figure
imagesc(mapAct, [0 140])
title('Activity map - m15')

figure
imagesc(mapTime, [0 2])
title('Delay in response - m15')

figure
imagesc(mapAmp, [0 8])
title('Amplitude map - m15')

load('MovMask.mat')
act = aF .* MovMask;

dF = reshape(dF, 512, 512, []);
dH = reshape(dH, 512, 512, []);
pixeldF = reshape(dF(300,200,:), 1, []);
pixeldH = reshape(dH(300,200,:), 1, []);

figure
tiledlayout(3,1)
nexttile
plot(pixeldF)
nexttile
plot(pixeldH)
nexttile
plot(MovMask)

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