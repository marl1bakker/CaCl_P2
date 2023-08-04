% GCaMP_vs_HbO_plot
% Make a plot with GCaMP act on one axis and HbO on the other. One point
% corresponds to either one pixel or one ROI at one frame.
% DataFolder =  '/media/mbakker/GDrive/P2/GCaMP/M13/A1-R2/CtxImg'

function GCaMPvsHbOPlot(DataFolder, Option)


if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

if matches(Option, 'Marleen')
    %% Open data
fid = fopen([DataFolder 'hemoCorr_fluo.dat']);
dF = fread(fid, inf, '*single');
dF = reshape(dF, 512,512,[]);
fclose(fid);
snapshot = dF(:,:,1);

fid = fopen([DataFolder 'HbONoFilt.dat']);
dH = fread(fid, inf, '*single');
dH = reshape(dH,size(dF));
fclose(fid);
clear fid
%% Get activations of fluo
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

%% Get ROI data and mask
fileseps = strfind(DataFolder, filesep);
load([DataFolder(1:fileseps(end-2)) 'ROImasks_data.mat']);
mask = img_info.logical_mask;
clear img_info 

%% normalize fluo and dat
dF = reshape(dF, 512*512, []);
dF = dF./mean(dF, 2);

%% Reducing memory usage
dH = reshape(dH,[], size(dH,3));
aF = reshape(aF,[], size(aF,3));
dH = dH(mask(:),:);
aF = aF(mask(:),:);
dF = dF(mask(:),:);

aF(:,1:151) = 0;
aF(:,(end-150):end) = 0;

load([DataFolder 'MovMask.mat']);
aF = aF .* MovMask;
clear MovMask

%% whole brain
% hbomax = [];
% fluomax = [];
% 
% for ind = 1:(size(aF, 1)*size(aF, 2))
%     if aF(ind) == 1
%         curvehbo = dH(ind-30:ind+74);
%         [maxcurvehbo, ~] = max(curvehbo(30:end));
%         
%         curvefluo = dF(ind-30:ind+74);
%         [maxcurvefluo, ~] = max(curvefluo(30:end));
%         
%         hbomax(end+1) = maxcurvehbo;
%         fluomax(end+1) = maxcurvefluo;
%     end
% end
% 
% figure
% scatter(fluomax, hbomax)
% clear curvefluo curvehbo fluomax hbomax ind maxcurveflow maxcurvehbo

%% get dF & dH in 512, 512 again
map = zeros(512,512,size(aF, 2));
map = reshape(map, 512*512,[]);
map(mask,:) = dF;
dF = reshape(map, 512,512,[]);

map = zeros(512,512,size(aF, 2));
map = reshape(map, 512*512,[]);
map(mask,:) = dH;
dH = reshape(map, 512,512,[]);

map = zeros(512,512,size(aF, 2));
map = reshape(map, 512*512,[]);
map(mask, :) = aF;
aF = reshape(map, 512,512,[]);
clear map

%% single pixel
% figure
% imagesc(snapshot) % choose pixel
% pixel = [330 160]; % X Y
% 
% aFpixel = reshape(aF(pixel(2), pixel(1), :), 1,[]); %get crossing of threshold for this single pixel, have to get it with Y X
% activations = find(aFpixel == 1); % get the indices of the activations
% 
% hbomaxpixel = zeros(size(activations,2),1);
% fluomaxpixel = zeros(size(activations,2),1);
% 
% for ind = 1:size(activations, 2)
%     timepoint = activations(ind);
%     
%     curve = dH(pixel(2), pixel(1), timepoint-30:timepoint+74);
%     curve = reshape(curve, 1, []);
%     hbomaxpixel(ind) = max(curve(30:end));
%     
%     curve = dF(pixel(2), pixel(1), timepoint-30:timepoint+74);
%     curve = reshape(curve, 1, []);
%     fluomaxpixel(ind) = max(curve(30:end));
% end
% 
% figure
% scatter(hbomaxpixel, fluomaxpixel)
% 
% clear aFpixel activations curve fluomaxpixel hbomaxpixel ind hbomaxpixel fluomaxpixel pixel timepoint

%% Regions of interest
hbomaxpixel = [];
fluomaxpixel = [];

figure
hold on
for ind = 1:size(ROI_info,2)
    pixel = ROI_info(ind).PixelSeeds;
    aFpixel = reshape(aF(pixel(2), pixel(1), :), 1,[]);
    activations = find(aFpixel == 1);
    
    for index = 1:size(activations, 2)
        timepoint = activations(index);
        
        curve = dH(pixel(2), pixel(1), timepoint-30:timepoint+74);
        curve = reshape(curve, 1, []);
        hbomaxpixel(end+1) = max(curve(30:end));
        
        curve = dF(pixel(2), pixel(1), timepoint-30:timepoint+74);
        curve = reshape(curve, 1, []);
        fluomaxpixel(end+1) = max(curve(30:end));
    end
    
    scatter(fluomaxpixel, hbomaxpixel);
    
    hbomaxpixel = [];
    fluomaxpixel = [];
end
    
title([DataFolder(fileseps(end-3)+1:fileseps(end-1)-1)])
    
    
    
    
    
    
    
 


    %% Make maps, work with parfor
    
elseif matches(Option, 'Sam')
    
    
if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end
    
    fid = fopen([DataFolder 'hemoCorr_fluo.dat']);
    dF = fread(fid, inf, '*single');
    dF = reshape(dF, 512,512,[]);
    
    fid = fopen([DataFolder 'HbONoFilt.dat']);
    dH = fread(fid, inf, '*single');
    dH = reshape(dH,size(dF));
    
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
    
    clear zF dF
    
    %% Segment the brain area
    load('/media/mbakker/GDrive/P2/GCaMP/M15/ROImasks_data.mat')
    mask = img_info.logical_mask;
    
    %% Reducing memory usage
    dH = reshape(dH,[], size(dH,3));
    aF = reshape(aF,[], size(aF,3));
    dH = dH(mask(:),:);
    aF = aF(mask(:),:);
    
    %% hbo
    aF(:,1:151) = 0;
    aF(:,(end-150):end) = 0;
    
    mapAct = zeros(size(aF,1),1, 'single');
    for indI = 1:size(aF,1)
        tmp = find(aF(indI, :));
        act{indI} = tmp;
        mapAct(indI) = length(tmp); %hoe vaak gaat pixel over threshold?
    end
    
    idx = find(mask);
    mapAmp = zeros(size(aF,1),1, 'single');
    mapTime = zeros(size(aF,1),1, 'single');
    mapAcc = zeros(size(idx,1),105, 'single');
%     allevents = zeros(size(aF,1), size(aF,2), 105, 'single');
    clear aF
    
    parfor( indI = 1:size(dH,1), 4) %5 is nr cores
        
        indexes = idx;
        Pixels = dH; %to make parfor fast
        Acc = zeros(1,105); %105 frames = 7 sec (2 before, 5 after)
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
            St = tAct(iAct) - 30;
            End = St + 104;
            tmp = subPix(:, St:End);
            tmp = mean(tmp- mean(tmp(:,1:10),2), 1); %2 sec before act.
            Acc = Acc + tmp;
        end
        Acc = Acc./iAct;
        mapAcc(indI,:) = Acc;
        [A, T] = max(Acc(30:end)); %onset to 5 sec after
        mapAmp(indI) = A;
        mapTime(indI) = T/15;
        
    end
    
    
    
    
    
    
    % map = zeros(512, 512);
    % map(mask) = mapAct;
    % mapAct = map;
    % map(mask) = mapTime;
    % mapTime = map;
    % map(mask) = mapAmp;
    % mapAmp = map;
    
    %in MapAcc is all the times a pixel went over the ca threshold on the first
    %dimension. The second dimension is the 7 seconds (2 before, 5 after) of
    %that event. This is the activity of hbo.
    plot(mean(mapAcc, 1))
    clear mapAct mapAmp mapTime ROI_info tmp fid idx ind indI img_info
    
    
    %% fluo
    fid = fopen([DataFolder 'hemoCorr_fluo.dat']);
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
    
    clear zF
    
    %% Segment the brain area
    load('/media/mbakker/GDrive/P2/GCaMP/M15/ROImasks_data.mat')
    mask = img_info.logical_mask;
    
    %% Reducing memory usage
    dF = reshape(dF,[], size(dF,3));
    aF = reshape(aF,[], size(aF,3));
    dF = dF(mask(:),:);
    aF = aF(mask(:),:);
    
    %% fluo
    aF(:,1:151) = 0;
    aF(:,(end-150):end) = 0;
    
    mapAct = zeros(size(aF,1),1, 'single');
    for indI = 1:size(aF,1)
        tmp = find(aF(indI, :));
        act{indI} = tmp;
        mapAct(indI) = length(tmp); %hoe vaak gaat pixel over threshold?
    end
    
    idx = find(mask);
    mapAmp = zeros(size(aF,1),1, 'single');
    mapTime = zeros(size(aF,1),1, 'single');
    mapAcc = zeros(size(idx,1),105, 'single');
    clear aF
    
    parfor( indI = 1:size(dF,1), 4) %5 is nr cores
        
        indexes = idx;
        Pixels = dF; %to make parfor fast
        Acc = zeros(1,105); %105 frames = 7 sec (2 before, 5 after)
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
            St = tAct(iAct) - 30;
            End = St + 104;
            tmp = subPix(:, St:End);
            tmp = mean(tmp- mean(tmp(:,1:10),2), 1); %2 sec before act.
            Acc = Acc + tmp;
        end
        Acc = Acc./iAct;
        mapAcc(indI,:) = Acc;
        [A, T] = max(Acc(30:end)); %onset to 5 sec after
        mapAmp(indI) = A;
        mapTime(indI) = T/15;
        
    end
    
    % map = zeros(512, 512);
    % map(mask) = mapAct;
    % mapAct = map;
    % map(mask) = mapTime;
    % mapTime = map;
    % map(mask) = mapAmp;
    % mapAmp = map;
    
    %in MapAcc is all the times a pixel went over the ca threshold on the first
    %dimension. The second dimension is the 7 seconds (2 before, 5 after) of
    %that event. This is the activity of hbo.
    plot(mean(mapAcc, 1))
    clear mapAct mapAmp mapTime ROI_info tmp fid idx ind indI img_info
    % save(
    mapAccHbO = load([DataFolder 'mapAccHbO.mat']);
    
    scatter(mapAccfluo, mapAccHbO);
    %
    % clear dH
    % fid = fopen([DataFolder 'hemoCorr_fluo.dat']);
    % dF = fread(fid, inf, '*single');
    % dF = reshape(dF, 512,512,[]);
    % fclose(fid)
    %
    % zF = (dF - mean(dF, 3))./std(dF,0,3);
    % aF = zF >= 1.95;
    % for ind = 1:size(aF,3)
    % aF(:,:,ind) = bwmorph(bwmorph(aF(:,:,ind), 'close', inf),'open',inf);
    % end
    % aF = aF(:,:,2:end)&~aF(:,:,1:(end-1));
    % aF = cat(3, false(size(aF,1),size(aF,2)), aF);
    %
    % dF = reshape(dF,[], size(dF,3));
    % dF = dF(mask(:),:);
    %
    % aF(:,1:151) = 0;
    % aF(:,(end-150):end) = 0;
    %
    % idx = find(mask);
    % mapAmpF = zeros(size(aF,1),1, 'single'); % 512 * 1
    % mapTimeF = zeros(size(aF,1),1, 'single'); % 512 * 1
    % mapAccF = zeros(size(idx,1),105, 'single'); %158650 * 105
    % clear aF
    %
    % parfor( indI = 1:size(dF,1), 4) %5 is nr cores
    %
    %     indexes = idx;
    %     Pixels = dF; %to make parfor fast
    %     Acc = zeros(1,105); %105 frames = 7 sec (2 before, 5 after)
    %     tAct = act{indI};
    %     [iY, iX] = ind2sub([512, 512], indexes(indI));
    %     subsY = iY + (-5:5); %neighbour pixels
    %     subsY = repelem(subsY,11);
    %     subsX = iX + (-5:5);
    %     subsX = repmat(subsX,1, 11);
    %     indices = sub2ind([512, 512], subsY, subsX);
    %     indices = find(ismember(indexes, indices));
    %     subPix = Pixels(indices,:);
    %     for iAct = 1:length(tAct)
    %         St = tAct(iAct) - 30;
    %         End = St + 104;
    %         tmp = subPix(:, St:End);
    %         tmp = mean(tmp- mean(tmp(:,1:10),2), 1); %2 sec before act.
    %         Acc = Acc + tmp;
    %     end
    %     Acc = Acc./iAct;
    %     mapAccF(indI,:) = Acc;
    %     [A, T] = max(Acc(30:end)); %onset to 5 sec after
    %     mapAmpF(indI) = A;
    %     mapTimeF(indI) = T/15;
    %
    % end
    
    Time = linspace(-2,5,105);
    
    
    
    
end

