% to plot periphery of seed spread

% DataFolder = pwd;
DataFolder = '/media/mbakker/GDrive2/P2/GCaMP/M15/A1-R2/CtxImg';

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

dataname = 'hemoCorr_fluo.dat';

seps = strfind(DataFolder, filesep);
load([DataFolder(1:seps(end-2)) 'BigROI.mat'], 'AtlasMask', 'regions');
Mouse = DataFolder(seps(end-3)+1:seps(end-2)-1);

% Get data
fid = fopen([DataFolder dataname]);
dat = fread(fid, 512*512, '*single');
fclose(fid);
dat = reshape(dat, 512, 512);
dims = size(dat);

% %% pics
% for ind = 1:size(regions,2) %go per region
%     if ind == 3 || ind == 8 %if it's the auditory cortex, skip
%         continue
%     end
%     
%     Mask = ismember(AtlasMask,ind); %pak nummers van atlas van ind waar je nu bent
% 
%     % Get centroid of ROI based on weight
%     [X, Y] = meshgrid(1:dims(1), 1:dims(2));
%     iX = sum(reshape(X.*Mask, [], 1))/sum(Mask(:));
%     iY = sum(reshape(Y.*Mask, [], 1))/sum(Mask(:));
%     iX = round(iX);
%     iY = round(iY);
%     
%     Seed = zeros(dims(1), dims(2));
%     Seed(iY, iX) = 1;
%     radius = 3;
%     SE = strel('sphere', radius);
%     Seed = imdilate(Seed, SE);
%     
%     %     Orig_seed = Seed; %to have later?
%     Seed = logical(Seed);
%     %     Orig_Seed = Seed;
%     %     timecourse_seed = mean(dat(Seed(:), :),1);
%     %     timecourse_seed = dat(iY, iX, :);
%     
%     periphery = 12;
%     per_corrs = [];
%     iteration = 1; 
%     %     while any(periphery, 'all') && radius < 60
%     while sum(periphery, 'all') > 11 && radius < 50
%         
%         %take growing circle around seed
%         SE = strel('sphere', radius);
%         radius = radius + 5; % for next round
%         periphery = imdilate(Seed, SE);
%         periphery = periphery - Seed;
%         periphery = periphery.* Mask;
%         
%         % plot
%         figure
%         imshowpair(dat, periphery)
%         hold on
%         plot(iX, iY, '.', 'color', 'red', 'MarkerSize', 5, 'LineWidth', 3)
%         
%         saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/SeedSpread/ImagePeriphery/'...
%             regions{ind} '_' num2str(iteration) '.tiff']);
%         
%         close all
%         
%         Seed = periphery + Seed;
%         iteration = iteration + 1;
%     end
% end

%% video

v = VideoWriter(['/media/mbakker/GDrive/P2/GCaMP/SeedSpread/ImagePeriphery/' Mouse ...
    'all-periphery.avi'], 'Uncompressed AVI');
v.FrameRate = 2;
open(v)

framenumber = 0;
for ind = 1:size(regions,2) %go per region
    if ind == 3 || ind == 8 %if it's the auditory cortex, skip
        continue
    end
    
%     v = VideoWriter(['/media/mbakker/GDrive/P2/GCaMP/SeedSpread/ImagePeriphery/' Mouse ...
%         regions{ind} '-periphery.avi'], 'Uncompressed AVI');
%     v.FrameRate = 2;
%     open(v)

    Mask = ismember(AtlasMask,ind); %pak nummers van atlas van ind waar je nu bent

    % Get centroid of ROI based on weight
    [X, Y] = meshgrid(1:dims(1), 1:dims(2));
    iX = sum(reshape(X.*Mask, [], 1))/sum(Mask(:));
    iY = sum(reshape(Y.*Mask, [], 1))/sum(Mask(:));
    iX = round(iX);
    iY = round(iY);
    
    Seed = zeros(dims(1), dims(2));
    Seed(iY, iX) = 1;
    radius = 3;
    SE = strel('sphere', radius);
    Seed = imdilate(Seed, SE);
    
    %     Orig_seed = Seed; %to have later?
    Seed = logical(Seed);
    %     Orig_Seed = Seed;
    %     timecourse_seed = mean(dat(Seed(:), :),1);
    %     timecourse_seed = dat(iY, iX, :);
    
    periphery = 12;
    per_corrs = [];
    iteration = 1; 
    
    radius = 10;
    
    %     while any(periphery, 'all') && radius < 60
    while sum(periphery, 'all') > 11 && radius < 50
        
        %take growing circle around seed
        SE = strel('sphere', radius);
%         radius = radius + 5; % for next round % oud
        periphery = imdilate(Seed, SE);
        periphery = periphery - Seed;
        periphery = periphery.* Mask;
        
        % plot
        figure
        imshowpair(dat, periphery)
        hold on
        plot(iX, iY, '.', 'color', 'red', 'MarkerSize', 5, 'LineWidth', 3)
        axis equal
        pause(0.05)
        
        drawnow
        F(framenumber + iteration) = getframe(gcf);
%         saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/SeedSpread/ImagePeriphery/'...
%             regions{ind} '_' num2str(iteration) '.tiff']);
        if ~isequal(size(F(framenumber + iteration).cdata), [600 692 3])
            disp('s')
        end
        
        close all
        
        Seed = periphery + Seed;
        iteration = iteration + 1;
        
    end
    
    framenumber = framenumber + iteration - 1;
%     writeVideo(v,F)

%     writeVideo(v,F)
%     close(v)
%     clear v F
end

writeVideo(v,F)
close(v)
clear v F

%% image

for ind = 1:size(regions,2) %go per region
    if ind == 3 || ind == 8 %if it's the auditory cortex, skip
        continue
    end

    Mask = ismember(AtlasMask,ind); %pak nummers van atlas van ind waar je nu bent

    % Get centroid of ROI based on weight
    [X, Y] = meshgrid(1:dims(1), 1:dims(2));
    iX = sum(reshape(X.*Mask, [], 1))/sum(Mask(:));
    iY = sum(reshape(Y.*Mask, [], 1))/sum(Mask(:));
    iX = round(iX);
    iY = round(iY);
    
    Seed = zeros(dims(1), dims(2));
    Seed(iY, iX) = 1;
%     radius = 3; % oud
    radius = 10; %new
    SE = strel('sphere', radius);
    Seed = imdilate(Seed, SE);
    
    Seed = logical(Seed);

    periphery = 12;
    per_corrs = [];
    iteration = 1; 
    %     while any(periphery, 'all') && radius < 60
    while sum(periphery, 'all') > 11 && radius < 50
        
        %take growing circle around seed
        SE = strel('sphere', radius);
%         radius = radius + 5; % for next round, oud
        periphery = imdilate(Seed, SE);
        periphery = periphery - Seed;
        periphery = periphery.* Mask;
        
        eval(['periphery' num2str(iteration) ' = periphery.*iteration;']);
        
        Seed = periphery + Seed;
        iteration = iteration + 1;
        
    end
    
    imageseedspread = Seed + periphery1 +periphery2+periphery3+periphery4+periphery5+periphery6;
    imagesc(imageseedspread)
    saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/SeedSpread/ImagePeriphery/'...
        regions{ind} '.tiff']);
    
    close all
    eval(['imageseedspread' num2str(ind) ' = Seed + periphery1 +periphery2+periphery3+periphery4+periphery5+periphery6;'])

end

