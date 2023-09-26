% choose one pixel in the brain, and make a correlationmap of that. 
% ManualInput 0 - dont ignore existing files
% ManualInput 1 - overwrite existing files

%% P2
% Input Seedname can be several, give like {'seed1', 'seed2'} etc.
%{'VisualROI_R'} SensoryROI_R AuditoryROI_R UnknownROI_R MotorROI_R

%goes with BigROI standard, because the small rois will give way too many
%points. 

function SingleSubjectSPCM(DataFolder, dataname, Overwrite, GSR)

if ~exist('Overwrite', 'var')
    Overwrite = 0;
end

if ~exist('GSR', 'var')
    GSR = 0;
end

if ~exist('dataname', 'var')
    dataname = 'mov_aligned.dat';
end

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

% if it already exists, dont do it
if( exist([DataFolder 'xxx'], 'file') && Overwrite == 0 )
    disp('Seed pixel correlation map already done, function exited')
    return
elseif( exist([DataFolder 'xxx'], 'file') && Overwrite == 1 )
    disp('Seed pixel correlation map already done, OVERWRITING FILES')
end 

%check if you have the clustered ROI, load it
idx = strfind(DataFolder, filesep); %zoek alle plekken van fileseps in naam
pathROI = DataFolder(1:idx(end-2));

if ~exist([pathROI 'BigROI.mat'], 'file')
    disp(['BigROI does not exist for ' DataFolder])
    return
end

Mouse = DataFolder(idx(end-3)+1:idx(end-2)-1);
Acq = DataFolder(idx(end-2)+1:idx(end-1)-1);

clear idx

%% Get data
fid = fopen([DataFolder dataname]);
dat = fread(fid, inf, '*single');
dat = reshape(dat, 512,512, []);

load([pathROI 'BigROI.mat'], 'AtlasMask', 'BigROI', 'regions');
GenMask = logical(AtlasMask);
dat = dat.* GenMask;
dat(dat == 0) = NaN;

%% GSR
dims = size(dat);

if GSR == 1
    dat = reshape(dat,[], dims(3));
    mS = mean(dat,1, 'omitnan');
    
    X = [ones(size(mS)); mS];
    B = X'\dat';
    A = (X'*B)';
    dat = dat./A;
    dat = reshape(dat,dims);
    clear h Mask mS X B A;
end

%% Do everything per seedname you gave in
for ind = 1:size(regions, 2) 
    disp(regions{ind})
    
    %% Get middle of ROIs
    % Get centroid of ROI based on weight
    [X, Y] = meshgrid(1:dims(1), 1:dims(2));
    iX = sum(reshape(X.*BigROI.(regions{ind}), [], 1))/sum(BigROI.(regions{ind})(:));
    iY = sum(reshape(Y.*BigROI.(regions{ind}), [], 1))/sum(BigROI.(regions{ind})(:));
    iX = round(iX);
    iY = round(iY);
    
    
    %% Calculate the seed pixel correlation map
    dat = reshape(dat,dims);
    Seeddat = dat(iY, iX, :);
    Seeddat = reshape(Seeddat, 1, []);
    dat = reshape(dat, dims(1)*dims(2), []);
    [rho, pval] = corr(Seeddat', dat(:,:)');
    rho = reshape(rho, dims(1), dims(2));
    
    figure()
    if GSR == 0
        imagesc(rho, [0 1]) %without gsr
    else
        imagesc(rho, [-1 1])
    end
%     colormap jet
    load('/media/mbakker/data1/Hypoxia/SeedPixelCorrMap/NL.mat');
    colormap(NL)
    colorbar
    
    axis image
%     hold on
%     line([5 20.7], [5 5]);
%     line([5 5], [5 20.7]); %for scale, it's 157 pix for 10 mm
    line([5 55], [5 5], 'color', 'yellow'); %50 pix per mm
    line([5 5], [5 55], 'color', 'yellow'); 
    
%     idx = strfind(DataFolder, filesep); 
%     MouseAcq = [DataFolder(idx(end-2):end) 'Frames_' num2str(StartFrame) '_to_' num2str(EndFrame) Seedname{ind}]; 
%     MouseAcq = replace(MouseAcq, filesep, '_');
%     if GSR == 1
%         saveas(gcf, ['/media/mbakker/data1/Hypoxia/SeedPixelCorrMap/fluo/fluo' MouseAcq '.tiff']);
%     else
%         saveas(gcf, ['/media/mbakker/data1/Hypoxia/SeedPixelCorrMap/fluo/fluo' MouseAcq '_NoGSR.tiff']);
%     end

    if GSR == 1
        saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/SPCM/Individual/' Mouse '-' Acq ...
            '_' regions{ind} '_GSR.tiff'], 'tiff');
        saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/SPCM/Individual/' Mouse '-' Acq ...
            '_' regions{ind} '_GSR.eps'], 'epsc');
    else
        saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/SPCM/Individual/' Mouse '-' Acq ...
            '_' regions{ind} '_noGSR.tiff'], 'tiff');
        saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/SPCM/Individual/' Mouse '-' Acq ...
            '_' regions{ind} '_noGSR.eps'], 'epsc');
    end
    close gcf
    
end
end