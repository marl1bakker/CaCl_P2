% attention: only works for gcamp, not for speckle

function GetSeeds(DataFolder, DataFolderA1, overwrite)

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end
if( exist('DataFolderA1', 'var') && ~strcmp(DataFolderA1(end), filesep) )
    DataFolderA1 = [DataFolderA1 filesep];
end

if ~exist('overwrite', 'var')
    overwrite = 0;
end

idx = strfind(DataFolder, filesep); 
pathROI = DataFolder(1:idx(end-2));
if ~exist([pathROI 'BigROI.mat'], 'file')
    disp(['BigROI does not exist for ' DataFolder])
    return
end

if exist([DataFolder 'Seeds.mat'], 'file') && overwrite == 1
    disp(['OVERWRITING SEEDS' DataFolder])
elseif exist([DataFolder 'Seeds.mat'], 'file') && overwrite == 0
    disp('Seeds already made')
    return
end

% Mouse = DataFolder(idx(end-3)+1:idx(end-2)-1);
Acq = DataFolder(idx(end-2)+1:idx(end-2)+2);

%% Get data
dims = [512, 512];
fid = fopen([DataFolder 'fluo_567.dat']);
fl = fread(fid, dims(1)*dims(2), '*single');
fl = reshape(fl, dims);
fclose(fid);
fid = fopen([DataFolder 'green.dat']);
gr = fread(fid, dims(1)*dims(2), '*single');
gr = reshape(gr, dims);
fclose(fid);
% fid = fopen([DataFolder 'red.dat']);
% red = fread(fid, dims(1)*dims(2), '*single');
% red = reshape(red, dims);
% fclose(fid);

%% Get roi
load([pathROI 'BigROI.mat'], 'AtlasMask', 'BigROI', 'regions');

% shift atlas if it's another acquisition
if ~matches(Acq, 'A1')
    load([DataFolder 'tform_acq.mat'], 'tform')
    invtform = invert(tform);
    Rois = fieldnames(BigROI);
    for indroi = 1:size(Rois, 1)
        eval(['BigROI.' Rois{indroi} '= imwarp(BigROI.' Rois{indroi} ', invtform, ''OutputView'', imref2d(size(AtlasMask)), ''interp'', ''nearest'');'])
    end
    AtlasMask = imwarp(AtlasMask, invtform, 'OutputView', imref2d(size(AtlasMask)), 'interp', 'nearest');
end

GenMask = logical(AtlasMask);
MaskAllCentroids = zeros(dims);

%% if seeds were made for A1, get those as starting out point
if ~matches(Acq, 'A1') && exist('DataFolderA1', 'var') && ...
        exist([DataFolderA1 'Seeds.mat'], 'file')

    load([DataFolderA1 'Seeds.mat'], 'Mask', 'XY', 'adjusted')

    for indroi = 1:size(Rois, 1)
        eval(['regionmask = imwarp(Mask.' Rois{indroi} ', invtform, ''OutputView'', imref2d(size(AtlasMask)), ''interp'', ''nearest'');'])
        [r, c] = find(regionmask);
        if ~isempty(r) && ~isempty(c)
            r = r(round(size(r,1)/2)); %get middle pixel
            c = c(round(size(c,1)/2));
            eval(['XY.' Rois{indroi} ' = [c r];'])

            regionmask = zeros(dims);
            regionmask(r, c) = 1;
            regionmask = conv2(regionmask, ones(3,3), 'same');
            MaskAllCentroids = MaskAllCentroids + regionmask;
        end

        eval(['Mask.' Rois{indroi} '= regionmask;'])
    end
else
    %% get centroids
    for ind = 1:size(regions, 2)
        % disp(regions{ind})

        % Get middle of ROIs
        % Get centroid of ROI based on weight
        [X, Y] = meshgrid(1:dims(1), 1:dims(2));
        iX = sum(reshape(X.*BigROI.(regions{ind}), [], 1))/sum(BigROI.(regions{ind})(:));
        iY = sum(reshape(Y.*BigROI.(regions{ind}), [], 1))/sum(BigROI.(regions{ind})(:));
        iX = round(iX);
        iY = round(iY);

        eval(['XY.' regions{ind} ' = [iX iY];'])
        regionmask = zeros(dims);

        if ~isnan(iX) && ~isnan(iY)
            regionmask(iY, iX) = 1;
            regionmask = conv2(regionmask, ones(3,3), 'same');
            % eval(['Mask.' regions{ind} ' = regionmask;'])
            MaskAllCentroids = MaskAllCentroids + regionmask;
        end

        eval(['Mask.' regions{ind} ' = regionmask;'])
    end
end

%% check
answer = 'No';

while matches(answer, 'No')
    f1 = figure;
    imagesc(fl)
    f1.Position = [10 10 600 500];
    f2 = figure;
    imshowpair(fl, MaskAllCentroids)
    f2.Position = [610 10 600 500];
    f3 = figure;
    imshowpair(gr, MaskAllCentroids)
    axis on
    f3.Position = [1220 10 600 500];
    f4 = figure;
    imshowpair(fl, AtlasMask)
    f4.Position = [10 800 600 500];

    pause
    answer = questdlg('Are seeds correctly placed?');
    switch answer
        case 'Yes'
            close(f1, f2, f3, f4)
            adjusted = 0;
        case 'No'
            % give new coordinates
            coords = cell2mat(struct2cell(XY));
            Xs_old = cellstr(num2str(coords(:,1)));
            Ys_old = cellstr(num2str(coords(:,2)));

            fieldsize = repmat([1 20], size(fields(XY)));
            Xs_new = inputdlg(fields(XY), 'Give new X coords', fieldsize, Xs_old);
            Ys_new = inputdlg(fields(XY), 'Give new Y coords', fieldsize, Ys_old);

            close(f1, f2, f3, f4)
            MaskAllCentroids = zeros(dims);

            % get new mask to check
            for ind = 1:size(regions, 2)
                regionmask = zeros(dims);
                iX = str2double(Xs_new(ind));
                iY = str2double(Ys_new(ind));
                eval(['XY.' regions{ind} ' = [iX iY];'])

                if ~isnan(iX) && ~isnan(iY)
                    regionmask(iY, iX) = 1;
                    regionmask = conv2(regionmask, ones(3,3), 'same');
                    eval(['Mask.' regions{ind} ' = regionmask;'])
                    MaskAllCentroids = MaskAllCentroids + regionmask;
                end
            end
            adjusted = [~matches(Xs_old, Xs_new) ~matches(Ys_old, Ys_new)]; 

        case 'Cancel'
            error('Get seeds operation canceled')
    end
end

%% save
save([DataFolder 'Seeds.mat'], 'Mask', 'XY', 'adjusted')

end
