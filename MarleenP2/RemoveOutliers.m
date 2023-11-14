function RemoveOutliers(DataFolder, dataname)

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

if ~exist('dataname', 'var')
    dataname = 'hemoCorr_fluo';
end

if ~exist([DataFolder dataname '.dat'], 'file')
    disp([dataname ' could not be found, function exited'])
    return
end

%% get data
fid = fopen([DataFolder dataname '.dat']);
dat = fread(fid, inf, '*single');
dat = reshape(dat,512,512,[]);
fclose(fid);
dims = size(dat);

% imagesc(dat(:,:,23))






end