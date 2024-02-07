% Calculates BFI. Does not take into account the movement of the mouse yet.
% 

function CalculateBFI(DataFolder)

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

% Get data
fid = fopen([DataFolder 'speckle.dat']);
dat = fread(fid, inf, '*single');
fclose(fid);
dat = reshape(dat, 1024, 1024, []);

% Caluclate Blood Flow Index
sdat = movstd(dat, 10, 0, 3, 'omitnan');
mean(sdat(:,:,23), 'all', 'omitnan'); %330.92

mdat = movmean(dat, 10,3, 'omitnan');
mean(mdat(:,:,23), 'all', 'omitnan'); %9.8e+03

K = sdat./mdat;
mean(K(:,:,23),'all','omitnan'); %0.0488

clear sdat mdat
BFI = 1./(K.^2);
% imagesc(mean(BFI,3,'omitnan'))
clear K dat

% Save
fid = fopen([DataFolder 'BFI.dat'], 'w');
fwrite(fid,BFI,'single');
fclose(fid);

figure
imagesc(mean(BFI, 3),[1 1500])
colormap(jet)

end