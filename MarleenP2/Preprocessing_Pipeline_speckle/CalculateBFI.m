% Calculates BFI. Does not take into account the movement of the mouse yet.

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
sdat = movstd(dat, 10, 0, 3, 'omitnan'); % over time
% mean(sdat(:,:,23), 'all', 'omitnan'); %330.92
% sdat = stdfilt(dat); %over space?

mdat = movmean(dat, 10,3, 'omitnan');
% mean(mdat(:,:,23), 'all', 'omitnan'); %9.8e+03
% kernel = ones(3,3);
% mdat = conv2(dat(:,:,1), kernel, 'same')./conv2(ones(1024,1024), kernel, 'same');
clear dat

K = sdat./mdat;

clear sdat mdat
BFI = 1./(K.^2);
clear K dat

% Save
fid = fopen([DataFolder 'BFI.dat'], 'w');
% fid = fopen([DataFolder 'speckle.dat'], 'w');
fwrite(fid,BFI,'single');
fclose(fid);

meanBFI = mean(BFI, 3);
save([DataFolder 'meanBFI.mat'], 'meanBFI');

delete([DataFolder 'speckle.dat']);

end

