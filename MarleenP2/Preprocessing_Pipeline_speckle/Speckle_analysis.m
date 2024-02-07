fid = fopen('speckle.dat');
dat = fread(fid, inf, '*single');
fclose(fid);
dat = reshape(dat, 1024, 1024, []);

dat = dat./mean(dat,3);

%% temporal speckle
dat = reshape(dat, 1024*1024, []);
kern = ones(1,5,'single');
parfor (indP = 1:1024*1024, 8)
    dat(indP,:) = stdfilt(dat(indP,:),kern);
end
    
dat = reshape(dat, 1024, 1024, []);
imagesc(mean(dat,3)); % for nice vessel view

bla = dat./mean(dat,3);

%% spatial speckle
dat = reshape(dat, 1024,1024, []);
kern = ones(5,'single');

for indF = 1:size(dat, 3)
    dat(:,:,indF) = stdfilt(dat(:,:,indF),kern);
end

bla = movmean(dat,5,3);
bla = dat./mean(dat,3);
bla = imgaussfilt(bla,2.5);

figure
for ind = 1:500
    imagesc(bla(:,:,ind), [0.5 3])
    title(num2str(ind))
    pause(0.05)
end

figure
for ind = 1:500
    imagesc(-log(bla(:,:,ind)), [2 5])
    title(num2str(ind))
    pause(0.05)
end

bla2 = -log(bla);
bla2 = bla2./mean(bla2,3);

figure
for ind = 1:500
    imagesc(bla2(:,:,ind), [0.9 1.1])
    title(num2str(ind))
    pause(0.05)
end




