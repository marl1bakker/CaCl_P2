fid = fopen('normLPF.dat');
% fid = fopen('HbO.dat');
dat = fread(fid, inf, '*single');
dat = reshape(dat, 512, 512, []);

for ind = 1:100
imagesc(dat(:,:,ind), [0.95 1.05])
colorbar
% imagesc(dat(:,:,ind), [-10 10])
pause(0.05)
end

v = VideoWriter('M15-A1-R2.avi', 'Uncompressed AVI');
v.FrameRate = 30;
open(v)

% for ind = 1:100
for ind = 1:500
    imagesc(dat(:,:,ind), [0.95 1.05])
%     imagesc(dat(:,:,ind), [-10 10])
    axis image
    drawnow
    F(ind) = getframe(gcf);
end
writeVideo(v,F)

close(v)

% Movement, M8
DataFolder =  '/media/mbakker/GDrive/P2/GCaMP/M8/A1-R1/CtxImg/';
dataname = 'hemoCorr_fluo';
fid = fopen([DataFolder dataname '.dat']);
dat = fread(fid, inf, '*single');
dat = reshape(dat,512,512,[]);
fclose(fid);
dims = size(dat);

cd '/media/mbakker/GDrive/P2/GCaMP/Videos'

v = VideoWriter('M8-A1-R1_mov.avi', 'Uncompressed AVI');
v.FrameRate = 15;
open(v)

f = figure;
f.Position = [10 10 500 800];
for ind = 200:1000
    tiledlayout(2,1)
    nexttile
    imagesc(dat(:,:,ind), [0.95 1.05])
    colorbar
    title(num2str(ind))
    nexttile
    imagesc(MovMask(ind), [0 1])
%     pause(0.05)
    drawnow
    F(ind) = getframe(gcf);
end
writeVideo(v,F)
close(v)

f = figure;
f.Position = [10 10 500 800];
for ind = 200:1000
    tiledlayout(2,1)
    nexttile
    imagesc(dat(:,:,ind), [0.95 1.05])
    title(num2str(ind))
    nexttile
    imagesc(MovMask(ind), [0 1])
    pause(0.0667)
end