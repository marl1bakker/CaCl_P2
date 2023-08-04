fid = fopen('normLPF.dat');
dat = fread(fid, inf, '*single');
dat = reshape(dat, 512, 512, []);

for ind = 1:1000
imagesc(dat(:,:,ind), [0.95 1.05])
pause(0.05)
end

v = VideoWriter('M15-A1-R2.avi', 'Uncompressed AVI');
v.FrameRate = 30;
open(v)

for ind = 1:100
    imagesc(dat(:,:,ind), [0.95 1.05])
    axis image
    drawnow
    F(ind) = getframe(gcf);
end
writeVideo(v,F)

close(v)
