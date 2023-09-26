%% periodogram
% fid = fopen('normLPF.dat');
fid = fopen('green.dat');
dat = fread(fid, inf, '*single');
dat = reshape(dat, 512, 512, []);

% On imagesc if you find a pixel [X,Y], put it as dat(y,x,:)
pix = dat(281, 208,:);
pix = reshape(pix, 1,[]);

pxx = periodogram(pix);
periodogram(pix)
[p, f, t] = pspectrum(pix, 15, 'spectrogram', 'OverlapPercent', 99);

[s,f,t] = spectrogram(pix,hann(300),299,256,15,"power","yaxis");

%% pca
% Principal component analysis.
% [C, S] = pca(dat)
% In matlab, waarbij dat xas * yas, [] is.
% C geeft de componenten. De eerste component gaat verreweg het meeste verklaren. 
% Daarna loopt het af. Als je wilt zien hoeveel je componenten verklaren kun 
% je [C,S, ~,~, E] = pca(dat) doen, en dan plot(cumsum(E), ‘.’). Je krijgt 
% dan te zien hoe veel procent van het signaal elke component verklaart.
%
% De componenten zijn timecourses. Er wordt in matlab het maximaal aantal 
% componenten gegeven, wat gelijk is aan het aantal tijdpunten. Hoe zwaar 
% een component meetelt voor een bepaalde pixel wordt aangegeven met S (score). 
% S is dus de hoeveelheid pixels keer de hoeveelheid componenten.
% 
% We willen groeperen op de score van de componenten. Als de scores van twee 
% pixels dicht bij elkaar liggen dan horen deze pixels waarschijnlijk bij 
% dezelfde groep. We groeperen met kmeans:
% 
% Idx = kmeans(S(:,1:10),16)
% 
% Op deze manier zou het betekenen dat je de eerste 10 componenten in 
% beschouwing neemt, en dat je vraagt om 16 clusters.

fid = fopen('fluo_567.dat');
dat = fread(fid, inf, '*single');
dat = reshape(dat, 512*512, []);

dat = dat./mean(dat,2); %normalize

dat = reshape(dat, 512, 512, []);
imagesc(dat(:,:,20))
h =  drawpolygon;
mask = h.createMask;
dat = reshape(dat, 512*512, []);

% movement remov 
load('MovMask.mat')
dat = dat.*MovMask;

%pca
[C,S,~,~,E] = pca(dat(mask(:),:)); %principal component analysis, C is coefficient 
plot(cumsum(E),'.')

figure
idx = kmeans(S(:,1:16),8); %get clusters 
map = zeros(512);
map(mask) = idx;
imagesc(map)
title('M30 PCA, 8 clusters, 16 components')

% register within acq

%% IVA

%% non negative matrix factorization
[W,H] = nnmf(dat(mask(:),:),20);

idx = kmeans(W(:,1:16),8,'Replicates',10);
% idx = kmedoids(W(:,1:16),8,'Replicates',10);

figure
map = zeros(512);
map(mask) = idx;
imagesc(map)
title('M26-A2 nnmf, 8 clusters, 16 components')




%% deconvolution
% normal deconvolution:
load('mapAccfluo.mat')
mapAccfluo = mapAcc;
load('mapAccHbO.mat')
mapAccHbO = mapAcc;
clear mapAcc

%single pixel:
figure
tiledlayout(2,1)
nexttile
plot(mapAccHbO(654,:))
nexttile
plot(mapAccfluo(654,:))
[q,r] = deconv(mapAccHbO(654,:), mapAccfluo(654,:));
figure
plot(r)
%is this something?

%average over all pixels:
avfluo = mean(mapAccfluo, 1);
avhbo = mean(mapAccHbO, 1);

[q,r] = deconv(avhbo, avfluo);
figure
plot(r)
% only get one value for q and a curve for r that's the same as fluo almost

%maybe the problem is that there are some values close to 0?
avfluo = avfluo+100;
avhbo = avhbo +100;
[q,r] = deconv(avhbo, avfluo);
% same problem

%maybe it's because of the noise?
%try wiener hehe
bla1 = deconvwnr(avhbo, avfluo);
bla2 = deconvwnr(double(mapAccHbO), double(mapAccfluo));
%???????


%%

figure;
for ind = 1:9003
    imagesc(reshape(dat(:,ind),512,512),[0.95 1.05]);
    axis image;
    pause(0.01);
end

