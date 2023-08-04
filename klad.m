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

