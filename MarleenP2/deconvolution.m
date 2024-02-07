%% Deconvolution
% for proof of concept, see after function

% should we deconvolve our fluo response first by the GCaMP response curve?

% method can be: 'singlepixel', 'wholebrain', 'allpixels', 'wiener'

function deconvolution(DataFolder, method)

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

% load data
fid = fopen([DataFolder 'hemoCorr_fluo.dat']);
df = fread(fid, inf, '*single');
df = reshape(df, 512, 512, []);

if exist([DataFolder 'HbONoFilt.dat'], 'file')
    fid = fopen([DataFolder 'HbONoFilt.dat']);
    dh = fread(fid, inf, '*single');
    dh = reshape(dh, 512, 512, []);
else
    [dh, ~] = HemoComputeNoFilter(DataFolder, DataFolder, 'GCaMP', {'red','green'}, 1);
    dh = reshape(dh, 512, 512, []);
end

map = df(:,:,1);

%make sure you dont include frames with movement
load([DataFolder 'MovMask.mat']);
MovMask = ~(conv(MovMask<1, [zeros(1,50) ones(1,100)],'same')>0);
df = df(:,:,MovMask>0);
dh = dh(:,:,MovMask>0);
% The moved frames are marked but we expand them a bit so we have all the
% movement. We take the frames with movement out completely. This shouldnt
% mess with the fourier transform. Even though there's a cut in the signal,
% it's still the same frequency.

%load the mask, make sure you dont have anything outside of the brain
fileseps = strfind(DataFolder, filesep);
load([DataFolder(1:fileseps(end-2)) 'ROImasks_data.mat']);
mask = img_info.logical_mask;
df = df .* mask;
dh = dh .* mask;

% Normalize
df = df./mean(df, 3);

if matches(method, 'singlepixel')
    %choose a single pixel
    imagesc(map)
    df = squeeze(mean(mean(df(350:360, 295:305, 1:7800), 1), 2)); %m15?
    dh = squeeze(mean(mean(dh(350:360, 295:305, 1:7800), 1), 2));
    %[160 330]
    % df = squeeze(mean(mean(df(325:335, 155:165, 1:end), 1), 2));
    % dh = squeeze(mean(mean(dh(325:335, 155:165, 1:end), 1), 2));
    % df = squeeze(mean(mean(df(335:340, 315:325, 1:end), 1), 2)); %m17
    % dh = squeeze(mean(mean(dh(335:340, 315:325, 1:end), 1), 2));
    
    if rem(size(df,1), 2) == 1 %to make sure you have an even nr, otherwise fft doesn't work
        df = df(1:end-1,:);
        dh = dh(1:end-1,:);
    end
    
    % transfer to fourier space
    df = fft(df-1);
    dh = fft(dh);
    ha = padarray(hann(2000), (size(df,1)-2000)/2, 0, 'both');
    % hann array is soort bell curve. We multiply it by this to reduce the
    % high frequencies, which is where most of the noise is.
    hr = ifftshift(fftshift(dh).*ha)./df; %ask sam
    % first shift dh to be centered, then multiply by the hann array and divide
    % by fluo signal.
    
    %plot
    plot(real(hr))
    title([DataFolder(fileseps(end-3)+1:fileseps(end-2)-1) ' whole brain'])
    
elseif matches(method, 'wholebrain')
    df(df == 0) = NaN;
    dh(dh == 0) = NaN;
    %take the whole brain
    df = squeeze(mean(mean(df,1, 'omitnan'),2, 'omitnan'));
    dh = squeeze(mean(mean(dh,1, 'omitnan'),2, 'omitnan'));
    
    if rem(size(df,1), 2) == 1
        df = df(1:end-1,:);
        dh = dh(1:end-1,:);
    end
    
    df = fft(df-1);
    dh = fft(dh);
    ha = padarray(hann(2000), (size(df,1)-2000)/2, 0, 'both');
    hr = ifftshift(fftshift(dh).*ha)./df;
    
    %plot
    figure
    plot(abs(hr))
    
    hr = ifft(hr);
    hr = fftshift(hr);
    hr = real(hr);
    
    figure
    plot(hr)
    
    title([DataFolder(fileseps(end-3)+1:fileseps(end-2)-1) ' whole brain'])
    
elseif matches(method, 'allpixels') %do fft on each pixel of the brain, but not averaged.
    if rem(size(df,3), 2) == 1
        df = df(:,:,1:end-1);
        dh = dh(:,:,1:end-1);
    end
    
    % transfer to fourier space
    df = fft(df-1, [], 3);
    dh = fft(dh, [], 3);
    ha = padarray(hann(3000), (size(df,3)-3000)/2, 0, 'both');
    hr = ifftshift(bsxfun(@times,fftshift(dh,3),permute(ha, [3 2 1])), 3);
    
    %hr = ifftshift(fftshift(hr).*ha);
    hr = ifft(hr, [], 3); %do inverse fourier transform
    hr = fftshift(hr, 3); %
    hr = real(hr);
    
    %     figure
    %     plot(squeeze(hr(300,300,:)))
    figure
    plot(squeeze(mean(mean(hr,2, 'omitnan'),1,'omitnan')))
    
    
    %     time = linspace(-500/15, 500/15, 1000);
    %     figure
    %     for ind = 3500:4500
    %         imagesc(hr(:,:,ind),[-5 5])
    %         title(num2str(time(ind-3499)))
    %         pause(0.02)
    %     end
    
    title([DataFolder(fileseps(end-3)+1:fileseps(end-2)-1) 'av all pixels'])
    
elseif matches(method, 'wiener')

    % dF = dF./mean(dF,2);
    nbpts = size(df,3);
    
    df = reshape(df,512*512, []);
    dh = reshape(dh, 512*512, []);
    dh = dh(mask(:),:); %make sure you dont have nan!
    df = df(mask(:),:);
    
    df = fft(df-1, 2^nextpow2(size(df,2)),2); %fft optimization
    df = fftshift(df,2);
    
    dh = fft(dh, 2^nextpow2(size(dh,2)),2);
    dh = fftshift(dh,2);

    freq = linspace(-7.5, 7.5, size(dh,2));
    figure
%     plot(freq, squeeze(abs(mean(dh,1, 'omitnan'))));
%     figure
%     plot(freq, squeeze(abs(mean(df,1,'omitnan'))));
    
    sigma = 0.1; %can play with this
    WC = (abs(df).^2)./(abs(df).^2 + 1/sigma);
    fH = (dh./df).*WC;
    clear df dh
    fH = ifftshift(fH,2);
    dH = ifft(fH,nbpts, 2); %give different name
    dH = fftshift(real(dH),2);
    
    newdH = NaN(512*512,size(dH,2),'single');
    for ind = 1:size(dH,2)
        newdH(mask(:),ind) = dH(:,ind);
    end
    newdH = reshape(newdH, 512,512,[]);
    
    %%
    
%     figure;
%     for ind= 600:900
%         imagesc(newdH(:,:,ind),[-3 3]);
%         pause(0.1);
%     end
    
%     figure
%     plot(mean(dH(:,(size(dH,2)/2-200):(size(dH,2)/2+200)),1))
%% M14 weird peak in graph!!!
figure
plot(squeeze(mean(mean(newdH,2,'omitnan'),1,'omitnan')))

    Titel = [DataFolder(fileseps(end-3)+1:fileseps(end-2)-1) '-wiener'];
    title(Titel)
    saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/Deconvolution/' Titel '.png'])
    
end



end
























%% explanation/proof
% a = zeros(1000,1);
% a(100:200:end) = 1;
% a(200:201) = 1;
% b = sin(0:0.1:pi);
%
% c = conv(a,b,'same');
% plot(c)
%
% %Zscore:
% zF = (c - mean(c, 1))./std(c,0,1);
% %Threshold on Zscore:
% aF = zF >= 1.95;
% %Now, we want only the beginning of activations:
% aF = aF(2:end)&~aF(1:(end-1));
% aF = cat(1, false, aF);
%
% Acc = zeros(60,1);
% tAct = find(aF>0);
% for iAct = 1:length(tAct)
%     St = tAct(iAct) - 30; %50 is for 10s at 5 Hz
%     En = St + 59;
%     tmp = c(St:En);
%     Acc = Acc + tmp;
% end
% Acc = Acc./iAct;
%
% %%
% % The idea is that a and b convoluted will give c. A multiplication in the
% % fourier space is the same as a convolution in the time space. Let's test
% % this. First we shift a into the fourier space. Since a and b are
% % different lenghts, we have to pad b with zeros around. Both a and b in
% % the fourier space give a bit of an odd plot, so we shift them to the
% % middle. Then, we multiply a and b in the fourier space. If we then shift
% % it back (ifftshift = inverse fft shift) and do an inverse of the fourier
% % transform to get it back to timespace (ifft), we get the same as if we
% % did a convolution of a and b! So d, as calculated in fourier space, is
% % the same as c, as calculated with a direct convolution.
% afft = (fft(a));
% b = padarray(b, [0 (size(a,1)-size(b,2))/2], 0, 'both');
% bfft = (fft(b));
%
% dfft = afft.*bfft';
% d = fftshift(ifft(dfft));
%
% cfft = fft(c);
%
% b2 = cfft./afft;
% b2 = ifft(b2);
% b2 = fftshift(b2);
%
%

