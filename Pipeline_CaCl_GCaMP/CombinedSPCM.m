% choose one pixel in the brain, and make a correlationmap of that.
% ManualInput 0 - dont ignore existing files
% ManualInput 1 - overwrite existing files

%% P2
% Input Seedname can be several, give like {'seed1', 'seed2'} etc.
%{'VisualROI_R'} SensoryROI_R AuditoryROI_R UnknownROI_R MotorROI_R

%goes with BigROI standard, because the small rois will give way too many
%points.

function CombinedSPCM(dataname, Acquisition, GSR, SaveDir, Overwrite)


if ~exist('SaveDir', 'var')
    SaveDir = '/media/mbakker/GDrive/P2/GCaMP';
end

if ~exist('Acquisition', 'var')
    Acquisition = 'A1';
end

if ~exist('Overwrite', 'var')
    Overwrite = 0;
end

if ~exist('GSR', 'var')
    GSR = 0;
end

if ~exist('dataname', 'var')
    dataname = 'hemoCorr_fluo';
end
dataname = [dataname '.dat'];

% if it already exists, dont do it
% if( exist([DataFolder 'xxx'], 'file') && Overwrite == 0 )
%     disp('Seed pixel correlation map already done, function exited')
%     return
% elseif( exist([DataFolder 'xxx'], 'file') && Overwrite == 1 )
%     disp('Seed pixel correlation map already done, OVERWRITING FILES')
% end

% Make groups for cacl/nacl
groups = {'CaCl', 'NaCl'};

load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat', 'RecordingOverview');
RecordingOverview.Group = categorical(RecordingOverview.Group);
indcacl = RecordingOverview.Group == 'CaCl';

load('/media/mbakker/GDrive/P2/GCaMP/BigROI_ReferenceMask.mat', 'AtlasMask');
BigROIRef = AtlasMask;
BigROIRefsmall = BigROIRef(135:345, 135:390);
clear AtlasMask;
% roi x muis x spcm (512*512);

for indgroup = 1:size(groups,2) % Go per group, nacl/cacl
    group = groups{indgroup};
    
    if matches(group, 'CaCl')
        Mousegroup = RecordingOverview(indcacl,:);
        disp('cacl group')
    elseif matches(group, 'NaCl')
        Mousegroup = RecordingOverview;
        Mousegroup(indcacl, :) = [];
        disp('nacl group')
    else
        disp('something is going wrong, group not recognized')
        return
    end
    
    AllSPCM = NaN(10, size(Mousegroup,1), 512*512);
    
    for indmouse = 1:size(Mousegroup, 1) %go per mouse
        Mouse = Mousegroup.Mouse{indmouse};
        disp(Mouse)
        eval(['DataFolder = [Mousegroup.' Acquisition '{indmouse} filesep];']);
        SaveFolder = [SaveDir filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];
        
        
        %check if you have the clustered ROI, load it
        indfileseps = strfind(SaveFolder, filesep); %zoek alle plekken van fileseps in naam
        pathROI = SaveFolder(1:indfileseps(end-2));
        
        if ~exist([pathROI 'BigROI.mat'], 'file')
            %             ALLSPCM(:, indmouse, :) = nan;
            disp(['No bigroi found for ' Mouse])
            continue
            %     disp(['BigROI does not exist for ' DataFolder])
        elseif ~exist([SaveFolder dataname], 'file')
            disp(['No ' dataname ' found for ' Mouse])
            continue
        end
            
        %% coregistration
        % Make sure the mice brains are in the same location by matching them to
        % your BigROIRef
        
        % Get ROI and mask
        load([pathROI 'BigROI.mat'], 'AtlasMask', 'BigROI', 'regions');
        GenMask = logical(AtlasMask);        
        
        AtlasMasksmall = AtlasMask(135:345, 135:390);
        
        %compute how much to shift
        [optimizer, metric] = imregconfig('monomodal');
        optimizer.MaximumIterations = 1e10;
        optimizer.MinimumStepLength = 1e-4;
        optimizer.MaximumStepLength = 1e-2;
        optimizer.GradientMagnitudeTolerance = 1e-3;
        tform = imregtform(AtlasMasksmall, BigROIRefsmall, 'similarity', optimizer, metric,...
            'DisplayOptimization', false, 'PyramidLevels', 3);
        
        % Confirmation
        AtlasMasktemp = imwarp(AtlasMask,tform,'OutputView',imref2d(size(BigROIRef)));
        ssimval = ssim(AtlasMasktemp(135:345, 135:390), BigROIRef(135:345, 135:390));
        
        if ssimval < 0.80 %if it seems bad, check manually
            f1 = figure;
            imshowpair(AtlasMasktemp, BigROIRef);
            title(Mouse)
            answer = questdlg('Does it make sense?', ...
                'Coregistration', ...
                'Yes','No','Yes');
            % Handle response
            switch answer
                case 'Yes'
%                     disp('≧◠‿●‿◠≦    ᕙ(^▿^-ᕙ)');
                    close(f1)
                case 'No'
                    disp(['Coreg. not applied, mouse ' Mouse ' skipped'])
                    for indregion = 1:size(regions, 2)
                        AllSPCM(indregion, indmouse, :) = NaN(512*512,1); %give nans for mouse
                    end
                    continue
            end
        end
        
        %% get data
        fid = fopen([SaveFolder dataname]);
        dat = fread(fid, inf, '*single');
        dat = reshape(dat, 512,512, []);
        fclose(fid);
        dat = dat.* GenMask;
        dat(dat == 0) = NaN;
        dims = size(dat);
        
        %% GSR
        if GSR == 1
            dat = reshape(dat,[], dims(3));
            mS = mean(dat,1, 'omitnan');
            
            X = [ones(size(mS)); mS];
            B = X'\dat';
            A = (X'*B)';
            dat = dat./A;
            dat = reshape(dat,dims);
            clear h Mask mS X B A;
        end        

        %% Do everything per seedname you gave in
        for indregion = 1:size(regions, 2)
            
            % Get centroid of ROI based on weight
            [X, Y] = meshgrid(1:dims(1), 1:dims(2));
            iX = sum(reshape(X.*BigROI.(regions{indregion}), [], 1))/sum(BigROI.(regions{indregion})(:));
            iY = sum(reshape(Y.*BigROI.(regions{indregion}), [], 1))/sum(BigROI.(regions{indregion})(:));
            iX = round(iX);
            iY = round(iY);
            
            % Calculate the seed pixel correlation map
            dat = reshape(dat,dims);
            if isnan(iY) || isnan(iX) %if the region falls outside of the mask
                rho = NaN(1, dims(1)*dims(2), 'single');
            
            else
                Seeddat = dat(iY, iX, :);
                Seeddat = reshape(Seeddat, 1, []);
                dat = reshape(dat, dims(1)*dims(2), []);
                [rho, ~] = corr(Seeddat', dat(:,:)');
                rho = reshape(rho, dims(1), dims(2));
                
                %coregister with tform that you calculated before
                rho = imwarp(rho,tform,'OutputView',imref2d([dims(1) dims(2)]));
                rho = reshape(rho, 1, dims(1)*dims(2));
                
                rho(rho == 0) = NaN;
            end
            
            % save it in the overview variable AllSPCM
            AllSPCM(indregion, indmouse, :) = rho;

        end %end indregion
    end %end indmouse
    
    clear dat GenMask metric rho tform Seeddat iX iY X Y %clean some stuff
    eval(['SPCM_' group ' = mean(AllSPCM, 2, ''omitnan'');']);
    
end %end indgroup
clear AllSPCM AtlasMask AtlasMasksmall AtlasMasktemp BigROIRefsmall fid group ...
    indcacl indfileseps ingroup indmouse indregion Mouse Mousegroup optimizer 

%% save spcm's per group
for indgroup = 1:size(groups,2) % Go per group, nacl/cacl
    group = groups{indgroup};
    
    for indregion = 1:size(regions, 2)
        f = figure;
        eval(['imagesc(reshape(SPCM_' group '(indregion,1,:), 512, 512), [0 1]);'])
        title([group ' ' Acquisition ' ' regions{indregion}], 'interpreter', 'none')
        axis off
        colorbar
        if GSR == 0
            f.CurrentAxes.CLim = [0 1];
            saveas(gcf, [SaveDir '/SPCM/Combined/withoutGSR/' dataname '_' group '_' Acquisition '_' regions{indregion} '.tiff'])
        else
            f.CurrentAxes.CLim = [-1 1];
            saveas(gcf, [SaveDir '/SPCM/Combined/withGSR/' dataname '_' group '_' Acquisition '_' regions{indregion} '.tiff'])
        end
            %         pause
        close(f)
    end %of regions
    
end %end of saving spcm's per group

%% Save SPCM, difference
for indregion = 1:size(regions,2)
    f = figure;
    SPCM_diff = SPCM_NaCl(indregion, 1, :) - SPCM_CaCl(indregion, 1, :);
    SPCM_diff = reshape(SPCM_diff, 512, 512);
    imagesc(SPCM_diff)
    title(['Difference ' Acquisition ' ' regions{indregion}], 'interpreter', 'none')
    axis off
    colorbar
    if GSR == 0
        f.CurrentAxes.CLim = [0 1];
        saveas(gcf, [SaveDir '/SPCM/Combined/withoutGSR/' dataname '_Difference_' Acquisition '_' regions{indregion} '.tiff'])
    else
        f.CurrentAxes.CLim = [0 1];
        saveas(gcf, [SaveDir '/SPCM/Combined/withGSR/' dataname '_Difference_' Acquisition '_' regions{indregion} '.tiff'])
    end
        close(f)
end

end % end function