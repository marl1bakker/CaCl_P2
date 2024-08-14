% note: tried this after everything else, so it's not embedded as nicely.
% just load the overview table and work from there

% type = 'normal' or 'random'

function NumberOfActivations(type, Overwrite)

% set up
if~exist('type','var')
    type = 'normal';
end

if ~exist('Overwrite', 'var')
    Overwrite = 0;
end

load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat', 'RecordingOverview');
Recordings = [RecordingOverview.A1; RecordingOverview.A2; RecordingOverview.A3];
Mice = [RecordingOverview.Mouse; RecordingOverview.Mouse; RecordingOverview.Mouse];
% Acqs = [repmat({'A1'}, size(Recordings,1)/3, 1); repmat({'A2'}, size(Recordings,1)/3, 1); repmat({'A3'}, size(Recordings,1)/3, 1)];
SaveDirs = [RecordingOverview.SaveDirectory; RecordingOverview.SaveDirectory;RecordingOverview.SaveDirectory;];

if matches(type,'normal') && exist('/media/mbakker/GDrive/P2/GCaMP/NVC/NrOfActivations.mat', 'file') 
    load('/media/mbakker/GDrive/P2/GCaMP/NVC/NrOfActivations.mat', 'ActsTable')
    ActsTable.Acq = cellstr(ActsTable.Acq);
elseif matches(type,'random') && exist('/media/mbakker/GDrive/P2/GCaMP/NVC/NrOfActivations_random.mat', 'file') 
    load('/media/mbakker/GDrive/P2/GCaMP/NVC/NrOfActivations_random.mat', 'ActsTable')
    ActsTable.Acq = cellstr(ActsTable.Acq);
end

if matches(type, 'normal') && exist('/media/mbakker/GDrive/P2/GCaMP/NVC/PercentageMatches.mat', 'file')
    load('/media/mbakker/GDrive/P2/GCaMP/NVC/PercentageMatches.mat', 'PercTable')
    PercTable.Acq = cellstr(PercTable.Acq);
elseif matches(type,'random') && exist('/media/mbakker/GDrive/P2/GCaMP/NVC/PercentageMatches_random.mat', 'file')
    load('/media/mbakker/GDrive/P2/GCaMP/NVC/PercentageMatches_random.mat', 'PercTable')
    PercTable.Acq = cellstr(PercTable.Acq);
    % else
    %     PercTable = table;
    %     PercTable.Mouse = {'dummy'};
    %     PercTable.Acq = {'xx'};
    %     PercTable.Group = categorical({'empty'});
    %     PercTable.DetectedOn = {'xx'};
end

%% go per recording
for indRec = 1:size(Recordings,1)
    % set up
    DataFolder = Recordings{indRec};
    Mouse = Mice{indRec};
    fs = strfind(DataFolder, '-');
    Acq = {DataFolder(fs(end-1)+1:fs(end)-1)};
    clear fs
    if( ~strcmp(DataFolder(end), filesep) )
        DataFolder = [DataFolder filesep];
    end
    DataFolder = [SaveDirs{indRec} filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];

    if matches(type, 'normal')
        do_nr_act = 1;
        do_act_perc = 1;
    elseif matches(type, 'random')
        do_nr_act = 0; 
        do_act_perc = 1;
    end

    disp([Mouse ' ' Acq{1}])

    % %% temp
    % if exist([DataFolder 'ActivationPercentages.mat'], 'file')
    %     delete([DataFolder 'ActivationPercentages.mat'])
    % end

    try
    %check if already done
    if exist('ActsTable', 'var')
        indmouse = find(ismember(ActsTable.Mouse, Mouse).*...
            ismember(ActsTable.Acq, Acq).*ismember(ActsTable.DetectedOn, 'HbO'));
        temp_nrofacts = ActsTable(indmouse,:);
        if Overwrite == 0 && size(temp_nrofacts,1)>0
            disp('nr of acts mouse already done')
            do_nr_act = 0;
        end
    end

    if exist("PercTable",'var')
        indmouse = find(ismember(PercTable.Mouse, Mouse).*...
            ismember(PercTable.Acq, Acq).*ismember(PercTable.DetectedOn, 'HbO'));
        temp_perc = PercTable(indmouse,:);
        if Overwrite == 0 && size(temp_perc,1)>0
            disp('percentage match mouse already done')
            do_act_perc = 0;
        end
    end

    if matches(Mouse, 'M23')
        do_nr_act = 0;
        do_act_perc = 0;
        % continue %M23 has not enough fluorescence, so skip
    elseif matches(Mouse, 'M14') && matches(Acq, 'A3')
        do_nr_act = 0;
        do_act_perc = 0;
        % continue %M14 has really weird outliers for hbo/hbr data
    elseif matches(Mouse, 'M32') && (matches(Acq, 'A2') || matches(Acq, 'A3'))
        do_nr_act = 0;
        do_act_perc = 0;
        % continue %M32 has damaged window, hbo/hbr data looks very bizarre
    end
    clear temp* indmouse


    %% start
    if do_nr_act == 1 || do_act_perc == 1
        datatypes = {'hemoCorr_fluo', 'HbO'};

        for inddata = 1:size(datatypes,2)
            datanameactivations = datatypes{inddata};

            % load data
            fid = fopen([DataFolder datanameactivations '.dat']);
            dF = fread(fid, inf, '*single');
            dF = reshape(dF, 512,512,[]);
            fclose(fid);

            load([DataFolder 'OutlierMask.mat'], 'OutlierFrames', 'OutlierPixels');
            seps = strfind(DataFolder, filesep);
            load([DataFolder(1:seps(end-2)) 'ROImasks_data.mat'], 'img_info');
            mask = img_info.logical_mask;

            dF = reshape(dF, 512*512, []);
            dF = dF .* OutlierFrames;
            dF = reshape(dF, 512, 512, []);
            dF = dF .* OutlierPixels.(datanameactivations);
            dF = dF .* mask;
            dF(dF == 0) = NaN;

            %Zscore:
            zF = (dF - mean(dF, 3, 'omitmissing'))./std(dF,0,3, 'omitmissing');
            %Threshold on Zscore:
            aF = zF >= 1.95;
            % Removing noise:
            for ind = 1:size(aF,3)
                aF(:,:,ind) = bwmorph(bwmorph(aF(:,:,ind), 'close', inf),'open',inf);
            end
            %Now, we want only the beginning of activations:
            aF = aF(:,:,2:end)&~aF(:,:,1:(end-1));
            aF = cat(3, false(size(aF,1),size(aF,2)), aF);

            clear zF ind

            %% Get mask for brain and movement and outliers
            OutlierPix = OutlierPixels.HbO + OutlierPixels.HbR + OutlierPixels.hemoCorr_fluo;
            OutlierPix(OutlierPix<3) = 0;
            OutlierPix(OutlierPix == 3) = 1;

            mask = mask.*logical(OutlierPix);
            mask = logical(mask);

            %take only brain
            aF = reshape(aF,[], size(aF,3));
            aF = aF(mask(:),:);

            %exclude movement and outliers
            load([DataFolder 'MovMask.mat'], 'MovMask');
            aF = aF.* MovMask;
            aF = aF.*OutlierFrames;

            clear dF OutlierPix

            if do_nr_act == 1
                %% Number of activations
                MouseTable = table;
                MouseTable.Mouse = Mice(indRec);
                MouseTable.Acq = Acq;
                MouseTable.Group = RecordingOverview.Group(matches(RecordingOverview.Mouse, Mice{indRec}), :);
                MouseTable.DetectedOn = {datanameactivations};

                % get nr of activations per pixel
                actmap = zeros(512, 512);
                actmap(mask(:)) = sum(aF,2);

                clear img_info MovMask seps OutlierMask dF fid MovMask OutlierFrames OutlierPix OutlierPixels 

                %% go per brain area
                MouseTable.WholeBrain = sum(actmap, 'all');

                % Specific areas (ROI)
                seps = strfind(DataFolder, filesep);
                load([DataFolder(1:seps(end-2)) 'BigROI.mat'], 'regions', 'BigROI', 'AtlasMask');

                for indroi = 1:size(regions,2) %go per roi
                    mapROI = BigROI.(regions{indroi});
                    MouseTable.(regions{indroi}) = sum(actmap.*mapROI, 'all');
                end

                if ~exist('ActsTable', 'var')
                    ActsTable = MouseTable;
                else
                    ActsTable = [ActsTable; MouseTable];
                end
            end

            %% save activations for percentage match
            temp = zeros(512*512,size(aF,2));
            temp(mask(:),:) = aF;
            detected_activations.(datanameactivations) = temp;

            clear aF temp
        
            %% if case is random
            if matches(type, 'random')
                aF = detected_activations.(datanameactivations);
                nrofacts = sum(aF, 2);
                aF_random = zeros(size(aF));
                for ind_random = 1:size(aF,1) % go per pixel
                    if nrofacts(ind_random) == 0
                        continue
                    end
                    random_activation_indices = round(rand(nrofacts(ind_random),1)*(size(aF,2)-1))+1;
                    aF_random(ind_random,random_activation_indices) = 1;
                end
            end
            random_detected_activations.(datanameactivations) = aF_random; % save random aF 
            clear aF temp aF_random nrofacts random_activation_indices 
        end

        %% percentage match -- normal
        if do_act_perc == 1 

            % if you do random activations, you still want the other data
            % to be "real". For example, if you take random hbo
            % activations, you want to check how many times this matches
            % with fluo, so you want to keep the fluo data not random. This
            % is why there's going to be a used_activations.xxx variable.
            if matches(type, 'normal')
                used_activations = detected_activations;
            elseif matches(type, 'random')
                used_activations = random_detected_activations;
            end

            %% compare the activations of fluo and hbo
            fluo_detected.match = zeros(512*512,1);
            fluo_detected.dontmatch = zeros(512*512,1);

            hbo_detected.match = zeros(512*512,1);
            hbo_detected.dontmatch = zeros(512*512,1);

            imfreq = 15;

            for indpix = 1:512*512
                if sum(detected_activations.hemoCorr_fluo(indpix,:), 'all') == 0 ||...
                        sum(detected_activations.HbO(indpix,:), 'all') == 0
                    % either outlier or not a brain pixel, skip. Code is not perfect,
                    % could be a fluo activation and just no hbo activation for the
                    % whole acq, but seems unlikely
                    continue
                end

                %% detected on fluo
                for fluo_act_ind = find(used_activations.hemoCorr_fluo(indpix,:))

                    % define hbo response as detected hbo activation within 2 seconds
                    % after fluo activation
                    if fluo_act_ind + imfreq*2 > size(used_activations.hemoCorr_fluo,2) %if detected act is less than 2 sec before acq ends
                        continue
                    end

                    hbo_respond = detected_activations.HbO(indpix, fluo_act_ind:fluo_act_ind+imfreq*2);
                    if sum(hbo_respond)>0
                        fluo_detected.match(indpix) = fluo_detected.match(indpix) + 1;
                    else
                        fluo_detected.dontmatch(indpix) = fluo_detected.dontmatch(indpix) + 1;
                    end
                end

                %% detected on hbo
                for hbo_act_ind = find(used_activations.HbO(indpix,:))
                    % define fluo response as detected fluo activation within 2 seconds
                    % before hbo activation
                    if hbo_act_ind - imfreq*2 < 1 %if detected act is less than 2 sec into acquisition
                        continue
                    end

                    fluo_respond = detected_activations.hemoCorr_fluo(indpix, hbo_act_ind-imfreq*2:hbo_act_ind);
                    if sum(fluo_respond)>0
                        hbo_detected.match(indpix) = hbo_detected.match(indpix) + 1;
                    else
                        hbo_detected.dontmatch(indpix) = hbo_detected.dontmatch(indpix) + 1;
                    end
                end

            end


            % to exclude pixels without activations
            temp = zeros(512*512,1);
            temp(fluo_detected.match>0) = 1;
            temp(fluo_detected.dontmatch>0) = 1;
            fluo_detected.match(temp == 0) = NaN;

            temp = zeros(512*512,1);
            temp(hbo_detected.match>0) = 1;
            temp(hbo_detected.dontmatch>0) = 1;
            hbo_detected.match(temp == 0) = NaN;

            Percentage_fluo = reshape(fluo_detected.match./(fluo_detected.match+fluo_detected.dontmatch), 512, 512);
            Percentage_hbo = reshape(hbo_detected.match./(hbo_detected.match+hbo_detected.dontmatch), 512, 512);

            clear hbo_detected fluo_detected

            %% get values Whole brain
            MouseTablePerc = table;

            MouseTablePerc.Mouse = Mice(indRec);
            MouseTablePerc.Acq = Acq;
            MouseTablePerc.Group = RecordingOverview.Group(matches(RecordingOverview.Mouse, Mice{indRec}), :);
            MouseTablePerc.DetectedOn = {'hemoCorr_fluo'};
            MouseTablePerc.WholeBrain = mean(Percentage_fluo, 'all', 'omitnan');

            warning off
            MouseTablePerc.Mouse(2) = Mice(indRec);
            MouseTablePerc.Acq(2) = Acq;
            MouseTablePerc.Group(2) = RecordingOverview.Group(matches(RecordingOverview.Mouse, Mice{indRec}), :);
            MouseTablePerc.DetectedOn(2) = {'HbO'};
            MouseTablePerc.WholeBrain(2) = mean(Percentage_hbo, 'all', 'omitnan');
            warning on

            %% Specific areas (ROI)
            seps = strfind(DataFolder, filesep);
            load([DataFolder(1:seps(end-2)) 'BigROI.mat'], 'regions', 'BigROI', 'AtlasMask');

            for indroi = 1:size(regions,2) %go per roi
                mapROI = BigROI.(regions{indroi});
                Percentage_ROI_fluo = Percentage_fluo(mapROI==1);
                Percentage_ROI_hbo = Percentage_hbo(mapROI==1);
                % MouseTable.(regions{indroi}) = sum(actmap.*mapROI, 'all');
                MouseTablePerc.(regions{indroi})(1,:) = mean(Percentage_ROI_fluo, 'all', 'omitnan');
                MouseTablePerc.(regions{indroi})(2,:) = mean(Percentage_ROI_hbo, 'all', 'omitnan');

            end

            if ~exist('PercTable', 'var')
                PercTable = MouseTablePerc;
            else
                PercTable = [PercTable; MouseTablePerc];
            end
        end
    end

    catch
        disp(['Something went wrong with Mouse ' Mouse])
        pause
    end

    clear do_act_perc do_nr_act MouseTable MouseTablePerc MovMask OutlierPix temp_nrofacts temp_perc zF
end


%% save
if matches(type, 'normal')
    save('/media/mbakker/GDrive/P2/GCaMP/NVC/NrOfActivations.mat', 'ActsTable');
    save('/media/mbakker/GDrive/P2/GCaMP/NVC/PercentageMatches.mat', 'PercTable');
elseif matches(type, 'random')
    save('/media/mbakker/GDrive/P2/GCaMP/NVC/PercentageMatches_random.mat', 'PercTable');
end

end
