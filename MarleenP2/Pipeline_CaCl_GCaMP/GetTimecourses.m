%% timecourses for fluo and HbO and HbR
%Takes the HbO and Fluo data and the ROImasks_data and calculates the
%timecourses for the centroids of all the regions. Saves them in the folder
%of the mouse.
% option can be centroids or average over the whole region
%give dataname without .dat at the end

function GetTimecourses(DataFolder, dataname, manualinput)
% if ~exist('option', 'var')
    option = 'centroids';
% end

if ~exist('manualinput', 'var')
    manualinput = 0;
end

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

if ~exist('dataname', 'var')
    dataname = 'hemoCorr_fluo';
end

if ~exist([DataFolder dataname '.dat'], 'file')
    error([dataname ' could not be found, function exited'])
elseif ~exist([DataFolder 'Seeds.mat'], 'file') && manualinput == 0
    error('Cannot do timecourses, have to run GetSeeds first')
elseif ~exist([DataFolder 'Seeds.mat'], 'file')
    GetSeeds(DataFolder)
end

savename = ['timecourses_' dataname '_' option];

%% Get ROI
seps = strfind(DataFolder, filesep);
load([DataFolder 'Seeds.mat'], 'Mask')
regions = fields(Mask);
% load([DataFolder(1:seps(end-2)) 'BigROI.mat'], 'AtlasMask', 'regions','BigROI');
load([DataFolder 'OutlierMask.mat'], 'OutlierPixels');
eval(['OutlierPixels = OutlierPixels.' dataname ';']);

%% Get data
fid = fopen([DataFolder dataname '.dat']);
dat = fread(fid, inf, '*single');
dat = reshape(dat,512,512,[]);
fclose(fid);
dims = size(dat);

% Get data sanity check
fid = fopen([DataFolder 'fluo_567.dat']);
fl = fread(fid, 512*512, '*single');
fl = reshape(fl, 512, 512);
fclose(fid);
fid = fopen([DataFolder 'green.dat']);
gr = fread(fid, 512*512, '*single');
gr = reshape(gr, 512, 512);
fclose(fid);
MaskAllCentroids = zeros(512, 512);

%% Get timecourses
AllRois = {};
dat = reshape(dat,[],dims(3));

for ind = 1:size(regions, 1) 
    name = regions{ind};

    % disp(name)
    eval(['seedmask = Mask.' name ';'])
    seedmask = logical(seedmask.*OutlierPixels);

    if( sum(seedmask(:)) >= 1 )
        Signal = mean(dat(seedmask(:), :),1);
        AllRois(end+1,:) = {seedmask, Signal, name}; %voeg masker, timecourse (signaal) en naam toe aan matrix.

        % Sanity check
        if (matches(dataname, 'hemoCorr_fluo')) && ...
                (mean(Signal) < 0.99 || mean(Signal) > 1.01)
            disp(['*** GETTIMECOURSES NOT AV 1 FOR fluo ' DataFolder])
            error('Signal not centered on 1')
        elseif (matches(dataname, 'HbO') || matches(dataname, 'HbR')) && ...
                (mean(Signal) < -0.1|| mean(Signal) > 1.1)
            plot(Signal)
            disp(['*** GETTIMECOURSES NOT AV 1 FOR HbO or HbR ' DataFolder])
            %                 error('Signal not centered on 1')
        end

    else % if a certain region is missing, fill with nans
        Signal = NaN(1,size(dat,2), 'single');
        AllRois(end+1,:) = {seedmask, Signal, name};
    end
end

%% Save
save([DataFolder savename '.mat'], 'AllRois');

end




% OLD, changed 17-4-24
% for ind = 1:size(regions,2)
%     Mask = ismember(AtlasMask,ind); %pak nummers van atlas van ind waar je nu bent
% 
%     if matches(option, 'centroids')
%         Tmp = bwmorph(Mask,'shrink',inf); %maak mask kleiner, noem tmp
%         Tmp = conv2(Tmp, ones(3),'same')>=1; %zorg ervoor data je alleen de ROI hebt die binnen de mask vallen die je ook hebt aangegevne bij ROI_GUI
%         Mask = imerode(Mask, strel('diamond',1)) & Tmp; %krimp de ROI met 1 pixel, pak alleen pixels die ook binnen tmp vallen
%         Mask = logical(Mask.*OutlierPixels);
%         MaskAllCentroids = MaskAllCentroids+Mask;
%     end
% 
%     if( sum(Mask(:)) >= 1 )
%         Signal = mean(dat(Mask(:), :),1);
%         name = regions{ind};
%         AllRois(end+1,:) = {Mask, Signal, name}; %voeg masker, timecourse (signaal) en naam toe aan matrix.
% 
%         % Sanity check
%         if (matches(dataname, 'hemoCorr_fluo')) && ...
%                 (mean(Signal) < 0.99 || mean(Signal) > 1.01)
%             disp(['*** GETTIMECOURSES NOT AV 1 FOR fluo ' DataFolder])
%             error('Signal not centered on 1')
%         elseif (matches(dataname, 'HbO') || matches(dataname, 'HbR')) && ...
%                 (mean(Signal) < -0.1|| mean(Signal) > 1.1)
%             plot(Signal)
%             disp(['*** GETTIMECOURSES NOT AV 1 FOR HbO or HbR ' DataFolder])
%             %                 error('Signal not centered on 1')
%         end
% 
%     else % if a certain region is missing, fill with nans
%         Signal = NaN(1,size(dat,2), 'single');
%         name = regions{ind};
%         AllRois(end+1,:) = {Mask, Signal, name};
%     end
% end


%
% if Overwrite == 1 || ~exist([DataFolder 'AllRoisHbO.mat'], 'file')
%
%     % Get HbO data
%     if ~exist([DataFolder 'HbO.dat'], 'file')
%         disp('HbO.dat not found')
%         return
%     end
%
%     fid = fopen([DataFolder 'HbO.dat']);
%     HbO = fread(fid, inf, '*single');
%     HbO = reshape(HbO, [], dims(3));
%     fclose(fid);
%
%     % Get timecourses
%     AllRoisHbO = {};
%     HbO = reshape(HbO,[],dims(3));
%     % TimecoursesHbO = [];
%
%     for ind = 1:size(regions,2)
%         Mask = ismember(AtlasMask,ind); %pak nummers van atlas van i waar je nu bent
%
%         if matches(option, 'centroids')
%             Tmp = bwmorph(Mask,'shrink',inf); %maak mask kleiner, noem tmp
%             Tmp = conv2(Tmp, ones(3),'same')>=1; %zorg ervoor data je alleen de ROI hebt die binnen de mask vallen die je ook hebt aangegevne bij ROI_GUI
%             Mask = imerode(Mask, strel('diamond',1)) & Tmp; %krimp de ROI met 1 pixel, pak alleen pixels die ook binnen tmp vallen
%         end
%
%         if( sum(Mask(:)) >= 1 )
%             Signal = mean(HbO(Mask(:), :),1); %pak 5 punten om midden van ROI heen, bereken timecourse van deze seed
%             name = regions{ind};
%             AllRoisHbO(end+1,:) = {Mask, Signal, name}; %voeg masker, timecourse (signaal) en naam toe aan matrix.
%         end
%     end
%
%     save([DataFolder 'AllRoisHbO.mat'], 'AllRoisHbO');
% end
%
% end