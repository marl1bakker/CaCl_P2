%% timecourses for fluo and HbO
%Takes the HbO and Fluo data and the ROImasks_data and calculates the
%timecourses for the centroids of all the regions. Saves them in the folder
%of the mouse. 

function GetTimecourses(DataFolder)

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

% Get ROI
seps = strfind(DataFolder, filesep);
load([DataFolder(1:seps(end-2)) 'ROImasks_data.mat']);

% Get Fluo data
% fid = fopen([DataFolder 'fluo_567.dat']);
% fid = fopen([DataFolder 'normLPF.dat']);
fid = fopen([DataFolder 'mov_aligned.dat']);
dat = fread(fid, inf, '*single');
dat = reshape(dat,512,512,[]);
fclose(fid);
dims = size(dat);

% Combine ROI areas into a single mask (AtlasMask)
AtlasMask = zeros(dims(1),dims(2));
for ind = 1:size(ROI_info,2)
    AtlasMask(ROI_info(ind).Stats.ROI_binary_mask) = ind;
end

% Get timecourses Fluo
AllRoisFluo = {};
dat = reshape(dat,[],dims(3));
% TimecoursesFluo = [];

for i = unique(nonzeros(AtlasMask(:)))' %pak alleen waarden die ook echt op de atlas staan en binnen brein vallen
    Mask = ismember(AtlasMask,i); %pak nummers van atlas van i waar je nu bent
    Tmp = bwmorph(Mask,'shrink',inf); %maak mask kleiner, noem tmp
    Tmp = conv2(Tmp, ones(3),'same')>=1; %zorg ervoor data je alleen de ROI hebt die binnen de mask vallen die je ook hebt aangegevne bij ROI_GUI
    Mask = imerode(Mask, strel('diamond',1)) & Tmp; %krimp de ROI met 1 pixel, pak alleen pixels die ook binnen tmp vallen
    if( sum(Mask(:)) >= 1 )
        Signal = mean(dat(Mask(:), :),1); %pak 5 punten om midden van ROI heen, bereken timecourse van deze seed
        name = ROI_info(i).Name; %pak namen van ROI_GUI
        AllRoisFluo(end+1,:) = {Mask, Signal, name}; %voeg masker, timecourse (signaal) en naam toe aan matrix.
%         TimecoursesFluo(:,end+1) = Signal;
    end
end

save([DataFolder 'AllRoisFluo.mat'], 'AllRoisFluo');
clear dat TimecoursesFluo i Tmp Signal name AllRoisFluo fid img_info ind seps Mask


% Get HbO data
if ~exist([DataFolder 'HbO.dat'], 'file')
    return
end

fid = fopen([DataFolder 'HbO.dat']);
HbO = fread(fid, inf, '*single');
HbO = reshape(HbO, [], dims(3));
fclose(fid);

% Get timecourses Fluo
AllRoisHbO = {};
HbO = reshape(HbO,[],dims(3));
% TimecoursesHbO = [];

for i = unique(nonzeros(AtlasMask(:)))' %pak alleen waarden die ook echt op de atlas staan en binnen brein vallen
    Mask = ismember(AtlasMask,i); %pak nummers van atlas van i waar je nu bent
    Tmp = bwmorph(Mask,'shrink',inf); %maak mask kleiner, noem tmp
    Tmp = conv2(Tmp, ones(3),'same')>=1; %zorg ervoor data je alleen de ROI hebt die binnen de mask vallen die je ook hebt aangegevne bij ROI_GUI
    Mask = imerode(Mask, strel('diamond',1)) & Tmp; %krimp de ROI met 1 pixel, pak alleen pixels die ook binnen tmp vallen
    if( sum(Mask(:)) >= 1 )
        Signal = mean(HbO(Mask(:), :),1); %pak 5 punten om midden van ROI heen, bereken timecourse van deze seed
        name = ROI_info(i).Name; %pak namen van ROI_GUI
        AllRoisHbO(end+1,:) = {Mask, Signal, name}; %voeg masker, timecourse (signaal) en naam toe aan matrix.
%         TimecoursesHbO(:,end+1) = Signal;
    end
end

save([DataFolder 'AllRoisHbO.mat'], 'AllRoisHbO');

end