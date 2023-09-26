%% timecourses for fluo and HbO
%Takes the HbO and Fluo data and the ROImasks_data and calculates the
%timecourses for the centroids of all the regions. Saves them in the folder
%of the mouse.
% option can be centroids or average over the whole region

function GetTimecourses(DataFolder, option, Overwrite, dataname)
if ~exist('option', 'var')
    option = 'centroids';
end

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

if ~exist('dataname', 'var')
    dataname = 'hemoCorr_fluo';
end

if ~exist([DataFolder dataname '.dat'], 'file')
    disp([dataname ' could not be found, function exited'])
    return
end

savename = ['timecourses_' dataname];

if ~exist('Overwrite', 'var')
    Overwrite = 0;
end

if Overwrite == 1 || ~exist([DataFolder savename '.mat'], 'file')
    % Get ROI
    seps = strfind(DataFolder, filesep);
    % if you want smaller or different ROI, change this line and the regions
    % and names later on
    % load([DataFolder(1:seps(end-2)) 'ROImasks_data.mat']);
    load([DataFolder(1:seps(end-2)) 'BigROI.mat'], 'AtlasMask', 'regions');
    
    % Get Fluo data
    fid = fopen([DataFolder dataname '.dat']);
    dat = fread(fid, inf, '*single');
    dat = reshape(dat,512,512,[]);
    fclose(fid);
    dims = size(dat);
    
    % Get timecourses Fluo
    AllRois = {};
    dat = reshape(dat,[],dims(3));
    
    for ind = 1:size(regions,2)
        Mask = ismember(AtlasMask,ind); %pak nummers van atlas van ind waar je nu bent
        
        if matches(option, 'centroids')
            Tmp = bwmorph(Mask,'shrink',inf); %maak mask kleiner, noem tmp
            Tmp = conv2(Tmp, ones(3),'same')>=1; %zorg ervoor data je alleen de ROI hebt die binnen de mask vallen die je ook hebt aangegevne bij ROI_GUI
            Mask = imerode(Mask, strel('diamond',1)) & Tmp; %krimp de ROI met 1 pixel, pak alleen pixels die ook binnen tmp vallen
        end
        
        if( sum(Mask(:)) >= 1 )
            Signal = mean(dat(Mask(:), :),1);
            name = regions{ind};
            AllRois(end+1,:) = {Mask, Signal, name}; %voeg masker, timecourse (signaal) en naam toe aan matrix.
        else % if a certain region is missing, fill with nans
            Signal = NaN(1,size(dat,2), 'single');
            name = regions{ind};
            AllRois(end+1,:) = {Mask, Signal, name};
        end
    end
    
    eval(['save(''' DataFolder savename '.mat'', ''AllRois'');']);
%     save([DataFolder 'AllRoisFluo.mat'], 'AllRoisFluo');
    clear dat i Tmp Signal name AllRois fid img_info ind seps Mask
end
end


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