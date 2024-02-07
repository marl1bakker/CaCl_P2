%% timecourses for fluo and HbO and HbR
%Takes the HbO and Fluo data and the ROImasks_data and calculates the
%timecourses for the centroids of all the regions. Saves them in the folder
%of the mouse.
% option can be centroids or average over the whole region
%give dataname without .dat at the end

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

savename = ['timecourses_' dataname '_' option];

if ~exist('Overwrite', 'var')
    Overwrite = 0;
end

if Overwrite == 1 || ~exist([DataFolder savename '.mat'], 'file')
    
    % no longer have to do this, it's in umit_pipeline which checks if it's
    % normalized via tha anareg
%     % Open single frame to see if it's normalized
%     fid = fopen([DataFolder dataname '.dat']);
%     dat = fread(fid, 512*512, '*single');
%     fclose(fid);
%     dat = reshape(dat, 512*512, []);
%     if mean(dat, 'all', 'omitnan') > 1.5
% %         pause
%         disp('Data not normalized yet!! Exited function')
%         
%     else
%         disp(['running timecourses ' dataname])
%     end    
    
    % Get ROI
    seps = strfind(DataFolder, filesep);
    % if you want smaller or different ROI, change this line and the regions
    % and names later on
    load([DataFolder(1:seps(end-2)) 'BigROI.mat'], 'AtlasMask', 'regions');
    
    % Get data
%     tic
    fid = fopen([DataFolder dataname '.dat']);
    dat = fread(fid, inf, '*single');
%     toc
    dat = reshape(dat,512,512,[]);
    fclose(fid);
    dims = size(dat);
    
    % Get timecourses
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
            name = regions{ind};
            AllRois(end+1,:) = {Mask, Signal, name};
        end
    end
    
    eval(['save(''' DataFolder savename '.mat'', ''AllRois'');']);
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