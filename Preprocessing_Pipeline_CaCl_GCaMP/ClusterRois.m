% To make PooledROI.mat, see below function

function ClusterRois(DataFolder, Overwrite)

if ~exist('Overwrite', 'var')
    Overwrite = 0;
end

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

idx = strfind(DataFolder, filesep); %zoek alle plekken van fileseps in naam
pathROI = DataFolder(1:idx(end-2));

if exist([pathROI 'BigROI.mat'], 'file') && Overwrite == 0
    return
end


%% get ROI
load([pathROI 'ROImasks_data.mat']);
% load('/media/mbakker/GDrive/P2/GCaMP/LargerAreaAtlas.mat');
load('/home/mbakker/P2_scripts/MarleenP2/PooledROI');

BigROI = struct;
BigROI.VisualROI_R = zeros(size(img_info.imageData));
BigROI.SensoryROI_R = zeros(size(img_info.imageData));
BigROI.AuditoryROI_R = zeros(size(img_info.imageData));
BigROI.MotorROI_R = zeros(size(img_info.imageData));
BigROI.RetrosplenialROI_R = zeros(size(img_info.imageData));
BigROI.VisualROI_L = zeros(size(img_info.imageData));
BigROI.SensoryROI_L = zeros(size(img_info.imageData));
BigROI.AuditoryROI_L = zeros(size(img_info.imageData));
BigROI.MotorROI_L = zeros(size(img_info.imageData));
BigROI.RetrosplenialROI_L = zeros(size(img_info.imageData));

for idx = 1:size(ROI_info,2) % go per roi that is within mask for this mouse
    
    if contains(ROI_info(idx).Name, '_R') && ...
            sum( matches(Visual, ROI_info(idx).Name) ) > 0
        BigROI.VisualROI_R = BigROI.VisualROI_R + ROI_info(idx).Stats.ROI_binary_mask;
    elseif contains(ROI_info(idx).Name, '_L') && ...
            sum( matches(Visual, ROI_info(idx).Name) ) > 0
        BigROI.VisualROI_L = BigROI.VisualROI_L + ROI_info(idx).Stats.ROI_binary_mask;
        
    elseif contains(ROI_info(idx).Name, '_R') && ...
            sum( matches(Somatosensory, ROI_info(idx).Name) ) > 0
        BigROI.SensoryROI_R = BigROI.SensoryROI_R + ROI_info(idx).Stats.ROI_binary_mask;
    elseif contains(ROI_info(idx).Name, '_L') && ...
            sum( matches(Somatosensory, ROI_info(idx).Name) ) > 0
        BigROI.SensoryROI_L = BigROI.SensoryROI_L + ROI_info(idx).Stats.ROI_binary_mask;
        
    elseif contains(ROI_info(idx).Name, '_R') && ...
            sum( matches(Auditory, ROI_info(idx).Name) ) > 0
        BigROI.AuditoryROI_R = BigROI.AuditoryROI_R + ROI_info(idx).Stats.ROI_binary_mask;
    elseif contains(ROI_info(idx).Name, '_L') && ...
            sum( matches(Auditory, ROI_info(idx).Name) ) > 0
        BigROI.AuditoryROI_L = BigROI.AuditoryROI_L + ROI_info(idx).Stats.ROI_binary_mask;
        
    elseif contains(ROI_info(idx).Name, '_R') && ...
            sum( matches(Motor, ROI_info(idx).Name) ) > 0
        BigROI.MotorROI_R = BigROI.MotorROI_R + ROI_info(idx).Stats.ROI_binary_mask;
    elseif contains(ROI_info(idx).Name, '_L') && ...
            sum( matches(Motor, ROI_info(idx).Name) ) > 0
        BigROI.MotorROI_L = BigROI.MotorROI_L + ROI_info(idx).Stats.ROI_binary_mask;
        
    elseif contains(ROI_info(idx).Name, '_R') && ...
            sum( matches(Retrosplenial, ROI_info(idx).Name) ) > 0
        BigROI.RetrosplenialROI_R = BigROI.RetrosplenialROI_R + ROI_info(idx).Stats.ROI_binary_mask;
    elseif contains(ROI_info(idx).Name, '_L') && ...
            sum( matches(Retrosplenial, ROI_info(idx).Name) ) > 0
        BigROI.RetrosplenialROI_L = BigROI.RetrosplenialROI_L + ROI_info(idx).Stats.ROI_binary_mask;
    end
    
end

%make sure you don't have lines between originial, smaller regions.
% imagesc(imfill(bwmorph(BigROI.SensoryROI_L, 'close'), 'holes'))
% regions = fieldnames(BigROI);
regions = {'VisualROI_R', 'SensoryROI_R', 'AuditoryROI_R', 'MotorROI_R',...
    'RetrosplenialROI_R', 'VisualROI_L','SensoryROI_L', 'AuditoryROI_L',...
    'MotorROI_L','RetrosplenialROI_L' };

for ind = 1:length(regions)
    fieldname = regions{ind};
    eval(['BigROI.' fieldname ' = imfill(bwmorph(BigROI.' fieldname ', ''close''), ''holes'');']);
end

AtlasMask = zeros(size(img_info.imageData));
for ind = 1:length(regions)
    fieldname = regions{ind};
    eval(['AtlasMask(BigROI.' fieldname ') = ind;']);
end

save([pathROI 'BigROI.mat'], 'BigROI', 'AtlasMask', 'regions');

end


% % To make the clusters:
% load('/home/mbakker/P2_scripts/Umit-release_Astrocyte-v1.5/GUI/DataViz/mouse_ctx_borders.mat')
% allrois = [atlas.areatag atlas.area_longNames atlas.area_ClusterNames];
%
% ind1 = find(strcmp(allrois(:,3), 'Visual'));
% Visual = allrois(ind1,1);
% ind1 = find(strcmp(allrois(:,3), 'Motor'));
% Motor = allrois(ind1,1);
%
% % Somatosensory/multimodal
% ind1 = find(strcmp(allrois(:,3), 'Somatosensory'));
% ind2 = find(strcmp(allrois(:,3), 'Multimodal'));
% Somatosensory = [allrois(ind1,1); allrois(ind2,1)];
%
% % Retrosplenial/cingulate
% ind1 = find(strcmp(allrois(:,3), 'Retrosplenial'));
% ind2 = find(strcmp(allrois(:,3), 'Cingulate'));
% Retrosplenial = [allrois(ind1,1); allrois(ind2,1)];
%
% % Auditory/temporal association area
% ind1 = find(strcmp(allrois(:,3), 'Auditory'));
% ind2 = find(strcmp(allrois(:,3), 'Temporal association area'));
% Auditory = [allrois(ind1,1); allrois(ind2,1)];
%
% save('/home/mbakker/P2_scripts/MarleenP2/PooledROI.mat', 'allrois', 'Visual',...
%     'Auditory', 'Motor', 'Retrosplenial', 'Somatosensory');

