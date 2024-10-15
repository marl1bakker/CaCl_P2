% type is HbO or normal
% posthocoutcome comes from statsnvcinrs inside plotnvcinrs

function significance_matrix( Acq, type, SaveDir)

%% set up
if ~exist('Acq', 'var')
    Acq = 'A1';
end

if ~exist('type', 'var')
    type = 'normal';
end


if ~exist('SaveDir', 'var')
    SaveDir =  '/media/mbakker/GDrive/P2/GCaMP';
end

if ~exist([SaveDir '/NVC/Sham_RS/'], 'dir')
    mkdir([SaveDir '/NVC/Sham_RS/'])
end

load(['/media/mbakker/GDrive/P2/GCaMP/NVC/Sham_RS/Posthoc_outcome_' type '_' Acq '.mat'], 'posthocoutcome', 'cld_table')
ROI = unique(cld_table.ROI);
% names_specs = fields(posthocoutcome);
names_specs = {'GCaMPPeak', 'HbOPeak','HbRPeak', 'GCaMP_Increase', 'HbO_Increase', 'HbR_Decrease', ...
    'DelaySec', 'ResponseStrength', 'Strength_Increase'};
spectitles = {'GCaMP peak', 'HbO peak', 'HbR peak', 'GCaMP Increase', 'HbO Increase', 'HbR Decrease', 'Delay', 'Strength', 'Strength of Increase'};


%% for fig
f = figure('InvertHardcopy','off','Color',[1 1 1]);
t = tiledlayout(3,3, 'TileSpacing', 'tight', 'Padding', 'tight');

%% go per spec
for indspec = 1:length(names_specs)
    current_spec = names_specs{indspec};

    spec28 = [posthocoutcome.(current_spec).ROI_1, posthocoutcome.(current_spec).ROI_2, posthocoutcome.(current_spec).pValue];
    spec28(spec28(:,3)<0.05,3) = 0; % stat sign. is 0 -> black
    spec28(spec28(:,3)>0.05,3) = 1;

    spec64 = [1,1,NaN; spec28];
    for ind = 1:length(ROI)-1
        rowtoinsert = ind*8;
        roi1 = repmat(ind+1, ind+1, 1);
        roi2 = (1:ind+1)';
        newRows = [roi1, roi2, NaN(ind+1,1)];
        spec64 = [spec64(1:rowtoinsert,:); newRows; spec64(rowtoinsert+1:end,:)];
    end

    matrix = reshape(spec64(:,3), 8, 8);
    matrix(isnan(matrix)) = 1; % nan values should be set to not sign.

    % plot tile
    nexttile
    imagesc(matrix)
    colormap("gray")

    yticks(1:size(ROI,1));
    yticklabels(ROI);
    xticks(1:size(ROI,1));
    xticklabels(ROI);
    xtickangle(90)

    title(spectitles{indspec}, 'FontSize', 11)
    axis image
    axis equal;

end
f.Position = [100 70 800 800];
%% save
saveas(gcf, [SaveDir '/NVC/Sham_RS/' Acq '_' type '_BigROI_sign_matrix.svg'], 'svg');
saveas(gcf, [SaveDir '/NVC/Sham_RS/' Acq '_' type '_BigROI_sign_matrix.png'], 'png');
close(f)

end
