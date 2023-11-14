function SingleSubjectCorrMatrix(DataFolder, Overwrite)

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

if ~exist('Overwrite', 'var')
    Overwrite = 0 ;
end

if exist([DataFolder 'CorrMatrix.tiff'], 'file') && Overwrite == 0
    disp([DataFolder ' CorrMatrix already made, exited function'])
    return
end

% Make correlation matrix
if exist([DataFolder 'AllRoisFluo.mat'], 'file')
    load([DataFolder 'AllRoisFluo.mat']);
else
    disp('Timecourses not calculated yet. Run GetTimecourses for acquisition')
    return
end

Timecourses = cell2mat(AllRoisFluo(:,2));
Names = AllRoisFluo(:,3);

CorrMat = corr(Timecourses');

%plot
f = figure('InvertHardcopy','off','Color',[1 1 1]);
ax = gca;
data = tril(CorrMat);
data(data == 0 ) = NaN;
data(data == 1 ) = NaN; %if you want to get rid of the 1's in the diagonal row
imagesc(ax, data, 'AlphaData', ~isnan(data), [-1 1])
ax.Box = 'off';
axis image;

yticks(1:size(CorrMat,1));
yticklabels(Names);
ay = get(gca,'YTickLabel');
set(gca,'YTickLabel',ay,'FontSize', 10, 'FontWeight', 'bold', 'Linewidth', 2);
xticks(1:size(CorrMat,2));
xticklabels(Names);
xtickangle(90)
load('/media/mbakker/data1/Hypoxia/SeedPixelCorrMap/NL.mat');
colormap(NL)
colorbar

seps = strfind(DataFolder, filesep);
mouse = DataFolder(seps(end-3)+1:seps(end-2)-1);
acq = DataFolder(seps(end-2)+1:seps(end-1)-1);
load('/home/mbakker/P2_scripts/MarleenP2/MiceCodes.mat')

mouseinfo = Mice(ismember(Mice.CodeOfMouse, mouse),:);

title([mouseinfo.CodeOfMouse mouseinfo.CaClSham mouseinfo.MaleFemale])
subtitle(acq)
f.Position = [10 10 1500 1500]; %for size of screen before saving
% saveas(gcf, ['/media/mbakker/data1/Hypoxia/Fluctuations/' HypoxiaLevel '_GCaMP_mean_NoGSR.eps'], 'epsc');
saveas(gcf, [DataFolder 'CorrMatrix.tiff'], 'tiff');

close(f)

end
