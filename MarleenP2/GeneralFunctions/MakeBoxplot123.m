% Same as MakeBoxplot, but it takes all 3 acquisitions and plots them in a
% single figure, below each other. 

% Make sure in your overviewtable that you have:
% Mouse, some variable for x (usually ROI), Group, Sex, Combi (or at least the
% grouping that you want to use) and the y value that you want plotted.
%

function [f, t] = MakeBoxplot123(Grouping, groups, overviewtable, ylimvalues, x, y)

if ~exist('x','var')
    x = 'ROI';
end

%name the columns in your table that you want to use x and y
temp = find(matches(overviewtable.Properties.VariableNames, x));
overviewtable.Properties.VariableNames{temp} = 'x';

if exist('y', 'var')
    temp = find(matches(overviewtable.Properties.VariableNames, y));
    overviewtable.Properties.VariableNames{temp} = 'y';
end

overviewtable.idx = grp2idx(overviewtable.x); %this gives the groups based on alphabet, so sort the ROI labels as well:
labels = cellstr(unique(overviewtable.x))';

f = figure('InvertHardcopy','off','Color',[1 1 1]);
t = tiledlayout(3,1);
acquisitions = {'A1', 'A2', 'A3'};

for indacq = 1:size(acquisitions, 2)
    Acq = overviewtable(overviewtable.Acquisition == acquisitions{indacq},:);
    axes = nexttile(t);
    
    % BOXPLOT
    eval(['b = boxchart(Acq.idx, Acq.y, ''GroupByColor'', Acq.' Grouping ', ''LineWidth'', 2, ''MarkerStyle'', ''none'');'])
    hold on
    
    % SCATTER
    xaxisstep = 1/size(groups,1);
    xaxisplacement = 1 - 0.5*xaxisstep - (size(groups,1)/2-1)*xaxisstep - xaxisstep;
    
    for indroi = 1:length(labels)
        currentROI = Acq.x == labels{indroi}; 
        
        for indgroup = 1:size(groups,1)
            eval(['currentgroup = Acq.' Grouping ' == groups{indgroup};']);
            currentindex = currentgroup.*currentROI;
            Acq.idx(currentindex == 1) = xaxisplacement + xaxisstep*indgroup;
        end
        
        xaxisplacement = xaxisplacement + 1;
    end
    clear indroi indgroup currentgroup currentROI currentindex ind indmouse tablemouse xaxisplacement group idx Mousegroup toplot
    
    hs = scatter(Acq.idx, Acq.y, 30, 'filled', 'jitter','on','JitterAmount',0.02);
    
    % MAKE PRETTY
    hs.MarkerFaceColor = [0 0 0];
    hs.MarkerFaceAlpha = 0.3;
    xticks(1:length(labels));
    xticklabels(labels);
    xlim([0.2 length(labels)+0.7])
    
    axes1 = b.Parent;
    hold(axes1,'on');
%     set(axes1,'FontSize',15,'FontWeight','bold','LineWidth',2);
    set(axes1,'FontSize',15,'LineWidth',2);
%     legend(b, 'Location', 'northeast', 'NumColumns', 2);
    ylim(ylimvalues);
    
    b(1).SeriesIndex = 7;
    b(2).SeriesIndex = 1;
    if size(b,1)>2
        b(3).SeriesIndex = 2;
        b(4).SeriesIndex = 6;
    end
    
    title(['Acquisition ' num2str(indacq)], 'FontSize', 15, 'FontWeight', 'Normal');
%     t.Title(acquisitions{indacq}, 'interpreter', 'none','FontSize',20,'FontWeight','bold')

end

% title(t, 'xxx');

subtitleN = [];
for indgroup = 1:size(groups, 1)
    eval(['ngroup = sum(overviewtable.' Grouping ' == groups{indgroup});'])
    ngroup = ngroup/length(labels)/3; %because you have a row for each ROI and you count them all. divide by 8 and you have the nr of mice. Divide by 3 for 3 acq
    subtitleN = [subtitleN groups{indgroup} ' N = ' num2str(ngroup) ' -- '];
end

subtitle(t,subtitleN(1:end-4), 'FontSize', 10)
% title(t, 'Fluctuations', 'FontSize', 20, 'FontWeight', 'bold');

% t = b; %to give return

end