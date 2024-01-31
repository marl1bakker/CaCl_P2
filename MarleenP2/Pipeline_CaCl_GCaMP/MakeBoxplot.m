% Make sure in your overviewtable that you have:
% Mouse, some variable for x (usually ROI), Group, Sex, Combi (or at least the 
% grouping that you want to use) and the y value that you want plotted.
% 

function [f, t] = MakeBoxplot(Grouping, groups, overviewtable, ylimvalues, x, y)

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

if matches(Grouping, 'Sex') % if you're grouping by sex

    % BOXPLOT
    t = tiledlayout(2,1);
    conditions = {'CaCl','Sham'};
    
    %go per tile, first CaCl then Sham
    for indcondition = 1:size(conditions,2)
        Cond = overviewtable(overviewtable.Group == conditions{indcondition},:);
        
        axes = nexttile(t);
        b = boxchart(Cond.idx, Cond.y, 'GroupByColor', Cond.Sex,...
            'LineWidth', 2, 'MarkerStyle', 'none');
        hold on
        
        if indcondition == 1
            b(1).SeriesIndex = 7;
            b(2).SeriesIndex = 1;
        elseif indcondition == 2
            b(1).SeriesIndex = 2;
            b(2).SeriesIndex = 6;
        end
        
        % SCATTER
        xaxisstep = 1/size(groups,1);
        xaxisplacement = 1 - 0.5*xaxisstep - (size(groups,1)/2-1)*xaxisstep - xaxisstep;
        
        for indroi = 1:length(labels)
            currentROI = Cond.x == labels{indroi};
            
            for indgroup = 1:size(groups,1)
                currentgroup = Cond.Sex == groups{indgroup};
                currentindex = currentgroup.*currentROI;
                Cond.idx(currentindex == 1) = xaxisplacement + xaxisstep*indgroup;
            end
            
            xaxisplacement = xaxisplacement + 1;
        end
        clear indroi indgroup currentgroup currentROI currentindex xaxisplacement
        
        hs = scatter(Cond.idx, Cond.y, 70, 'filled', 'jitter','on','JitterAmount',0.02);
        
        % MAKE PRETTY
        ylim(ylimvalues);
        set(axes,'FontSize',20,'FontWeight','bold','LineWidth',2);
        legend('Female','Male','', 'Location', 'northeast', 'NumColumns', 2)
        hs.MarkerFaceColor = [0 0 0];
        hs.MarkerFaceAlpha = 0.3;
        hs.Marker;
        xticks(1:length(labels));
        xticklabels(labels)
        xlim([0.2 length(labels)+0.7])
        title(conditions{indcondition}, 'FontSize', 18)
        
        subtitleN = [];
        for indgroup = 1:size(groups, 1)
            ngroup = sum(Cond.Sex == groups{indgroup});
            ngroup = ngroup/length(labels);
            subtitleN = [subtitleN groups{indgroup} ' N = ' num2str(ngroup) ' -- '];
        end
        subtitle(axes, subtitleN(1:end-4), 'FontSize', 15)
    end
    
else
    
    % BOXPLOT
    %dont display outliers because we will do scatter that will show them
    eval(['b = boxchart(overviewtable.idx, overviewtable.y, ''GroupByColor'', overviewtable.' Grouping ', ''LineWidth'', 2, ''MarkerStyle'', ''none'');'])
    hold on
    
    % SCATTER
    xaxisstep = 1/size(groups,1);
    xaxisplacement = 1 - 0.5*xaxisstep - (size(groups,1)/2-1)*xaxisstep - xaxisstep;
    
    for indroi = 1:length(labels)
        currentROI = overviewtable.x == labels{indroi}; %is called currentroi because it's usually grouped by ROI, but can be something else
        
        for indgroup = 1:size(groups,1)
            eval(['currentgroup = overviewtable.' Grouping ' == groups{indgroup};']);
            currentindex = currentgroup.*currentROI;
            overviewtable.idx(currentindex == 1) = xaxisplacement + xaxisstep*indgroup;
        end
        
        xaxisplacement = xaxisplacement + 1;
    end
    clear indroi indgroup currentgroup currentROI currentindex ind indmouse tablemouse xaxisplacement group idx Mousegroup toplot
    
    hs = scatter(overviewtable.idx, overviewtable.y, 70, 'filled', 'jitter','on','JitterAmount',0.02);
    
    % MAKE PRETTY
    hs.MarkerFaceColor = [0 0 0];
    hs.MarkerFaceAlpha = 0.3;
    xticks(1:length(labels));
    xticklabels(labels);
    xlim([0.2 length(labels)+0.7])
    
    axes1 = b.Parent;
    hold(axes1,'on');
    set(axes1,'FontSize',20,'FontWeight','bold','LineWidth',2);
    legend(b, 'Location', 'northeast', 'NumColumns', 2);
%     ylim([0 100]);
    ylim(ylimvalues);
     
    b(1).SeriesIndex = 7;
    b(2).SeriesIndex = 1;
    if size(b,1)>2
        b(3).SeriesIndex = 2;
        b(4).SeriesIndex = 6;
    end
    
%     title(['Connectivity ' dataname(1:end-4) ' ' Acquisition], 'interpreter', 'none','FontSize',20,'FontWeight','bold')
    subtitleN = [];
    for indgroup = 1:size(groups, 1)
        eval(['ngroup = sum(overviewtable.' Grouping ' == groups{indgroup});'])
        ngroup = ngroup/length(labels); %because you have a row for each ROI and you count them all. divide by 8 and you have the nr of mice
        subtitleN = [subtitleN groups{indgroup} ' N = ' num2str(ngroup) ' -- '];
    end
    subtitle(subtitleN(1:end-4), 'FontSize', 15)
    t = b; %to give return
end
end