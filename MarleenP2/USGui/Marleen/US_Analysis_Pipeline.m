%% Analysis in the app
% First, do us_analysis via the GUI. Save the results in a CSV file. Import
% the CSV file into matlab, and call it "ResultsUS". This is important to
% make the table later. 

%% Make a nice table
load('C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\ResultsUS.mat');

load('C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\Marleen\Mice.mat')

ResultsUS.Mouse = cellstr(repmat('empty', size(ResultsUS,1),1));
ResultsUS.Side = cellstr(repmat('empty', size(ResultsUS,1), 1));
ResultsUS.CaClSham = cellstr(repmat('empty', size(ResultsUS,1),1));
ResultsUS.MaleFemale = cellstr(repmat('empty', size(ResultsUS,1),1));

for ind = 1:size(ResultsUS, 1)
    fullname = ResultsUS.Row{ind};
    ResultsUS.Mouse(ind,:) = cellstr(fullname(1:strfind(fullname, '-')-1));
    ResultsUS.Side(ind,:) = cellstr(fullname(strfind(fullname, '-')+1:end));

    indmouse = find(ismember(Mice.CodeOfMouse, ResultsUS.Mouse(ind))); %find mouse
    ResultsUS.CaClSham(ind,:) = cellstr(Mice.CaClSham(indmouse));
    ResultsUS.MaleFemale(ind,:) = cellstr(Mice.MaleFemale(indmouse));
end

clear ind indmouse fullname Mice

ResultsUS.X = single(ismember(ResultsUS.Side, 'R')); %L = 0, R = 1

save('C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\ResultsUS.mat');

%% Analysis "by hand"
% For acquisitions that cannot be neatly analysed via the GUI, go to the 
% file US_analysis_by_hand. Pick out the right recordings and add them
% manually to the table. 
load('C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\ResultsUS+byhand.mat', 'ResultsUS');

%% Plot the results. 
PlotUSResults(ResultsUS)

saveas(gcf, 'C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\Marleen\UltrasoundPlots.tiff', 'tiff');
saveas(gcf, 'C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\Marleen\UltrasoundPlots.epsc', 'epsc');





