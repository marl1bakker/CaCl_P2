% Analyse results from us_analysis
% open CSV file, name "ResultsUS"

ResultsUS.Mouse = cellstr(repmat('empty', 8,1));
ResultsUS.Side = cellstr(repmat('empty', 8, 1));
ResultsUS.Sidew = repmat('empty', 8, 1);


for ind = 1:size(ResultsUS, 1)
    fullname = ResultsUS.Row{ind};
    ResultsUS.Mouse(ind,:) = cellstr(fullname(1:strfind(fullname, '-')-1));
    ResultsUS.Side(ind,:) = cellstr(fullname(strfind(fullname, '-')+1:end));
end



