function [Stats] = StatsNVC(allSpecs)
load('/media/mbakker/GDrive/temp.mat', 'allSpecs');

% get groups
% CaCl = allSpecs(allSpecs.Group == 'CaCl',:);
% Sham = allSpecs(allSpecs.Group == 'Sham',:);

% continuous data
specs = {'ResponseStrength', 'GCaMPPeak', 'HbOPeak','HbRPeak'};

for indspec = 1:size(specs,2)
    eval(['spec = allSpecs.' specs{indspec} ';'])

    hnorm = adtest(spec);% test normalcy
    if hnorm == 0 % check homogeneity of variance
        hvar = vartestn(spec, 'Display', 'off');
        if hvar == 0
            [p, h] = ttest2(spec(allSpecs.Group == 'CaCl'), spec(allSpecs.Group == 'Sham'));
        end

    elseif hnorm == 1 || hvar == 1 % go to mann whitney u
        [p, h] = ranksum(spec(allSpecs.Group == 'CaCl'), spec(allSpecs.Group == 'Sham'));
    end
        
    eval(['Stats.' specs{indspec} ' = p;']);
    
end


% delay frames/categorical




end