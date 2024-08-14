
function [Stats] = CorrMatricesStats
Groups = {'CaCl', 'Sham'};
SaveFolder = '/media/mbakker/GDrive/P2/GCaMP/CombinedCorrMatrices/';
datatypes = {'GCaMP', 'HbO', 'HbR'};
Acquisitions = {'A1', 'A2', 'A3'};
allrois = {'VR', 'SR', 'MR', 'RR', 'VL', 'SL', 'ML', 'RL'};
[allrois, sortindex] = sort(allrois); % alfabetisch

for inddata = 1:size(datatypes, 2)
    for indAcq = 1:size(Acquisitions, 2)
        Acq = Acquisitions{indAcq};
        datatype = datatypes{inddata};

        % disp([Acq ' ' datatype])

        % Get data
        for indgroup = 1:size(Groups, 2)
            group = Groups{indgroup};
            load([SaveFolder datatype filesep group '_' Acq '_corrmatrices.mat'], 'allcorrmatrices')
            
            % z transform fisher
            allcorrmatrices = atanh(allcorrmatrices);

            allcorrmatrices(allcorrmatrices == inf) = NaN; %mouse x roi
            eval([group ' = allcorrmatrices;'])
        end

        %% Check assumptions
        norm_ad = NaN(64, 1);
        norm_ks = NaN(64, 1);
        var = NaN(64, 1);

        for indroi = 1:64
            % if it's roi with itself, take next
            if sum(~isnan([CaCl(:,indroi); Sham(:,indroi)])) < 4 
                continue
            end
            
            % check normal distribution (anderson darling & kolmogorov smirnov)
            % [h1, ~] = adtest([CaCl(:,indroi); Sham(:,indroi)]); %note: gives warning that p is larger than largest tabulated value...?
            [h2, ~] = kstest([CaCl(:,indroi); Sham(:,indroi)]);
            % norm_ad(indroi) = h1;
            norm_ks(indroi) = h2;

            % check homogeneity of variance (bartlett)
            [p] = vartestn([CaCl(:,indroi); Sham(:,indroi)], 'Display', 'off');
            if p < 0.05
                var(indroi) = 1;
            else
                var(indroi) = 0;
            end
        end

        % Note: for A1, fluo: all adtests said it was a normal
        % distribution, all ks tests said it was not... 
        % Want to use parametric test based on ks results. Kruskal wallis
        % is best for more than 2 groups. Only have 2 groups here, so will
        % use Mann whitney u test. 
        % Assumptions Mann Whitney U:
        % - sample drawn from population is random
        % - independence within samples & mutual independence
        % - ordinal measurement scale
        % Mann whitney u is same as wilcoxon rank sum in matlab
        clear norm_ad norm_ks p h1 h2 indroi var

        % % start t-test
        % [h, p, ci, stats] = ttest2(CaCl(:,indroi), Sham(:,indroi));
        % pvalues(indroi) = p;

        %% start wilcoxon rank sum test
        pvalues = NaN(64, 1);
        hvalues = NaN(64, 1);
        for indroi = 1:64
            if sum(~isnan([CaCl(:,indroi); Sham(:,indroi)])) < 4
                continue
            end
            [p, h] = ranksum(CaCl(:,indroi), Sham(:,indroi));
            pvalues(indroi) = p;
            hvalues(indroi) = h;
        end

        if any(hvalues, 1)
            disp(['p SIGN! ' Acq ' ' datatype])
        else
            disp(['Nothing sign... ' Acq ' ' datatype])
        end

        % eval(['pvalues_' Acq '_' datatype ' = pvalues;'])
        % eval(['hvalues_' Acq '_' datatype ' = hvalues;'])
        
        %% Do FDR
        pvalues = reshape(pvalues, 8, 8);
        pvalues = pvalues(sortindex, sortindex); %alfabetisch
        pvalues = tril(pvalues, -1);
        pvalues = reshape(pvalues, 64, 1);

        pvalues28 = [];
        for ind = 1:size(pvalues,1)
            if pvalues(ind) ~= 0
                pvalues28 = [pvalues28; pvalues(ind)];
            end
        end
        clear ind

        qvalues = mafdr(reshape(pvalues28, [], 1),'BHFDR', 'true');

        q = tril(ones(8), -1);
        q = reshape(q, 64, 1);
        ind2 = 1;
        for ind = 1:64
            if q(ind) ~= 0
                q(ind) = qvalues(ind2);
                ind2 = ind2 + 1;
            end
        end
        q = reshape(q, 8, 8);
        pvalues = reshape(pvalues, 8,8);
        pqcombi = pvalues' + q;
        q(q==0) = NaN;
        pvalues(pvalues==0) = NaN;

        if any(q<0.05)
            disp('HHHHHHHAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA HOLY SHIT')
        end

        % Keep results
        % eval(['pvalues_' Acq '_' datatype ' = pvalues;'])
        % eval(['qvalues_' Acq '_' datatype ' = qvalues;'])
        % eval(['Stats.pqcombi_' Acq '_' datatype ' = pqcombi;'])
        % writematrix(pqcombi, [SaveFolder 'pqvalues.xlsx'], 'Sheet', [datatype '_' Acq]);
        
        pqcombi(pqcombi == 0) = 1;
        pqcombi = [['ROIs'; allrois'], [allrois; num2cell(pqcombi)]];
        writecell(pqcombi, [SaveFolder 'pqvalues.xlsx'], 'Sheet', [datatype '_' Acq]);
        
    end

end

end     