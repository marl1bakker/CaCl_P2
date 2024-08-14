%% compact letter display
% Needs a table with a column ROI_1, ROI_2 and pValue. The rois that are
% compared with each other need to be double: if ROI_1 is A and ROI_2 is B,
% you also need ROI_1 B and ROI_2 A. If you have 8 roi, youll have 56 rows.
% The output is the group numbers and the letters according to CLD. If
% ROIs have one or more letters in common, they do not differ
% significantly.

% Note: if a group has letters that are not in alfabetical order, like acb,
% it is likely that the last letter is not necessary. for example, if you
% have a group ab and a group bc, and the unalfabetical letter is not sign.
% different from either of these groups, then ac is enough to indicate
% that, and the b would be unessecary. It's not wrong, but it's not
% optimal. So if you have a letter pair out of alfabetical order, check if
% the last one(s) are really necessary. 


function groups = cld(c)

gn=['a';'b';'c';'d';'e';'f';'g';'h';'i';'j';'k';'l';'m';'n';'o';'p';'q';'r';'t';'u';'v';'w';'x';'y';'z'];
groups = table;
groups.nr = unique([c.ROI_1; c.ROI_2]);
groups.letter = repmat({'0'}, size(groups,1),1);
gn_ind = 0;

for roi1 = 1:size(groups,1)
    % get all pvals comparisons of  a group
    gr1 = c(c.ROI_1 == roi1,:);
    nonsig_ind = find(gr1.pValue>0.05);

    % give group letter if it doesnt have one yet
    if matches((groups.letter(roi1)), '0')
        gn_ind = gn_ind+1;
        groups.letter(roi1) = {gn(gn_ind)};
    end

    if size(nonsig_ind) == 0
        % disp('Group sign. diff from all other groups')
        continue

    else
        % get non-significant pairings
        nonsiggroups = gr1.ROI_2(nonsig_ind); %all groups

        for ind2 = 1:size(nonsiggroups)
            roi2 = nonsiggroups(ind2);
            gr2 = c(c.ROI_1 == roi2,:);
            % ind_group_sigs = c(c.ROI_1 == ind,:);
            % othergroups = groups(groups.nr ~= ns_group_nr,:);
            % othergroups = othergroups(othergroups.nr ~= ind,:);

            % non significant pairing already has a letter in common, skip
            if contains(groups(groups.nr == roi1,:).letter{:}, groups(groups.nr == roi2,:).letter{:}) ||...
                    contains(groups(groups.nr == roi2,:).letter{:}, groups(groups.nr == roi1,:).letter{:}) ...
                    && ~matches(groups(groups.nr == roi1,:).letter{:}, '0') && ~matches(groups(groups.nr == roi2,:).letter{:}, '0')
                % disp('already letter in common')
                continue
            end

            % if it fits into a group:
            done = 0;
            for prev_let_ind = 1:gn_ind % check all already existing groups
                roi_in_group_ind = contains(groups.letter, gn(prev_let_ind));
                roi_in_group = groups.nr(roi_in_group_ind);
                roi_in_group(roi_in_group == roi2) = [];
                groupmatch = 1;

                % check if roi1 and roi2 match all prev groups
                for ind_3 = 1:size(roi_in_group,1)
                    roi3 = roi_in_group(ind_3);
                    p1 = gr1(gr1.ROI_2 == roi3,:).pValue;
                    p2 = gr2(gr2.ROI_2 == roi3,:).pValue;
                    if ~isempty(p1) && p1<0.05
                        groupmatch = 0;
                        break
                    end
                    if ~isempty(p2) && p2<0.05
                        groupmatch = 0;
                        break
                    end
                end

                if groupmatch == 1
                % if all(ns_group_sigs.pValue(roi_in_group_ind) > 0.05)

                % if all(ns_group_sigs.pValue(contains(othergroups.letter, gn(gn_ind_2))) > 0.05) && ...
                %     groups(groups.nr == gn(gn_ind_2),:).letter
                    % ind_group_sigs(ind_group_sigs.ROI_2 == ns_group_nr,:).pValue > 0.05
                    % all(ind_group_sigs.pValue(contains(othergroups.letter, gn(gn_ind_2))) > 0.05) &&...
                    

                    % fits: if it doesnt have a letter yet, make one
                    if matches(groups(groups.nr == roi2,:).letter, '0')
                        groups(groups.nr == roi2,:).letter = {gn(prev_let_ind)};

                        % fits: otherwise add a letter to one or both groups
                    else
                        if ~contains(groups(groups.nr == roi2,:).letter{:}, gn(prev_let_ind))
                            groups(groups.nr == roi2,:).letter{:} = [groups(groups.nr == roi2,:).letter{:}, gn(prev_let_ind)];
                        end
                        if ~contains(groups(groups.nr == roi1,:).letter{:}, gn(prev_let_ind))
                            groups(groups.nr == roi1,:).letter{:} = [groups(groups.nr == roi1,:).letter{:}, gn(prev_let_ind)];
                        end
                    end

                    done = done+1;
                    break
                end
            end

            % if it doesnt fit with the rest of the group:
            if done == 0
                gn_ind = gn_ind+1; % get a new letter
                % groups(groups.nr == ind,:).letter{:} = [groups(groups.nr == ind,:).letter{:}, gn(gn_ind)];
                if matches(groups(groups.nr == roi2,:).letter{:}, '0')
                    groups(groups.nr == roi2,:).letter = {gn(gn_ind)};
                else
                    groups(groups.nr == roi2,:).letter{:} = [groups(groups.nr == roi2,:).letter{:}, gn(gn_ind)];
                end

                if matches(groups(groups.nr == roi1,:).letter{:}, '0')
                    groups(groups.nr == roi1,:).letter = {gn(gn_ind)};
                else
                    groups(groups.nr == roi1,:).letter{:} = [groups(groups.nr == roi1,:).letter{:}, gn(gn_ind)];
                end
            end


        end
    end
end
end

