function [current_chemical_counts, current_reaction_propensities] = disperse(p, update_mask, diff_mask, current_chemical_counts, current_reaction_propensities, coordinate_list, prob_cloud, prop_prob_cloud)

    % For the moment, this is a placeholder. A more accurate dispersal
    % function will be added once the basic concept is tested.

    arguments (Input)
        p;
        update_mask;
        diff_mask;
        current_chemical_counts;
        current_reaction_propensities;
        coordinate_list;
        prob_cloud;
        prop_prob_cloud;
    end

    % Calculate where everything goes
    for coord = 1:size(current_chemical_counts,2)
        for species = 1:size(current_chemical_counts{1,coord},2)
            if (diff_mask(species) == 1) && (current_chemical_counts{1,coord}{2,species} > 0)
                % Must update prop_prob_cloud so it accounts for potential
                % multiple rings with different numbers of subcycles
                if contains(current_chemical_counts{1,coord}{1,species}, "prop")
                    sample_array = randsample(1:size(coordinate_list,2), current_chemical_counts{1,coord}{2,species}, true, cell2mat(prop_prob_cloud{coord}));
                else
                    sample_array = randsample(1:size(coordinate_list,2), current_chemical_counts{1,coord}{2,species}, true, cell2mat(prob_cloud{coord}));
                end
                [counts, ordering] = groupcounts(sample_array');
                for d = 1:length(ordering)
                    update_mask{ordering(d)}{2,species} = update_mask{ordering(d)}{2,species} + counts(d);
                    update_mask{coord}{2,species} = update_mask{coord}{2,species} - counts(d);
                end
            end
        end
    end

    % Move things around
    for coord = 1:size(current_chemical_counts,2)
        for species = 1:size(current_chemical_counts{1,coord},2)
            current_chemical_counts{1,coord}{2,species} = current_chemical_counts{1,coord}{2,species} + update_mask{1,coord}{2,species};
        end
    end

    % Update reaction propensities
    for coord = 1:size(coordinate_list, 2)
        all_reactions = 1:size(current_reaction_propensities{1,coord},2);
        current_reaction_propensities = update_propensities(all_reactions, coord, reactions, current_chemical_counts, current_reaction_propensities);
    end

end