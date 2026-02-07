function [updated_propensities] = linear_update_propensities(absolute_coord,reactions_to_update,current_chemical_counts,current_reaction_propensities,propensity_indices,I)

    arguments (Input)
        absolute_coord;
        reactions_to_update;
        current_chemical_counts;
        current_reaction_propensities;
        propensity_indices;
        I;
    end


    for i = 1:length(reactions_to_update)
        %reactions_to_update
        reaction = reactions_to_update(i);
        h_i = 1;    % The propensity
        for reactant = 1:length(I.reactions{reaction}{2})
            index = I.species_index_map{reaction}{reactant};
            reactant_count = current_chemical_counts{absolute_coord}{2,index};
            reactant_stoichiometry = I.reactions{reaction}{3}{reactant};

            if reactant_count < reactant_stoichiometry
                h_i = 0;
                break;
            else
                h_i = h_i * nchoosek(reactant_count, reactant_stoichiometry);    % how many combinations of n items (reactant_stoichiometry) can be taken from a set of size k (reactant_count)
            end
        end
        rate_constant = I.reactions{reaction}{1};
        propensity_index = propensity_indices(i);
        current_reaction_propensities(propensity_index) = h_i * rate_constant;  % the propensity of a particular reaction
    end
    clear rate_constant; clear propensity_index; clear h_i; clear index; clear reactant_count; clear reactant_stoichiometry; clear reaction;
    updated_propensities = current_reaction_propensities;
end