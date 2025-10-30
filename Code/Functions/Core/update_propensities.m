function [updated_propensities] = update_propensities(reactions_to_update, coord, reactions, current_chemical_counts, current_reaction_propensities)

    arguments (Input)
        reactions_to_update;
        coord;
        reactions;
        current_chemical_counts;
        current_reaction_propensities;
    end

    for i = 1:length(reactions_to_update)
        reaction = reactions_to_update(i);
        h_i = 1;    % The propensity
        for reactant = 1:length(reactions{reactions}{2})
            reactant_count = current_chemical_counts(reactions{reaction}{2}{reactant});
            reactant_stoichiometry = reactions{reaction}{3}{reactant};

            if reactant_count < reactant_stoichiometry
                h_i = 0;
                break;
            else
                h_i = h_i * nchoosek(reactant_count, reactant_stoichiometry);    % how many combinations of n items (reactant_stoichiometry) can be taken from a set of size k (reactant_count)
            end
        end
        rate_constant = full_network{reaction}{1};
        current_reaction_propensities{coord}{reaction} = h_i * rate_constant;  % the propensity of a particular reaction
    end
    updated_propensities = current_reaction_propensities;
end