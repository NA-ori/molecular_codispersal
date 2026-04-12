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
        reaction = reactions_to_update(i);
        if I.tags(reaction) == "prop_form"

            h_i = 1;

            reactant_counts = current_chemical_counts(absolute_coord,I.reactions(reaction,:)<0);
            reactant_stoichiometry = -(I.reactions(reaction,(I.reactions(reaction,:)<0)));

            [reactant_counts,min_indeces] = mink(reactant_counts,2);
            reactant_stoichiometry = reactant_stoichiometry(min_indeces);

            for reactant = 1:length(reactant_counts)
                if reactant_counts(reactant) < reactant_stoichiometry(reactant)
                    h_i = 0;
                    break;
                else
                    h_i = h_i * nchoosek(reactant_counts(reactant), reactant_stoichiometry(reactant));
                end
            end
            current_reaction_propensities(propensity_indices(i)) = h_i * I.rate_constants(reaction);  % the propensity of a particular reaction

        else

            h_i = 1;    % The propensity

            % reactant_stoichiometry = I.reactions(reaction,:);
            % reactant_counts = current_chemical_counts(absolute_coord,reactant_stoichiometry<0);
            % reactant_stoichiometry = -reactant_stoichiometry(reactant_stoichiometry<0);

            reactant_counts = current_chemical_counts(absolute_coord,I.reactions(reaction,:)<0);
            reactant_stoichiometry = -(I.reactions(reaction,(I.reactions(reaction,:)<0)));
            
            for reactant = 1:length(reactant_counts)

                if reactant_counts(reactant) < reactant_stoichiometry(reactant)
                    h_i = 0;
                    break;
                else
                    h_i = h_i * nchoosek(reactant_counts(reactant), reactant_stoichiometry(reactant));    % how many combinations of n items (reactant_stoichiometry) can be taken from a set of size k (reactant_count)
                end
            end
            current_reaction_propensities(propensity_indices(i)) = h_i * I.rate_constants(reaction);  % the propensity of a particular reaction

        end
    end
    clear h_i; clear reactant_counts; clear reactant_stoichiometry; clear reaction; clear reactant;
    updated_propensities = current_reaction_propensities;

end