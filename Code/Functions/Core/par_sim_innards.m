function [I, current_reaction_propensities, current_chemical_counts] = par_sim_innards(p, I, current_reaction_propensities, current_chemical_counts, coord)

    % Put parallel computation of catalyzed sites in here to try to fix
    % problems.

    arguments (Input)
        p;
        I;
        current_reaction_propensities;
        current_chemical_counts;
        coord;
    end

    error = false;

    absolute_coord = I.catalyzed_sites_mask(coord);
    local_time = 0;

    while local_time <= p.disperse_frequency    % Lesser loop

        r_1 = rand();
        r_2 = rand();

        a_0 = sum(cell2mat(current_reaction_propensities{coord}));

        tau = (1/a_0) * log(1/r_1);

        % now determine which reaction will occur
        mu = 0;
        summ = 0;
        target = r_2*a_0;   % target is some random percentage of the total propensities a_0
        for reaction = 1:size(I.reactions, 1)
            summ = summ + current_reaction_propensities{coord}{reaction};
            if summ > target
                mu = reaction;  % the reaction that happens to pass the target value first will occur.
                break
            end
        end

        % Change concentrations and update propensities
        for reactant = 1:size(I.reactions{mu}{2}, 2)
            % subtract the number of particles that react from the current counts
            index = strcmp(current_chemical_counts{absolute_coord}(1,:), I.reactions{mu}{2}{reactant});
            current_chemical_counts{absolute_coord}{2,index} = current_chemical_counts{absolute_coord}{2,index} - I.reactions{mu}{3}{reactant};
            if current_chemical_counts{absolute_coord}{2,index} < 0
                fprintf("Negative concentration detected >:3\n");
                fprintf("Offending reaction: " + num2str(mu) + "\n");
                fprintf("Offending propensity: " + current_reaction_propensities{absolute_coord}(mu) + "\n");
                error = true;
            end
        end
        for product = 1:size(I.reactions{mu}{4}, 2)
            % Add number of products formed to current counts
            index = strcmp(current_chemical_counts{absolute_coord}(1,:), I.reactions{mu}{4}{product});
            current_chemical_counts{absolute_coord}{2,index} = current_chemical_counts{absolute_coord}{2,index} + I.reactions{mu}{5}{product};
        end
        if error == true
            break
        end

        % Update the propensities
        reactions_to_update = I.reaction_update_key{mu};
        current_reaction_propensities = update_propensities(reactions_to_update,coord,I.reactions,current_chemical_counts,current_reaction_propensities,I);

        local_time = local_time + tau;
        if local_time > p.disperse_frequency
            local_time = p.disperse_frequency;
            break
        end
    end
end