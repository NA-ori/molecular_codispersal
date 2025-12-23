function [current_chemical_counts] = outflow(p, diff_mask, current_chemical_counts, I)

    % For the moment, this is a placeholder. A more accurate dispersal
    % function will be added once the basic concept is tested.

    arguments (Input)
        p;
        diff_mask;
        current_chemical_counts;
        I;
        %coordinate_list;
    end

    diff_mask = I.out_mask;

    % Calculate where everything goes
    for coord = 1:size(current_chemical_counts,2)
        for species = 1:size(current_chemical_counts{1,coord},2)
            if (diff_mask(species) == 1) && (current_chemical_counts{1,coord}{2,species} > 0)

                if contains(current_chemical_counts{1,coord}{1,species}, "prop")
                    leave_prob = (p.out_rate / (p.flow_rate*6))/sqrt(max(p.subcycles_per_ring));
                    stay_prob = 1 - leave_prob;       % THIS IS CHEATING, CORRECT THIS             
                    sample_array = randsample([0,-1], ceil((current_chemical_counts{1,coord}{2,species})*p.diffusion_divisor), true, [stay_prob, leave_prob]);                    
                    
                else
                    leave_prob = p.out_rate / (p.flow_rate*6);
                    stay_prob = 1 - leave_prob;
                    sample_array = randsample([0,-1], ceil((current_chemical_counts{1,coord}{2,species})*p.diffusion_divisor), true, [stay_prob, leave_prob]);
                end

                [counts, ordering] = groupcounts(sample_array');
                for d = 1:length(ordering)
                    if ordering(d) == -1
                        current_chemical_counts{coord}{2,species} = current_chemical_counts{coord}{2,species} - counts(d);
                    end
                end

            end
        end
    end

end