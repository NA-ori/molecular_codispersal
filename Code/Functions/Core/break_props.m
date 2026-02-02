function current_chemical_counts = break_props(p, I, current_chemical_counts)

    % Fake some Gillespie steps so propagules can break between catalyzed
    % pixels

    arguments (Input)
        p;
        I;
        current_chemical_counts;
    end

    index = contains(I.base_species, "prop_");
    food_surrogate = p.in_rate;
    mu = I.prop_break_reaction;

    for coord = 1:length(I.empty_coord_mask)
        absolute_coord = I.empty_coord_mask(coord);
        prop_count = current_chemical_counts{absolute_coord}{2, index};

        if prop_count > 0

            local_time = 0;
            while local_time <= p.disperse_frequency

                r_1 = rand();
                r_2 = rand();
    
                a_0 = food_surrogate + (prop_count*p.prop_break_rate);

                tau = (1/a_0) * log(1/r_1);

                % Determine whether a prop will break
                target = r_2*a_0;
                if target > food_surrogate && target <= a_0
                    for reactant = 1:size(I.reactions{mu}{2}, 2)
                        % subtract the number of particles that react from the current counts
                        index = strcmp(current_chemical_counts{absolute_coord}(1,:), I.reactions{mu}{2}{reactant});
                        current_chemical_counts{absolute_coord}{2,index} = current_chemical_counts{absolute_coord}{2,index} - I.reactions{mu}{3}{reactant};
                    end
                    for product = 1:size(I.reactions{mu}{4}, 2)
                        % Add number of products formed to current counts
                        index = strcmp(current_chemical_counts{absolute_coord}(1,:), I.reactions{mu}{4}{product});
                        current_chemical_counts{absolute_coord}{2,index} = current_chemical_counts{absolute_coord}{2,index} + I.reactions{mu}{5}{product};
                    end
                    if p.introspection == true
                         fprintf("A propagule burst in " + absolute_coord + " ^w^\n");
                    end
                end

                local_time = local_time + tau;
                if local_time > p.disperse_frequency
                    break
                end

            end
        end
    end

end