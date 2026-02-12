function [O] = summarize(p, I, O)

    % Takes raw output files and compiles them into the most useful bits as
    % a single table good for saving to a csv or whatnot

    arguments (Input)
        p;
        I;
        O;
    end

    % Return sums of species counts across all pixels
    global_species_counts = cellstr(I.base_species);
    global_species_counts(2:length(O.time)+1,:) = {0};
    
    for row = 2:size(global_species_counts,1)
        row_sum = cell2mat(global_species_counts(row,:));
        for coord = 1:size(O.concentration_tracker,2)
            row_sum = row_sum + cell2mat(O.concentration_tracker{1,coord}(row,:));
        end
        global_species_counts(row,:) = num2cell(row_sum);
    end

    O.global_species_counts = global_species_counts;

    % Return relative occupancy of cycle 1

    cycle_of_interest = 1;
    total_sum = 0;
    site_occupation = zeros(1, p.rings);

    member_species = I.member_species_record;
    all_species = I.ring_list;
    for cycle = 1:size(I.member_species_record,2)
        for spec = 1:size(I.member_species_record{cycle},2)
            member_species{cycle}(spec) = convertStringsToChars(member_species{cycle}(spec) + "_ad");
        end
        for spec = 1:size(I.ring_list{cycle},2)
            all_species{cycle}(spec) = convertStringsToChars(all_species{cycle}(spec) + "_ad");
        end
    end


    for cycle = 1:length(I.member_species_record)
        % Calculate relative concentrations
        indeces = matches(I.base_species, member_species{cycle});
        cycle_sum = sum(cell2mat(global_species_counts(end,indeces)));
        total_sum = total_sum + cycle_sum;
        if cycle == cycle_of_interest
            interest_sum = cycle_sum;
        end

        % Calculate percent site occupation
        all_indeces = matches(I.base_species, all_species{cycle});
        all_species_cycle_sum = sum(cell2mat(global_species_counts(end,all_indeces)));
        site_occupation(cycle) = all_species_cycle_sum / (p.site_concentration*length(I.coordinate_list));
    end

    relative_concentration = interest_sum / total_sum;
    
    O.relative_concentration = relative_concentration;
    O.site_occupation = site_occupation;

end