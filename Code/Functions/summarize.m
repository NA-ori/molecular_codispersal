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
            row_sum = row_sum + cell2mat(O.concentration_tracker{coord}(row,:));
        end
        global_species_counts(row,:) = num2cell(row_sum);
    end

    O.global_species_counts = global_species_counts;

    % Return relative occupancy of cycle 1

    cycle_of_interest = 1;
    total_sum = 0;

    member_species = I.member_species_record;
    for cycle = 1:size(I.member_species_record,2)
        for spec = 1:size(I.member_species_record{cycle},2)
            member_species{cycle}(spec) = convertStringsToChars(member_species{cycle}(spec) + p.seed_state(cycle));
        end
    end

    for cycle = 1:length(I.member_species_record)
        indeces = matches(I.base_species, member_species{cycle});
        cycle_sum = sum(cell2mat(global_species_counts(end,indeces)));
        total_sum = total_sum + cycle_sum;
        if cycle == cycle_of_interest
            interest_sum = cycle_sum;
        end
    end

    relative_concentration = interest_sum / total_sum;
    
    O.relative_concentration = relative_concentration;

end