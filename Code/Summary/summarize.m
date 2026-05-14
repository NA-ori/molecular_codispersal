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
            row_sum = row_sum + O.concentration_tracker{1,coord}(row-1,:);
        end
        global_species_counts(row,:) = num2cell(row_sum);
    end

    time_col = ["t"; double(O.time)];
    O.global_species_counts = [time_col, global_species_counts];


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

    O.member_species = member_species;
    O.all_species = all_species;

    % Calculate extinction

    O.adsorbed_A_extinct = 0;
    O.adsorbed_NA_extinct = 0;
    O.total_A_extinct = 0;
    O.total_NA_extinct = 0;

    O.adsorbed_none_extinct = 0;
    O.adsorbed_all_extinct = 0;
    O.total_none_extinct = 0;
    O.total_all_extinct = 0;

    A_count_ad = 0;
    NA_count_ad = 0;
    A_count_all = 0;
    NA_count_all = 0;

    for i = 1:length(I.base_species)
        d = i+1;

        if ( (contains(I.base_species(i), "_r1")) && (contains(I.base_species(i), "_ad")) ) || (contains(I.base_species(i), "prop"))
            A_count_ad = A_count_ad + str2double(O.global_species_counts(end,d));
            A_count_all = A_count_all + str2double(O.global_species_counts(end,d));
        end
        if ( (contains(I.base_species(i), "_r1")) && (contains(I.base_species(i), "_diff")) )
            A_count_all = A_count_all + str2double(O.global_species_counts(end,d));
        end        

        if ( (contains(I.base_species(i), "_r2")) && (contains(I.base_species(i), "_ad")) )
            NA_count_ad = NA_count_ad + str2double(O.global_species_counts(end,d));
            NA_count_all = NA_count_all + str2double(O.global_species_counts(end,d));
        end
        if ( (contains(I.base_species(i), "_r2")) && (contains(I.base_species(i), "_diff")) )
            NA_count_all = NA_count_all + str2double(O.global_species_counts(end,d));
        end 

    end

    if A_count_all == 0 && NA_count_all == 0
        O.total_all_extinct = 1;
    elseif A_count_all == 0 && NA_count_all > 0
        O.total_A_extinct = 1;
    elseif A_count_all > 0 && NA_count_all == 0
        O.total_NA_extinct = 1;
    elseif A_count_all > 0 && NA_count_all > 0
        O.total_none_extinct = 1;
    end

    if A_count_ad == 0 && NA_count_ad == 0
        O.adsorbed_all_extinct = 1;
    elseif A_count_ad == 0 && NA_count_ad > 0
        O.adsorbed_A_extinct = 1;
    elseif A_count_ad > 0 && NA_count_ad == 0
        O.adsorbed_NA_extinct = 1;
    elseif A_count_ad > 0 && NA_count_ad > 0
        O.adsorbed_none_extinct = 1;
    end


end