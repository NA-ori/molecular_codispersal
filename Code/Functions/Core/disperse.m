function [] = disperse(p, current_chemical_counts, coordinate_list, prob_cloud)

    % For the moment, this is a placeholder. A more accurate dispersal
    % function will be added once the basic concept is tested.

    arguments (Input)
        p;
        current_chemical_counts;
        coordinate_list;
        prob_cloud;
    end

    % Move chemicals around
    update_mask = cell(1,length(current_chemical_counts));

    for species = 1:current_chemical_counts{1,:}
        for particlecount = 1:current_chemical_counts{1,species}{1,particle}
            
        end
    end


    % Update reaction propensities

end