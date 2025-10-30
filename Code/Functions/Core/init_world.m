function [species_counts, coordinate_list] = init_world(p, reactions, ring_list, base_species, member_species_record)
    % Initialize the dimensions of the world and all the data structures
    % that will hold information

    arguments (Input)
        p;
        reactions;
        ring_list;
        base_species;
        member_species_record;
    end

    % Initialize coordinate system
    % This is a placeholder.
    coordinate_list = 1:p.dim;

    % Make separate data structures to keep track of species in every
    % coordinate
    species_counts = cell(1, p.dim);
    for coord = 1:size(species_counts, 2)
        species_counts{1, coord} = cellstr(base_species);
        species_counts{1, coord}(2,:) = {0};

        % Set starting concentrations of sites and food

        indeces = strcmp(base_species, 'F');
        species_counts{1, coord}(2,indeces) = {p.food_concentration};
        indeces = strcmp(base_species, 'site');
        species_counts{1, coord}(2,indeces) = {p.site_concentration};

        % Prepare seeds
        for ring = 1:length(ring_list)
            temp_member_species_record = member_species_record;
            for seed = 1:size(member_species_record{1,ring},2)
                temp_member_species_record{1,ring}{seed} = convertStringsToChars(temp_member_species_record{1,ring}{seed} + p.seed_state(ring));
            end
            if coord == p.seed_locations(ring)
                indeces = matches(base_species, temp_member_species_record{ring});
                species_counts{1, coord}(2,indeces) = {p.seed_concentration(ring)};
            end
        end

    end
    
end