function [species_counts, coordinate_list] = init_world(p, reactions, ring_list, base_species, member_species_record)
    % Initialize the dimensions of the world and all the data structures
    % that will hold concentration information

    % Also generate the probability distribution for the given coords

    arguments (Input)
        p;
        reactions;
        ring_list;
        base_species;
        member_species_record;
    end

    % Initialize coordinate system

    D = p.separation_distance;
    s = p.sites;
    R = s*D;

    coordinate_list = {};

    for x = 0:D:(D*(s-1))
        for y = 0:D:(D*(s-1))
            coordinate_list{1,end+1} = {};
            coordinate_list{1,end}{1,1} = x;
            coordinate_list{1,end}{1,2} = y;
            coordinate_list{1,end}{1,3} = (0-x-y);
        end
    end
    
    % Calculate probability distribution around each point
    prob_cloud = cell(1,size(coordinate_list,2));
    for n0 = 1:size(coordinate_list,2)
        for n1 = 1:size(coordinate_list,2)

        end
    end


    % Make separate data structures to keep track of species in every
    % coordinate
    species_counts = cell(1, size(coordinate_list,2));
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