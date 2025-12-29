function [I] = init_world(p, I)
    % Initialize the dimensions of the world and all the data structures
    % that will hold concentration information

    % Also generate the probability distribution for the given coords

    arguments (Input)
        p;
        I;
    end

    % Pop variables out of I
    reactions = I.reactions;
    ring_list = I.ring_list;
    base_species = I.base_species;
    member_species_record = I.member_species_record;


    % Initialize coordinate system
    coordinate_list = {};
    I.all_coordinates = {};

    if p.shape == "parallelogram"

        D = p.separation_distance;
        s = p.sites;
        R = s*D;
    
        for x = 0:D:(D*(s-1))
            for y = 0:D:(D*(s-1))
                coordinate_list{1,end+1} = {};
                coordinate_list{1,end}{1,1} = x;
                coordinate_list{1,end}{1,2} = y;
                coordinate_list{1,end}{1,3} = (0-x-y);
            end
        end

    elseif p.shape == "hex"

        D = p.separation_distance;
        s = p.sites;
        R = s*D + floor(D/2);   I.R = R;
        omega = 3*R*(R+1) + 1;  I.omega = omega;

        I.all_coordinates = cell(1,omega);
        I.empty_coord_mask = [];
        coordinate_list = {};
        catalyzed_sites_mask = [];
        origin = 0;

        % Set up pixels
        loops = 1;
        for x = -R:R
            for y = -R:R
                for z = -R:R
                    if (x+y+z) == 0
                        I.all_coordinates{1,loops}{1} = x; I.all_coordinates{1,loops}{2} = y; I.all_coordinates{1,loops}{3} = z;
                        if mod(abs(x),D) == 0 && mod(abs(y),D) == 0 && mod(abs(z),D) == 0
                            coordinate_list{1,end+1}{1} = x; coordinate_list{1,end}{2} = y; coordinate_list{1,end}{3} = z;
                            catalyzed_sites_mask = [catalyzed_sites_mask, loops];
                        else
                            I.empty_coord_mask = [I.empty_coord_mask, loops];
                        end
                        if x == 0 && y == 0 && z == 0
                            origin = loops;
                        end
                        loops = loops + 1;
                    end
                end
            end
        end
        if length(catalyzed_sites_mask) ~= length(coordinate_list)
            fprintf("PROBLEM!!!! >:3")
        end
       
    end

    I.R = R;
    I.map_name = p.map_prefix + "probmap_R" + num2str(R) + ".mat";
    I.prop_map_name = p.map_prefix + "probmap_R" + num2str(R) + "_" + num2str(max(p.subcycles_per_ring) + ".mat");

    if isfile(I.map_name)
        fprintf("Found " + I.map_name + "!\n");
        fprintf("Loading precalculated probability map... >w<\n");
        load(I.map_name, "prob_cloud");
    else
        fprintf("Couldn't find precalculated probability map, making one from scratch >w<\n");
        fprintf("Making " + I.map_name + "...\n")
        precalculate_occupation_p(p,I);
        load(I.map_name, "prob_cloud");
    end

    if isfile(I.prop_map_name)
        fprintf("Found " + I.prop_map_name + "!\n");
        fprintf("Loading precalculated prop probability map... >w<\n");
        load(I.prop_map_name, "prop_prob_cloud");
    else
        fprintf("Couldn't find precalculated prop probability map, making one from scratch >w<\n");
        fprintf("Making " + I.prop_map_name + "...\n")
        precalculate_occupation_p(p,I,prop=true);
        load(I.prop_map_name, "prop_prob_cloud");
    end


    % Make separate data structures to keep track of species in every
    % coordinate
    species_counts = cell(1, size(I.all_coordinates,2));
    for coord = 1:size(I.all_coordinates,2)
        species_counts{1, coord} = cellstr(base_species);
        species_counts{1, coord}(2,:) = {0};
    end


    for c = 1:length(catalyzed_sites_mask)
        % Set starting concentrations of sites and food
        coord = catalyzed_sites_mask(c);
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
            %if coord == p.seed_locations(ring)
            if coord == origin
                indeces = matches(base_species, temp_member_species_record{ring});
                species_counts{1, coord}(2,indeces) = {p.seed_concentration(ring)};
            end
        end

    end

    % Map species indeces to reactions for faster propensity updating
    species_index_map = {};

   for i = 1:length(reactions)
        reaction = i;
        species_index_map{1,end+1} = {};
        for reactant = 1:length(reactions{reaction}{2})
            index = strcmp(species_counts{1}(1,:), reactions{reaction}{2}{reactant});
            species_index_map{end}{end+1} = index;
        end
   end

    % Stick outputs into I

    I.species_counts = species_counts;
    I.coordinate_list = coordinate_list;
    I.prob_cloud = prob_cloud;
    I.prop_prob_cloud = prop_prob_cloud;
    I.catalyzed_sites_mask = catalyzed_sites_mask;
    I.species_index_map = species_index_map;
    I.origin = origin;
    
end