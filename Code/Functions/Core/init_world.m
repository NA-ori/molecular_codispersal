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
        R = s*D;

        % Make a list of catalyzed sites
        used_sites = ["000"];
        coordinate_list{1,end+1} = {};
        coordinate_list{1,end}{1,1} = 0;
        coordinate_list{1,end}{1,2} = 0;
        coordinate_list{1,end}{1,3} = 0;
        for ring = 1:s
            current_rings = size(coordinate_list,2);
            for coord = 1:current_rings
                x = coordinate_list{1,coord}{1}; y = coordinate_list{1,coord}{2}; z = coordinate_list{1,coord}{3};
                friends = neighbors(x,y,z,D);
                for friend = 1:size(friends,2)
                    new_coord = strcat(num2str(friends{1,friend}{1}), num2str(friends{1,friend}{2}), num2str(friends{1,friend}{3}));
                    if ~any(strcmp(used_sites, new_coord))
                        coordinate_list{1,end+1} = friends{1,friend};
                        used_sites = [used_sites, new_coord];
                    end
                end
            end
        end

        % Make a list of all sites, including empty space
        second_used_sites = ["000"];
        I.all_coordinates{1,end+1} = {};
        I.all_coordinates{1,end}{1,1} = 0;
        I.all_coordinates{1,end}{1,2} = 0;
        I.all_coordinates{1,end}{1,3} = 0;
        for ring = 1:R
            current_rings = size(I.all_coordinates,2);
            for coord = 1:current_rings
                x = I.all_coordinates{1,coord}{1}; y = I.all_coordinates{1,coord}{2}; z = I.all_coordinates{1,coord}{3};
                friends = neighbors(x,y,z,1);
                for friend = 1:size(friends,2)
                    new_coord = strcat(num2str(friends{1,friend}{1}), num2str(friends{1,friend}{2}), num2str(friends{1,friend}{3}));
                    if ~any(strcmp(second_used_sites, new_coord))
                        I.all_coordinates{1,end+1} = friends{1,friend};
                        second_used_sites = [second_used_sites, new_coord];
                    end
                end
            end
        end
        catalyzed_sites_mask = find(contains(second_used_sites, used_sites));

    end

    % Calculate probability distribution around each point
    prob_cloud = cell(1,size(I.all_coordinates,2));
    prop_prob_cloud = cell(1,size(I.all_coordinates,2));
    for nZero = 1:size(I.all_coordinates,2)
        missing_prob = 1;
        prop_missing_prob = 1;
        for nOne = 1:size(I.all_coordinates,2)

            t = R;  % Revise t to reflect diffusion rates (or make 2 versions of prob_cloud)
            propt = floor(R/sqrt(max(p.subcycles_per_ring)));   % THIS IS CHEATING, CORRECT THIS
            n01 = I.all_coordinates{nZero}{1}; n02 = I.all_coordinates{nZero}{2};
            n1 = I.all_coordinates{nOne}{1}; n2 = I.all_coordinates{nOne}{2};

            P = occupation_p(p.shape, R, 1, t, n01, n02, n1, n2);
            % ~~~ The following is a temporary approximation ~~~
            % ~~~ Figure out why probs are going negative ~~~~~~
            if P < 0, P = 0; end

            missing_prob = missing_prob - P;
            prob_cloud{nZero}{1,end+1} = P;

            propP = occupation_p(p.shape, R, 1, propt, n01, n02, n1, n2);
            % ~~~ The following is a temporary approximation ~~~
            % ~~~ Figure out why probs are going negative ~~~~~~
            if propP < 0, propP = 0; end

            prop_missing_prob = prop_missing_prob - propP;
            prop_prob_cloud{nZero}{1,end+1} = propP;
            
        end
        prob_cloud{nZero}{nZero} = prob_cloud{nZero}{nZero} + missing_prob; % Make total probs = 1. Higher probability of a particle staying in original site.
        prop_prob_cloud{nZero}{nZero} = prop_prob_cloud{nZero}{nZero} + prop_missing_prob;
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
            if coord == p.seed_locations(ring)
                indeces = matches(base_species, temp_member_species_record{ring});
                species_counts{1, coord}(2,indeces) = {p.seed_concentration(ring)};
            end
        end

    end

    % Stick outputs into I

    I.species_counts = species_counts;
    I.coordinate_list = coordinate_list;
    I.prob_cloud = prob_cloud;
    I.prop_prob_cloud = prop_prob_cloud;
    I.used_sites = used_sites;
    I.second_used_sites = second_used_sites;
    I.catalyzed_sites_mask = catalyzed_sites_mask;
    
end