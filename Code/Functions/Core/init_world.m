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

        coordinate_list{1,end+1} = {};
        coordinate_list{1,end}{1,1} = 0;
        coordinate_list{1,end}{1,2} = 0;
        coordinate_list{1,end}{1,3} = 0;

        used_sites = ["000"];

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
    end

    % Calculate probability distribution around each point
    prob_cloud = cell(1,size(coordinate_list,2));
    prop_prob_cloud = cell(1,size(coordinate_list,2));
    for nZero = 1:size(coordinate_list,2)
        missing_prob = 1;
        prop_missing_prob = 1;
        for nOne = 1:size(coordinate_list,2)

            t = R;  % Revise t to reflect diffusion rates (or make 2 versions of prob_cloud)
            propt = floor(R/sqrt(max(p.subcycles_per_ring)));   % THIS IS CHEATING, CORRECT THIS
            n01 = coordinate_list{nZero}{1}; n02 = coordinate_list{nZero}{2};
            n1 = coordinate_list{nOne}{1}; n2 = coordinate_list{nOne}{2};

            P = occupation_p(p.shape, R, 1, t, n01, n02, n1, n2);
            missing_prob = missing_prob - P;
            prob_cloud{nZero}{1,end+1} = P;

            propP = occupation_p(p.shape, R, 1, propt, n01, n02, n1, n2);
            prop_missing_prob = prop_missing_prob - propP;
            prop_prob_cloud{nZero}{1,end+1} = propP;
            
        end
        prob_cloud{nZero}{nZero} = prob_cloud{nZero}{nZero} + missing_prob; % Make total probs = 1. Higher probability of a particle staying in original site.
        prop_prob_cloud{nZero}{nZero} = prop_prob_cloud{nZero}{nZero} + prop_missing_prob;
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

    % Stick outputs into I

    I.species_counts = species_counts;
    I.coordinate_list = coordinate_list;
    I.prob_cloud = prob_cloud;
    I.prop_prob_cloud = prop_prob_cloud;
    
end