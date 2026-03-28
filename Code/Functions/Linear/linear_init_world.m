function [p, I] = linear_init_world(p, I)

    % For linear sims with no parallelization
    % Initialize the dimensions of the world and all the data structures
    % that will hold concentration information

    % Also generate the probability distribution for the given coords

    arguments (Input)
        p;
        I;
    end

    if p.reboot == false

    % Initialize coordinate system
    I.coordinate_list = {};
    I.all_coordinates = {};

    if p.shape == "parallelogram"

        D = p.separation_distance;
        s = p.sites;
        R = s*D;
    
        for x = 0:D:(D*(s-1))
            for y = 0:D:(D*(s-1))
                I.coordinate_list{1,end+1} = {};
                I.coordinate_list{1,end}{1,1} = x;
                I.coordinate_list{1,end}{1,2} = y;
                I.coordinate_list{1,end}{1,3} = (0-x-y);
            end
        end

    elseif p.shape == "hex"

        D = p.separation_distance;
        s = p.sites;
        R = s*D + floor(D/2);   I.R = R;
        omega = 3*R*(R+1) + 1;  I.omega = omega;

        I.all_coordinates = cell(1,omega);
        I.empty_coord_mask = [];
        I.coordinate_list = {};
        I.catalyzed_sites_mask = [];
        I.origin = 0;

        % Set up pixels
        loops = 1;
        for x = -R:R
            for y = -R:R
                for z = -R:R
                    if (x+y+z) == 0
                        % I.all_coordinates{1,loops}{1} = x; I.all_coordinates{1,loops}{2} = y; I.all_coordinates{1,loops}{3} = z;
                        I.all_coordinates{1,loops} = [x,y,z];
                        if mod(abs(x),D) == 0 && mod(abs(y),D) == 0 && mod(abs(z),D) == 0
                            % I.coordinate_list{1,end+1}{1} = x; I.coordinate_list{1,end}{2} = y; I.coordinate_list{1,end}{3} = z;
                            I.coordinate_list{1,end+1} = [x,y,z];
                            I.catalyzed_sites_mask = [I.catalyzed_sites_mask, loops];
                        else
                            I.empty_coord_mask = [I.empty_coord_mask, loops];
                        end
                        if x == 0 && y == 0 && z == 0
                            I.origin = loops;
                        end
                        loops = loops + 1;
                    end
                end
            end
        end
        if length(I.catalyzed_sites_mask) ~= length(I.coordinate_list)
            fprintf("PROBLEM!!!! >:3")
        end
        clear x; clear y; clear z; clear loops;
       
    end

    % Calculate mirrored centers for propagation across periodic boundaries
    I.centres = cell(1,6);
    q = (2*R)+1; r = -R; s = -R-1;
    I.centres{1} = [q,r,s];
    I.centres{2} = [-r,-s,-q];
    I.centres{3} = [s,q,r];
    I.centres{4} = [-s,-q,-r];
    I.centres{5} = [r,s,q];
    I.centres{6} = [-q,-r,-s];

    I.R = R;
    clear R;

    % set up seed locations
    % Fix this so it is alterable in params later
    if p.seed_start == "origin"
        p.seed_locations = [I.origin, I.origin];
    else
        p.seed_locations = randsample(I.catalyzed_sites_mask,2,false);
    end

    % Add the missing dispersal reactions (and outflow)
    for s = 1:length(I.base_species)
        spec = I.base_species(s);
            % If a diffuse reaction is selected, the sim will handle where
            % the particle is sent separately
        if contains(spec, "_diff")
            I.reactions{end+1} = {p.flow_rate, {spec}, {1}, {spec}, {0}, "diffuse"};
            I.reactions{end+1} = {p.out_rate, {spec}, {1}, {spec}, {0}, "outflow"};
        elseif contains(spec, "prop_")
            % Apply the right diffusion rate
            splitprop = split(spec,"_"); cycle_num = str2double(splitprop(2));
            I.reactions{end+1} = {p.flow_rate/sqrt(p.subcycles_per_ring(cycle_num)), {spec}, {1}, {spec}, {0}, "diffuse"};
            if p.consistent_out_rate == true
                I.reactions{end+1} = {p.out_rate/sqrt(p.subcycles_per_ring(cycle_num)), {spec}, {1}, {spec}, {0}, "outflow"};
            else
                I.reactions{end+1} = {p.out_rate, {spec}, {1}, {spec}, {0}, "outflow"};
            end
        end
    end
    if p.chemostat == false
        I.reactions{end+1} = {p.in_rate, {"F"}, {1}, {"F"}, {0}, "outflow"};
    end
    clear splitprop; clear s; clear spec; clear cycle_num;

    % Make separate data structures to keep track of species in every
    % coordinate

    I.species_counts = single(zeros(size(I.all_coordinates,2), length(I.base_species)));

    % I.species_counts = cell(1, size(I.all_coordinates,2));
    % for coord = 1:size(I.all_coordinates,2)
    %     I.species_counts{1, coord} = cellstr(I.base_species);
    %     I.species_counts{1, coord}(2,:) = {0};
    % end
    % clear coord;


    for c = 1:length(I.catalyzed_sites_mask)
        % Set starting concentrations of sites and food
        coord = I.catalyzed_sites_mask(c);
        indeces = strcmp(I.base_species, 'F');
        I.species_counts(coord,indeces) = p.food_concentration;
        indeces = strcmp(I.base_species, 'site');
        I.species_counts(coord, indeces) = p.site_concentration;

        % Prepare seeds
        for ring = 1:length(I.ring_list)
            temp_member_species_record = I.member_species_record;
            for seed = 1:size(I.member_species_record{1,ring},2)
                temp_member_species_record{1,ring}{seed} = convertStringsToChars(temp_member_species_record{1,ring}{seed} + p.seed_state(ring));
            end
            if coord == p.seed_locations(ring)
            %if coord == I.origin
                indeces = matches(I.base_species, temp_member_species_record{ring});
                I.species_counts(coord,indeces) = p.seed_concentration(ring);
                if p.seed_state(ring) == "_ad"
                    site_index = matches(I.base_species, "site");
                    I.species_counts(coord,site_index) = (I.species_counts(coord,site_index) - p.seed_concentration(ring)*p.subcycles_per_ring(ring));
                end
            end
        end
        clear temp_member_species_record; clear c; clear ring; clear coord; clear indeces; clear site_index;

    end

    % Map species indeces to reactions for faster propensity updating
    I.species_index_map = {};

   for i = 1:length(I.reactions)
        reaction = i;
        I.species_index_map{1,end+1} = {};
        for reactant = 1:length(I.reactions{reaction}{2})
            index = strcmp(I.base_species(1,:), I.reactions{reaction}{2}{reactant});
            I.species_index_map{end}{end+1} = index;
        end
   end
   clear index; clear i; clear reactant;

   % Precalculate reaction map for faster simulating
   % Set up the update key
    I.reaction_update_key = cell(1,size(I.reactions,1));
    reaction_updates_calculated = zeros(1, size(I.reactions,1));

    for mu = 1:size(I.reactions, 1)
        affected_species = [];
        affected_reaction_indices = [];
        for reactant = 1:length(I.reactions{mu}{2})
            affected_species = [affected_species, I.reactions{mu}{2}{reactant}];
        end
        for product = 1:length(I.reactions{mu}{4})
            affected_species = [affected_species, I.reactions{mu}{4}{product}];
        end
        for affected_reaction = 1:size(I.reactions, 1)
            if any(ismember(string(I.reactions{affected_reaction}{2}), affected_species))
                affected_reaction_indices = [affected_reaction_indices, affected_reaction];
            end
        end
        I.reaction_update_key{mu} = affected_reaction_indices;
        reaction_updates_calculated(mu) = 1;
        clear affected_reaction_indices; clear affected_species;
        clear reactant; clear product; clear affected_reaction;
    end


    % Test: put EVERYTHING into arrays and handle reactions with a simple
    % stoichiometric array. Cells take too long to handle

    I.cell_reactions = I.reactions;
    
    I.rate_constants = single(zeros(1, I.num_reactions));
    I.reactions = single(zeros(I.num_reactions, length(I.base_species)));
    I.tags = [];

    for i = 1:size(I.cell_reactions,1)
        I.rate_constants(i) = single(I.cell_reactions{i}{1});
        I.tags = [I.tags, I.cell_reactions{i}{6}];
        % Make stoich arrays
        for reactant = 1:size(I.cell_reactions{i}{2}, 2)
            index = strcmp(I.base_species(1,:), I.cell_reactions{i}{2}{reactant});
            I.reactions(i,index) = int8(-I.cell_reactions{i}{3}{reactant});
        end
        for product = 1:size(I.cell_reactions{i}{4}, 2)
            index = strcmp(I.base_species(1,:), I.cell_reactions{i}{4}{product});
            I.reactions(i,index) = int8(I.cell_reactions{i}{5}{product});
        end
    end


    else
        fprintf("Rebooting, skipped world generation <w<\n");
    end
    
end