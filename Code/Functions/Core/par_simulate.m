function [O, I] = par_simulate(p, I, O, prob_cloud, prop_prob_cloud)

    arguments (Input)
        p;
        I;
        O;
        prob_cloud;
        prop_prob_cloud;
    end

    % Extract variables from I
    base_species = I.base_species;
    species_counts = I.species_counts;
    reactions = I.reactions;
    coordinate_list = I.coordinate_list;


    % Initial variables
    sample_interval = p.t_max / p.sample_number-1;
    next_sample_time = sample_interval;

    if p.disturb_freq ~= 0
        next_disturbance = exprnd(p.disturb_freq);
    else
        next_disturbance = Inf;
    end

    error = false;
    t = 0;
    time = [t];
    next_recording = p.recording_freq;

    % % Set up the update key
    % reaction_update_key = cell(1,size(reactions,1));
    % reaction_updates_calculated = zeros(1, size(reactions,1));
    reaction_update_key = I.reaction_update_key;

    % Initialize propensity tracker (and initial reaction propensities)
    current_chemical_counts = species_counts;
    current_reaction_propensities = cell(1, length(coordinate_list));
    
    for coord = 1:size(coordinate_list, 2)

        % Set up data structure
        current_reaction_propensities{1,coord} = cell(1,size(reactions,1));
        current_reaction_propensities{1,coord}(1,:) = {0};

        % Set the starting reaction propensities
        all_reactions = 1:size(current_reaction_propensities{1,coord},2);
        current_reaction_propensities = update_propensities(all_reactions, coord, reactions, current_chemical_counts, current_reaction_propensities,I);

    end

    concentration_tracker = current_chemical_counts;

    % Some stuff for handling diffusion
    update_mask = cell(1,length(I.all_coordinates));
    for coord = 1:size(update_mask,2)
        update_mask{1, coord} = cellstr(base_species);
        update_mask{1, coord}(2,:) = {0};
    end
    diff_mask = contains(base_species, "_diff");
    diff_mask(find(contains(base_species, "prop_"),1)) = 1;
    I.diff_mask = diff_mask;
    I.out_mask = diff_mask; I.out_mask(find(contains(base_species, "F"),1)) = 1;
    out_mask = I.out_mask;

    for reaction = 1:size(I.reactions,1)
        if I.reactions{reaction}{6} == "prop_break"
            I.prop_break_reaction = reaction;
            break
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%
    % THE ACTUAL SIMULATION %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    while t <= p.t_max      % Greater loop

        parfor coord = 1:length(coordinate_list)
            [I, current_reaction_propensities, current_chemical_counts] = par_sim_innards(p, I, current_reaction_propensities, current_chemical_counts, coord)
        end

        % if error == true
        %     fprintf("Simulation ended because an error was encountered! ;__;\n")
        %     break
        % end
        
        t = t + p.disperse_frequency;
        if t >= next_recording
            for c = 1:size(I.all_coordinates,2)
                concentration_tracker{c}(end+1,:) = current_chemical_counts{c}(2,:);
            end
            next_recording = next_recording + p.recording_freq;
            time = [time; t];
            if p.introspection == true, fprintf(t + "/" + p.t_max + "\n"); end
        end

        % ~~~~~~~~~~~~~ Break Propagules! ~~~~~~~~~~~~~

        current_chemical_counts = break_props(p, I, current_chemical_counts);

        % ~~~~~~~~~~~~~~~~~ Outflow! ~~~~~~~~~~~~~~~~~~

        current_chemical_counts = outflow(p, out_mask, current_chemical_counts, I);

        % Update propensities
        for coordin = 1:size(coordinate_list, 2)
            all_reactions = 1:size(current_reaction_propensities{1,coordin},2);
            current_reaction_propensities = update_propensities(all_reactions, coordin, reactions, current_chemical_counts, current_reaction_propensities,I);
        end

        % ~~~~~~~~~~~~~~~~~ Disperse! ~~~~~~~~~~~~~~~~~

        current_chemical_counts = disperse(p, update_mask, diff_mask, current_chemical_counts, prob_cloud, prop_prob_cloud);

        % Update propensities
        for coordin = 1:size(coordinate_list, 2)
            all_reactions = 1:size(current_reaction_propensities{1,coordin},2);
            current_reaction_propensities = update_propensities(all_reactions, coordin, reactions, current_chemical_counts, current_reaction_propensities,I);
        end

        % ~~~~~~~~~~~~~~~~~ Disturb! ~~~~~~~~~~~~~~~~~~

        if t > next_disturbance

            disturbed_coord = randsample(I.catalyzed_sites_mask, 1);
            current_reaction_propensities = disturb(p, current_chemical_counts, disturbed_coord, base_species);

            % Update propensities
            all_reactions = 1:size(current_reaction_propensities{1,disturbed_coord},2);
            current_reaction_propensities = update_propensities(all_reactions, disturbed_coord, reactions, current_chemical_counts, current_reaction_propensities,I);

        end


        % ~~~~~~~~~~~~~~~ Time to stop? ~~~~~~~~~~~~~~~

        if t >= p.t_max
            break
        end

    end

    % Add variables to output structures

    I.species_counts = species_counts;
    I.current_reaction_propensities = current_reaction_propensities;
    I.current_chemical_counts = current_chemical_counts;

    O.time = time;
    O.concentration_tracker = concentration_tracker;

end