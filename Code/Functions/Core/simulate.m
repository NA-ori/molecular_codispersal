function [O, I] = simulate(p, I, O, prob_cloud, prop_prob_cloud)

    arguments (Input)
        p;
        I;
        O;
        prob_cloud;
        prop_prob_cloud;
    end


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

    % Set up the update key
    reaction_update_key = cell(1,size(I.reactions,1));
    reaction_updates_calculated = zeros(1, size(I.reactions,1));

    % Initialize propensity tracker (and initial reaction propensities)
    current_chemical_counts = I.species_counts;
    current_reaction_propensities = cell(1, length(I.coordinate_list));
    
    for coord = 1:size(I.coordinate_list, 2)

        % Set up data structure
        current_reaction_propensities{1,coord} = cell(1,size(I.reactions,1));
        current_reaction_propensities{1,coord}(1,:) = {0};

        % Set the starting reaction propensities
        all_reactions = 1:size(current_reaction_propensities{1,coord},2);
        current_reaction_propensities = update_propensities(all_reactions, coord, I.reactions, current_chemical_counts, current_reaction_propensities,I);

    end

    concentration_tracker = current_chemical_counts;

    % Some stuff for handling diffusion
    update_mask = cell(1,length(I.all_coordinates));
    for coord = 1:size(update_mask,2)
        update_mask{1, coord} = cellstr(I.base_species);
        update_mask{1, coord}(2,:) = {0};
    end
    diff_mask = contains(I.base_species, "_diff");
    diff_mask(find(contains(I.base_species, "prop_"),1)) = 1;
    I.diff_mask = diff_mask;
    I.out_mask = diff_mask; I.out_mask(find(contains(I.base_species, "F"),1)) = 1;

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

        for coord = 1:length(I.coordinate_list)
            absolute_coord = I.catalyzed_sites_mask(coord);
            local_time = 0;

            while local_time <= p.disperse_frequency    % Lesser loop
    
                r_1 = rand();
                r_2 = rand();

                a_0 = sum(cell2mat(current_reaction_propensities{coord}));
    
                tau = (1/a_0) * log(1/r_1);
    
                % now determine which reaction will occur
                mu = 0;
                summ = 0;
                target = r_2*a_0;   % target is some random percentage of the total propensities a_0
                for reaction = 1:size(I.reactions, 1)
                    summ = summ + current_reaction_propensities{coord}{reaction};
                    if summ > target
                        mu = reaction;  % the reaction that happens to pass the target value first will occur.
                        break
                    end
                end

                % Change concentrations and update propensities
                for reactant = 1:size(I.reactions{mu}{2}, 2)
                    % subtract the number of particles that react from the current counts
                    index = strcmp(current_chemical_counts{absolute_coord}(1,:), I.reactions{mu}{2}{reactant});
                    current_chemical_counts{absolute_coord}{2,index} = current_chemical_counts{absolute_coord}{2,index} - I.reactions{mu}{3}{reactant};
                    if current_chemical_counts{absolute_coord}{2,index} < 0
                        fprintf("Negative concentration detected >:3\n");
                        fprintf("Offending reaction: " + num2str(mu) + "\n");
                        fprintf("Offending propensity: " + current_reaction_propensities{absolute_coord}(mu) + "\n");
                        error = true;
                    end
                end
                for product = 1:size(I.reactions{mu}{4}, 2)
                    % Add number of products formed to current counts
                    index = strcmp(current_chemical_counts{absolute_coord}(1,:), I.reactions{mu}{4}{product});
                    current_chemical_counts{absolute_coord}{2,index} = current_chemical_counts{absolute_coord}{2,index} + I.reactions{mu}{5}{product};
                end
                if error == true
                    break
                end
                

                % Update the propensities
                if reaction_updates_calculated(mu) == 0
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
                    reaction_update_key{mu} = affected_reaction_indices;
                    reaction_updates_calculated(mu) = 1;
                    reactions_to_update = reaction_update_key{mu};
                elseif reaction_updates_calculated(mu) == 1
                    reactions_to_update = reaction_update_key{mu};
                else
                    error = true;
                end
                if error == true
                    break;
                end
                current_reaction_propensities = update_propensities(reactions_to_update,coord,I.reactions,current_chemical_counts,current_reaction_propensities,I);

                local_time = local_time + tau;
                if local_time > p.disperse_frequency
                    local_time = p.disperse_frequency;
                    % concentration_tracker{absolute_coord}(end+1,:) = current_chemical_counts{absolute_coord}(2,:);
                    break
                end

            end
        end

        if error == true
            fprintf("Simulation ended because an error was encountered! ;__;\n")
            break
        end
        
        t = t + p.disperse_frequency;
        % if p.introspection == true, fprintf(t + "/" + p.t_max + "\n"); end
        % time = [time; t];
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

        current_chemical_counts = outflow(p, I.out_mask, current_chemical_counts, I);

        % Update propensities
        for coord = 1:size(I.coordinate_list, 2)
            all_reactions = 1:size(current_reaction_propensities{1,coord},2);
            current_reaction_propensities = update_propensities(all_reactions, coord, I.reactions, current_chemical_counts, current_reaction_propensities,I);
        end


        % ~~~~~~~~~~~~~~~~~ Disperse! ~~~~~~~~~~~~~~~~~
        
        current_chemical_counts = disperse(p, update_mask, diff_mask, current_chemical_counts, prob_cloud, prop_prob_cloud);

        % Update propensities
        for coord = 1:size(I.coordinate_list, 2)
            all_reactions = 1:size(current_reaction_propensities{1,coord},2);
            current_reaction_propensities = update_propensities(all_reactions, coord, I.reactions, current_chemical_counts, current_reaction_propensities,I);
        end


        % ~~~~~~~~~~~~~~~~~ Disturb! ~~~~~~~~~~~~~~~~~~

        if t > next_disturbance

            disturbed_coord = randsample(I.catalyzed_sites_mask, 1);
            current_reaction_propensities = disturb(p, current_chemical_counts, disturbed_coord, I.base_species);

            % Update propensities
            all_reactions = 1:size(current_reaction_propensities{1,disturbed_coord},2);
            current_reaction_propensities = update_propensities(all_reactions, disturbed_coord, I.reactions, current_chemical_counts, current_reaction_propensities,I);

        end


        % ~~~~~~~~~~~~~~~ Time to stop? ~~~~~~~~~~~~~~~

        if t >= p.t_max
            break
        end

    end

    % Add variables to output structures

    I.current_reaction_propensities = current_reaction_propensities;
    I.current_chemical_counts = current_chemical_counts;
    clear current_reaction_propensities;
    clear current_chemical_counts;

    O.time = time;
    O.concentration_tracker = concentration_tracker;

end