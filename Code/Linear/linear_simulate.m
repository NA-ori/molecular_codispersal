function [O, I] = linear_simulate(p, I, O)

    % Simulation with no parallelization of reactions, everything runs in
    % sequence with dispersal handled as reactions

    arguments (Input)
        p;
        I;
        O;
    end


    %~~~~~ Moving on to the Simulation ~~~~~%

    if p.reboot == false

    % Initial variables
    O.sim_id = randi([1,99999999]);
    save_name = "savefile_" + O.sim_id + ".mat";
    
    if p.disturb_freq ~= 0
        next_disturbance = exprnd(p.disturb_freq);
    else
        next_disturbance = Inf;
    end

    error = false;
    t = 0;
    time = [t];
    next_recording = p.recording_freq;
    sitedex = strcmp(I.base_species, "site");
    foodex = strcmp(I.base_species, "F");

    num_reactions = size(I.reactions,1);
    num_coords = size(I.all_coordinates,2);

    % Set up propensity and concentration trackers
    current_chemical_counts = I.species_counts;

    current_reaction_propensities = zeros(1,num_coords*num_reactions);
    all_reactions = 1:num_reactions;
    for coord = 1:num_coords
        propensity_indices = all_reactions;
        for index = 1:length(propensity_indices)
            propensity_indices(index) = (coord-1)*num_reactions + propensity_indices(index);
        end
        current_reaction_propensities = linear_update_propensities(coord,all_reactions,current_chemical_counts,current_reaction_propensities,propensity_indices,I);
    end
    
    concentration_tracker = cell(1,num_coords);
    for cel = 1:size(concentration_tracker,2)
        concentration_tracker(1,cel) = {current_chemical_counts(cel,:)};
    end

    else

        fprintf("Attempting to restart from " + p.reboot_from + "^w^\n");
        load(p.reboot_from);
        fprintf("Restart successful! ~W~\n");
    end

    O.run_time = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % THE ACTUAL SIMULATION %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    propensity_length = 1:length(current_reaction_propensities);

    while t <= p.t_max

        tic;
        tau = (1/(sum(current_reaction_propensities))) * log(1/rand());

        % Try an alternative method of choosing a reaction

        try
            mu = randsample(propensity_length, 1, true, current_reaction_propensities);
        catch
            if sum(current_reaction_propensities) == 0
                fprintf("Ended simulations because everything went extinct >~<\n");
                O.end_time = t;
                O.incomplete_sim = false;
                break
            else
                error = true;
            end
        end

        % Now figure out which coord things are happening in to update the
        % appropriate values in the concentration tracker
        absolute_coord = ceil(mu/size(I.reactions,1));
        reac = mod(mu,size(I.reactions,1));
        if reac == 0
            reac = num_reactions;
        end

        if I.tags(reac) == "diffuse"

            % Remove the diffused particle from its home :(
            current_chemical_counts(absolute_coord,:) = current_chemical_counts(absolute_coord,:) + I.reactions(reac,:);

            if any(current_chemical_counts(absolute_coord,:)<0)
                fprintf("Negative concentration detected in coord " + num2str(absolute_coord) + " >:3\n");
                %fprintf("Species " + I.reactions{reac}{2}{reactant} + " reduced to " + current_chemical_counts{absolute_coord}{2,index} + "\n");
                fprintf("Offending reaction: " + num2str(reac) + "\n");
                O.current_reaction_propensities = current_reaction_propensities;
                O.current_chemical_counts = current_chemical_counts;
                O.time = time;
                error = true;
            end

            reactions_to_update = I.reaction_update_key{reac};
            propensity_indices = reactions_to_update;
            for index = 1:length(propensity_indices)
                propensity_indices(index) = (absolute_coord-1)*size(I.reactions,1) + propensity_indices(index);
            end
            clear index;
            current_reaction_propensities = linear_update_propensities(absolute_coord,reactions_to_update,current_chemical_counts,current_reaction_propensities,propensity_indices,I);

            % Randomly determine which neighboring pixel it will be sent to
            
            friends = neighbors(I.all_coordinates{absolute_coord}(1), I.all_coordinates{absolute_coord}(2), I.all_coordinates{absolute_coord}(3), 1, I);
            sampled_coord = randsample(friends,1); sampled_coord = sampled_coord{1};
            % compfunction = @(x) isequal(x, sampled_coord); index = cellfun(compfunction, I.all_coordinates); clear sampled_coord;
            for coord = 1:num_coords
                if isequal(I.all_coordinates{coord}, sampled_coord)
                    index = coord;
                end
            end
            new_coord = index; clear index; clear sampled_coord; clear friends;

            % Add the particle to the new pixel
            current_chemical_counts(new_coord,:) = current_chemical_counts(new_coord,:) - (I.reactions(reac,:));


            % Update propensities in the new pixel

            reactions_to_update = I.reaction_update_key{reac};
            propensity_indices = reactions_to_update;
            for index = 1:length(propensity_indices)
                propensity_indices(index) = (new_coord-1)*size(I.reactions,1) + propensity_indices(index);
            end
            current_reaction_propensities = linear_update_propensities(new_coord,reactions_to_update,current_chemical_counts,current_reaction_propensities,propensity_indices,I);

            clear reactions_to_update; clear propensity_indices; clear index;

        else

            % Change concentrations and update propensities
            current_chemical_counts(absolute_coord,:) = current_chemical_counts(absolute_coord,:) + I.reactions(reac,:);

            if any(current_chemical_counts(absolute_coord,:)<0)
                fprintf("Negative concentration detected in coord " + num2str(absolute_coord) + " >:3\n");
                %fprintf("Species " + I.reactions{reac}{2}{reactant} + " reduced to " + current_chemical_counts{absolute_coord}{2,index} + "\n");
                fprintf("Offending reaction: " + num2str(reac) + "\n");
                O.current_reaction_propensities = current_reaction_propensities;
                O.current_chemical_counts = current_chemical_counts;
                O.time = time;
                error = true;
            end
                
            % Reset food if chemostatted
            if any(ismember(I.catalyzed_sites_mask, absolute_coord))
                current_chemical_counts(absolute_coord,foodex) = p.food_concentration;
            end

            if error == true
                break
            end

            reactions_to_update = I.reaction_update_key{reac};
            propensity_indices = reactions_to_update;
            for index = 1:length(propensity_indices)
                propensity_indices(index) = (absolute_coord-1)*size(I.reactions,1) + propensity_indices(index);
            end
            current_reaction_propensities = linear_update_propensities(absolute_coord,reactions_to_update,current_chemical_counts,current_reaction_propensities,propensity_indices,I);
            clear propensity_indices; clear reactions_to_update; clear index;

        end

        if error == true
            fprintf("Simulation ended because an error was encountered! ;__;\n")
            break
        end

        % ~~~~~~~~~~~~ Disturb ~~~~~~~~~~~~ %

        if t > next_disturbance

            disturbed_coord = randsample(I.catalyzed_sites_mask, 1);

            current_chemical_counts = disturb(p, I, current_chemical_counts, disturbed_coord);

            % Update propensities
            reactions_to_update = 1:num_reactions;
            propensity_indices = reactions_to_update;
            for index = 1:length(propensity_indices)
                propensity_indices(index) = (disturbed_coord-1)*size(I.reactions,1) + propensity_indices(index);
            end
            current_reaction_propensities = linear_update_propensities(disturbed_coord,reactions_to_update,current_chemical_counts,current_reaction_propensities,propensity_indices,I);
            if p.introspection == true, fprintf("Disturbance in " + disturbed_coord + " >w<!\n"); end
            clear propensity_indices; clear reactions_to_update; clear index; clear disturbed_coord;

            next_disturbance = next_disturbance + exprnd(p.disturb_freq);

        end

        
        % ~~~~~~~ End of disturbance ~~~~~~ %

        t = t + tau;

        if t >= next_recording
            for c = 1:size(I.all_coordinates,2)
                concentration_tracker{1,c}(end+1,:) = current_chemical_counts(c,:);
            end
            clear c;
            next_recording = next_recording + p.recording_freq;
            time = [time; t];
            if p.introspection == true, fprintf(t + "/" + p.t_max + "\n"); end
        end

        clear mu; clear tau;

        O.run_time = O.run_time + toc;

        if t > p.t_max
            O.end_time = t;
            O.incomplete_sim = false;
            break
        elseif O.run_time / 60 / 60 >= p.stop_time_hrs
            fprintf("Ended simulation because it took " + p.stop_time_hrs + " hours! >u<\n");
            O.end_time = t;
            save(save_name, "p", "I", "O", "next_disturbance", "error", "t", "time", "next_recording", ...
                "sitedex", "num_reactions", "num_coords", "current_reaction_propensities", "current_chemical_counts", "concentration_tracker", "save_name");
            O.incomplete_sim = true;
            break
        end

        % Stop if everything is done, but only under certain parameters
        % Sims with no disturbance don't change once everything is filled,
        % so end them once everything is filled

        if p.disturb_freq == 0 && p.allow_stopping == true
            sitesum = 0;
            for coord = 1:size(I.all_coordinates,2)
                sitesum = sitesum + current_chemical_counts(coord,sitedex);
            end
            unoccupied = sitesum / (p.site_concentration*size(I.coordinate_list,2));
            if unoccupied == 0
                fprintf("Simulation ended because all sites were filled! >W<\n");
                for c = 1:size(I.all_coordinates,2)
                    concentration_tracker{1,c}(end+1,:) = current_chemical_counts(c,:);
                end
                clear c;
                time = [time; t];
                O.end_time = t;
                O.incomplete_sim = false;
                break
            end
            clear unoccupied; clear sitesum; clear coord;
        end

    end

    O.time = time; clear time;
    O.concentration_tracker = concentration_tracker; clear concentration_tracker;

    clear current_reaction_propensities;
    clear current_chemical_counts;

end