function [O, I] = par_simulate(p, I, c, O)

    % Simulation with parallelization of reactions
    % Autocatalytic reactions can run simultaneously; then sim will enter a
    % dispersal/propagule breaking phase

    arguments (Input)
        p;
        I;
        c;
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

        % Some stuff for diffusion and easier indexing

        diff_mask = contains(I.base_species, "_diff"); diff_ind = find(diff_mask);
        prop_mask = contains(I.base_species, "prop_"); prop_ind = find(prop_mask);
        comb_mask = diff_mask+prop_mask; comb_ind = find(comb_mask);

        prop_break_reaction = find(strcmp(I.tags, "prop_break"));
    
        num_reactions = size(I.reactions,1);
        num_coords = size(I.all_coordinates,2);

        % leave_prob = p.out_rate / p.flow_rate + p.out_rate;
        leave_prob = 0.5;
        stay_prob = 1 - leave_prob;
    
        % Set up propensity and concentration trackers
        current_chemical_counts = I.species_counts;
    
        % Set up a 2-dimensional reaction propensity array
    
        current_reaction_propensities = zeros(num_coords,num_reactions);
        all_reactions = 1:num_reactions;
    
        for coord = 1:num_coords
            current_reaction_propensities = update_propensities(coord,all_reactions,current_chemical_counts,current_reaction_propensities,I);
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

    % propensity_length = 1:length(current_reaction_propensities);

    while t <= p.t_max      % Greater loop

        tic;
        for co = 1:length(I.catalyzed_sites_mask)
            coord = I.catalyzed_sites_mask(co);
            local_time = 0;

            while local_time <= p.chunk_time    % Lesser loop
        
                tau = (1/(sum(current_reaction_propensities(coord,:)))) * log(1/rand());
        
                % Try an alternative method of choosing a reaction
        
                try
                    reac = randsample(num_reactions, 1, true, current_reaction_propensities(coord,:));
                catch
                    if sum(current_reaction_propensities) == 0
                        fprintf("Stopped running chunk because there's nothing in it! >~<\n");
                        break
                    else
                        error = true;
                    end
                end
                
                % Change concentrations and update propensities
                current_chemical_counts(coord,:) = current_chemical_counts(coord,:) + I.reactions(reac,:);
    
                if any(current_chemical_counts(coord,:)<0)
                    fprintf("Negative concentration detected in coord " + num2str(absolute_coord) + " >:3\n");
                    fprintf("Offending reaction: " + num2str(reac) + "\n");
                    O.current_reaction_propensities = current_reaction_propensities;
                    O.current_chemical_counts = current_chemical_counts;
                    O.time = time;
                    error = true;
                end
                    
                % Reset food if chemostatted
                if p.chemostat == true
                    current_chemical_counts(coord,foodex) = p.food_concentration;
                end
    
                if error == true
                    break
                end
    
                reactions_to_update = I.reaction_update_key{reac};
                current_reaction_propensities = update_propensities(coord,reactions_to_update,current_chemical_counts,current_reaction_propensities,I);
                clear reactions_to_update;

                local_time = local_time + tau;
                if local_time > p.chunk_time
                    break
                end

            end
        end

        if error == true
            fprintf("Simulation ended because an error was encountered! ;__;\n")
            break

        end

        % ~~~~~ Disperse! ~~~~~ %
        
        current_chemical_counts = disperse(I, num_coords, diff_ind, current_chemical_counts, c.cloud);   % Regular particles
        current_chemical_counts = disperse(I, num_coords, prop_ind, current_chemical_counts, c.prop_cloud);   % Propagules

        % ~~~~~~ Outflow! ~~~~~~ %
        
        current_chemical_counts = outflow(p, comb_ind, leave_prob, stay_prob, num_coords, current_chemical_counts, I);

        % ~~~~~~ Break Propagules! ~~~~~~ %

        current_chemical_counts = break_props(p, I, prop_ind, prop_break_reaction, current_chemical_counts);

        % ~~~~~~~~~~~~ Disturb ~~~~~~~~~~~~ %

        if t > next_disturbance

            disturbed_coord = randsample(I.catalyzed_sites_mask, 1);

            current_chemical_counts = disturb(p, I, current_chemical_counts, disturbed_coord);

            if p.introspection == true, fprintf("Disturbance in " + disturbed_coord + " >w<!\n"); end

            clear disturbed_coord;

            next_disturbance = next_disturbance + exprnd(p.disturb_freq);

        end

        % ~~~~~~~ End of disturbance ~~~~~~ %


        % Update propensities after all that nonsense
        for coord = 1:num_coords
            current_reaction_propensities = update_propensities(coord,1:num_reactions,current_chemical_counts,current_reaction_propensities,I);
        end


        % ~~~~ Handle end-loop stuff ~~~~ %

        t = t + p.chunk_time;

        if t >= next_recording
            for co = 1:num_coords
                concentration_tracker{1,co}(end+1,:) = current_chemical_counts(co,:);
            end
            clear co;
            next_recording = next_recording + p.recording_freq;
            time = [time; t];
            if p.introspection == true, fprintf(t + "/" + p.t_max + "\n"); end
        end

        clear reac; clear tau;

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
        % 
        if p.disturb_freq == 0 && p.allow_stopping == true
            sitesum = 0;
            for coord = I.catalyzed_sites_mask
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

    O.time = time;
    O.concentration_tracker = concentration_tracker;

end