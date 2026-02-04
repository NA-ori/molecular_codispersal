function [O, I] = linear_simulate(p, I, O)

    % Simulation with no parallelization of reactions, everything runs in
    % sequence with dispersal handled as reactions

    arguments (Input)
        p;
        I;
        O;
    end

    %~~~~~ Moving on to the Simulation ~~~~~%

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
    
    concentration_tracker = current_chemical_counts;

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % THE ACTUAL SIMULATION %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    while t <= p.t_max
        r_1 = rand();
        r_2 = rand();

        a_0 = sum(current_reaction_propensities);
        tau = (1/a_0) * log(1/r_1);

        % Try an alternative method of choosing a reaction
        % mu = randsample(1:length(current_reaction_propensities), 1, true, current_reaction_propensities);

        % now determine which reaction will occur
        mu = 0;
        summ = 0;
        target = r_2*a_0;   % target is some random percentage of the total propensities a_0
        for reaction = 1:length(current_reaction_propensities)
            summ = summ + current_reaction_propensities(reaction);
            if summ > target
                mu = reaction;  % the reaction that happens to pass the target value first will occur.
                break
            end
        end

        % Now figure out which coord things are happening in to update the
        % appropriate values in the concentration tracker
        absolute_coord = ceil(mu/size(I.reactions,1));
        reac = mod(mu,size(I.reactions,1));
        if reac == 0
            reac = num_reactions;
        end

        if I.reactions{reac}{6} == "diffuse"

            % Remove the diffused particle from its home :(
            for reactant = 1:size(I.reactions{reac}{2}, 2)
                % subtract the number of particles that react from the current counts
                index = strcmp(current_chemical_counts{absolute_coord}(1,:), I.reactions{reac}{2}{reactant});
                current_chemical_counts{absolute_coord}{2,index} = current_chemical_counts{absolute_coord}{2,index} - 1;
                if current_chemical_counts{absolute_coord}{2,index} < 0
                    fprintf("Negative concentration detected >:3\n");
                    fprintf("Offending reaction: " + num2str(reac) + "\n");
                    O.current_reaction_propensities = current_reaction_propensities;
                    O.current_chemical_counts = current_chemical_counts;
                    error = true;
                end
            end

            reactions_to_update = I.reaction_update_key{reac};
            propensity_indices = reactions_to_update;
            for index = 1:length(propensity_indices)
                propensity_indices(index) = (absolute_coord-1)*size(I.reactions,1) + propensity_indices(index);
            end
            current_reaction_propensities = linear_update_propensities(absolute_coord,reactions_to_update,current_chemical_counts,current_reaction_propensities,propensity_indices,I);

            % Randomly determine which neighboring pixel it will be sent to
            
            friends = neighbors(I.all_coordinates{absolute_coord}(1), I.all_coordinates{absolute_coord}(2), I.all_coordinates{absolute_coord}(3), 1, I);
            sampled_coord = randsample(friends,1); sampled_coord = sampled_coord{1};
            compfunction = @(x) isequal(x, sampled_coord); index = cellfun(compfunction, I.all_coordinates); clear sampled_coord;
            new_coord = find(index);

            % Add the particle to the new pixel

            for product = 1:size(I.reactions{reac}{4}, 2)
                % Add number of products formed to current counts
                index = strcmp(current_chemical_counts{absolute_coord}(1,:), I.reactions{reac}{4}{product});
                current_chemical_counts{new_coord}{2,index} = current_chemical_counts{new_coord}{2,index} + 1;
            end

            % Update propensities in the new pixel

            reactions_to_update = I.reaction_update_key{reac};
            propensity_indices = reactions_to_update;
            for index = 1:length(propensity_indices)
                propensity_indices(index) = (new_coord-1)*size(I.reactions,1) + propensity_indices(index);
            end
            current_reaction_propensities = linear_update_propensities(new_coord,reactions_to_update,current_chemical_counts,current_reaction_propensities,propensity_indices,I);


        else

            % Change concentrations and update propensities
            for reactant = 1:size(I.reactions{reac}{2}, 2)
                % subtract the number of particles that react from the current counts
                index = strcmp(current_chemical_counts{absolute_coord}(1,:), I.reactions{reac}{2}{reactant});
                current_chemical_counts{absolute_coord}{2,index} = current_chemical_counts{absolute_coord}{2,index} - I.reactions{reac}{3}{reactant};
                if current_chemical_counts{absolute_coord}{2,index} < 0
                    fprintf("Negative concentration detected in coord " + num2str(absolute_coord) + " >:3\n");
                    fprintf("Species " + I.reactions{reac}{2}{reactant} + " reduced to " + current_chemical_counts{absolute_coord}{2,index} + "\n");
                    fprintf("Offending reaction: " + num2str(reac) + "\n");
                    O.current_reaction_propensities = current_reaction_propensities;
                    O.current_chemical_counts = current_chemical_counts;
                    O.time = time;
                    error = true;
                end
            end
            for product = 1:size(I.reactions{reac}{4}, 2)
                % Add number of products formed to current counts
                index = strcmp(current_chemical_counts{absolute_coord}(1,:), I.reactions{reac}{4}{product});
                current_chemical_counts{absolute_coord}{2,index} = current_chemical_counts{absolute_coord}{2,index} + I.reactions{reac}{5}{product};
            end
            % Reset food concentration if chemostatted
            if p.chemostat == true
                index = strcmp(I.base_species(1,:), "F");
                current_chemical_counts{absolute_coord}{2,index} = p.food_concentration;
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

        end


        if error == true
            fprintf("Simulation ended because an error was encountered! ;__;\n")
            break
        end

        % ~~~~~~~~~~~~ Disturb ~~~~~~~~~~~~ %

        
        % ~~~~~~~ End of disturbance ~~~~~~ %

        t = t + tau;

        if t >= next_recording
            for c = 1:size(I.all_coordinates,2)
                concentration_tracker{c}(end+1,:) = current_chemical_counts{c}(2,:);
            end
            next_recording = next_recording + p.recording_freq;
            time = [time; t];
            if p.introspection == true, fprintf(t + "/" + p.t_max + "\n"); end
        end

        if t > p.t_max
            break
        end

    end

    O.time = time; clear time;
    O.concentration_tracker = concentration_tracker; clear concentration_tracker;

    clear current_reaction_propensities;
    clear current_chemical_counts;

end