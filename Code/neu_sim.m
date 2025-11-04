function [reactions, ring_list, base_species, species_counts, member_species_record, current_reaction_propensities, current_chemical_counts, concentration_tracker, p] = neu_sim(options)

% Smarter version of the original CHTC sim.
% Reduces hard-coding and therefore hopefully mistakes.
% Also, put the whole thing in a function so every value is easily
% changeable.

% Ideally near-everything should be alterable just by changing the function
% inputs.

    arguments (Input)
        % Inputs are modified into lists to allow more flexibility
        options.dim = 10;
        options.t_max = 10000;
        options.default_concentration = 0;
        options.food_concentration = 500;
        options.seed_concentration = [25,25];
        options.site_concentration = 500;
        options.independence_disadvantage = [Inf];
        options.flow_rate = 1;
        options.in_rate = 0.1;
        options.out_rate = 1;
        options.adsorb_rate = 0.01;
        options.fierce_reaction_rate = 0.01;
        options.prop_formation_rate = 0.01/5;
        options.prop_funnel_rate = 0.1;
        options.seed_locations = [4,6];
        options.seed_state = ["_diff", "_diff"];
        options.disturb_freq = Inf;
        options.sample_number = 1000;
        %options.disperse_frequency = options.t_max / options.sample_number; % How frequently dispersal will happen

    end
    
    rng("shuffle");     % Ensure stochastics sims are different every time
    
    % Simply convert the options into the famed parameter struct p

    p = options;
    p.disperse_frequency = p.t_max / p.sample_number;

    fprintf("Initializing reactions >w<\n");
    [reactions, ring_list, member_species_record, base_species] = init_reactions(p);
    fprintf("Initializing world >w<\n");
    [species_counts, coordinate_list] = init_world(p, reactions, ring_list, base_species, member_species_record);

    [species_counts, time, current_reaction_propensities, current_chemical_counts, concentration_tracker] = simulate(p, species_counts, reactions, coordinate_list);

end