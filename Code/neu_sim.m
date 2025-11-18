function [reactions, ring_list, base_species, species_counts, member_species_record, current_reaction_propensities, current_chemical_counts, concentration_tracker, p] = neu_sim(options)

% Smarter version of the original CHTC sim.
% Reduces hard-coding and therefore hopefully mistakes.
% Also, put the whole thing in a function so every value is easily
% changeable.

% Ideally near-everything should be alterable just by changing the function
% inputs.

    arguments (Input)

        % World parameters
        options.shape {mustBeMember(options.shape,{'parallelogram','hex'})} = "hex";
        options.separation_distance {mustBeInteger,mustBePositive} = 2;
        options.sites {mustBeInteger,mustBePositive} = 2;

        % Network parameters
        options.rings {mustBeInteger,mustBePositive} = 2;
        options.subcycles_per_ring = [3,3];
        options.prop_forms = [1,0];
        options.formation_type = ["split", "split"];
        options.reac_rate = [0.01,0.01];
        options.must_adsorb = [1,1];
        options.fac_rings = [0,0];

        options.flow_rate = 1;
        options.in_rate = 0.1;
        options.out_rate = 1;
        options.fierce_reaction_rate = 0.01;
        options.prop_formation_rate = 0.01/5;
        options.prop_funnel_rate = 0.1;
        options.adsorb_rate = 0.01;
        options.independence_disadvantage = [Inf];

        % Starting condition parameters
        options.default_concentration = 0;
        options.food_concentration = 500;
        options.seed_concentration = [25,25];
        options.site_concentration = 500;
        options.seed_locations = [1,1];
        options.seed_state = ["_diff", "_diff"];

        % Simulation parameters
        options.t_max = 500;
        options.disturb_freq = 0;
        options.sample_number = 1000;

        % Debug
        options.introspection = true;

    end
    
    rng("shuffle");     % Ensure stochastic sims are different every time
    
    % Simply convert the options into the famed parameter struct p

    p = options;
    p.disperse_frequency = p.t_max / p.sample_number;

    O = struct;     % Output value structure
    I = struct;     % Internal value structure

    fprintf("Initializing reactions >w<\n");
    [reactions, ring_list, member_species_record, base_species] = init_reactions(p);
    fprintf("Initializing world >w<\n");
    [species_counts, coordinate_list, prob_cloud, prop_prob_cloud] = init_world(p, reactions, ring_list, base_species, member_species_record);

    [species_counts, time, current_reaction_propensities, current_chemical_counts, concentration_tracker] = simulate(p, base_species, species_counts, reactions, coordinate_list, prob_cloud, prop_prob_cloud);

end