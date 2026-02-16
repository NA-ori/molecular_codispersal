function [p, I, O] = neu_sim(options)

% Smarter version of the original CHTC sim.
% Reduces hard-coding and therefore hopefully mistakes.
% Also, put the whole thing in a function so every value is easily
% changeable.

% Ideally near-everything should be alterable just by changing the function
% inputs.

    arguments (Input)

        options.shape = "hex";
        options.separation_distance = 7;
        options.sites = 1;

        % Network parameters
        options.rings = 2;
        options.subcycles_per_ring = [3,3];
        options.prop_forms = [1,0];
        options.formation_type = ["split", "split"];
        options.reac_rate = [0.01,0.01];
        options.must_adsorb = [1,1];
        options.fac_rings = [0,0];

        options.flow_rate = 0.6;
        options.in_rate = 50;
        options.out_rate = 0.1;
        options.fierce_reaction_rate = 0.01;
        options.prop_formation_rate = 0.01 / 5;
        options.prop_funnel_rate = 1500;
        options.prop_break_rate = 0.01 / 5;
        options.adsorb_rate = 0.01;
        options.independence_disadvantage = [Inf];

        % Starting condition parameters
        options.seed_start = "origin";
        options.default_concentration = 0;
        options.food_concentration = 1000;
        options.seed_concentration = [25,25];
        options.site_concentration = 500;
        options.seed_locations = [159,173];
        options.seed_state = ["_ad", "_ad"];

        % Simulation parameters
        options.t_max = 5000;
        options.disturb_freq = 0;
        options.sample_number = 50;
        options.recording_freq = 50;
        options.diffusion_divisor = 0.5;
        options.chemostat = true;

        % Config
        options.introspection = true;
        options.map_prefix = "";

    end
    
    rng("shuffle");     % Ensure stochastic sims are different every time
    
    % Simply convert the options into the famed parameter struct p

    p = options;

    % Do the simulation
    p.disperse_frequency = p.sample_number;
    I = struct;
    O = struct;
    tic;

    fprintf("Initializing reactions >w<\n");
    [I] = init_reactions(p, I);
    fprintf("Initializing world >w<\n");
    [p, I] = linear_init_world(p, I);
    O.build_time = toc;

    tic;
    [O, I] = linear_simulate(p, I, O);
    O.run_time = toc;

    [O] = summarize(p, I, O);
    [O] = approx_unoccupied_pix(p, I, O);

    print_output(p,O);

    % Write data to output file
    % matrixname = "results_" + p.separation_distance + "_" + p.prop_break_rate + "_" + randi([1,9999999]) + ".txt";
    % writematrix(["Relative_concentration", "Approx_unoccupied_pixels", "separation_distance"], matrixname);
    % 
    % current_data = [O.relative_concentration, O.approx_unoccupied_pixels, p.separation_distance];
    % fid = fopen(matrixname, 'a+'); 
    % fprintf(fid, "%d,%d,%d\n", current_data);
    % fclose(fid); 
    % fprintf("Data saved. Happy days! >w<\n");

end
