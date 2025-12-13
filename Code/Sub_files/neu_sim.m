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
        options.sites = 2;

        % Network parameters
        options.rings = 2;
        options.subcycles_per_ring = [3,3];
        options.prop_forms = [1,0];
        options.formation_type = ["split", "split"];
        options.reac_rate = [0.01,0.01];
        options.must_adsorb = [1,1];
        options.fac_rings = [0,0];

        options.flow_rate = 1;
        options.in_rate = 10;   % THIS SEEMS GOOD, KEEP IT FOR NOW!
        options.out_rate = 1;
        options.fierce_reaction_rate = 0.01;
        options.prop_formation_rate = 0.0005;   % SEEMS GOOD?
        options.prop_funnel_rate = 150;
        options.adsorb_rate = 0.01;
        options.independence_disadvantage = [Inf];

        % Starting condition parameters
        options.default_concentration = 0;
        options.food_concentration = 1000;
        options.seed_concentration = [25,25];
        options.site_concentration = 500;
        options.seed_locations = [1,1];
        options.seed_state = ["_diff", "_diff"];

        % Simulation parameters
        options.t_max = 1000;
        options.disturb_freq = 0;
        options.sample_number = 10;     % Try messing with this next!!!

        % Debug
        options.introspection = true;

    end
    
    rng("shuffle");     % Ensure stochastic sims are different every time
    
    % Simply convert the options into the famed parameter struct p

    p = options;

    % Initialize the output file

    %matrixname = "results_" + chtc_batch + "_" + run + "_" + separration + "_" + cycles + "_" + flipseeds + "_" + independence_disadvantage + "_" + disturb_freq + ".txt";
    %writematrix(["Relative_concentration", "lrc_variance_prop", "lrc_variance_noprop", "P_extinct", "NP_extinct", "P_tot_extinct", "NP_tot_extinct", "separration", "cycles", "independence_disadvantage", "disturb_freq", "flipseeds"], matrixname);
    matrixname = "results_" + p.prop_formation_rate + "_" + p.separation_distance + "_" + p.sites + "_" + randi([1,9999999]) + ".txt";
    writematrix(["Relative_concentration", "separration", "sites", "prop_formation_rate"], matrixname);


    % Do the simulation
    %p.disperse_frequency = p.t_max / (p.sample_number/10);
    p.disperse_frequency = p.sample_number;
    I = struct;
    O = struct;
    tic;

    fprintf("Initializing reactions >w<\n");
    [I] = init_reactions(p, I);
    fprintf("Initializing world >w<\n");
    [I] = init_world(p, I);
    O.build_time = toc;

    tic;
    [O, I] = simulate(p, I, O);
    O.run_time = toc;

    [O] = summarize(p, I, O);

    % Write data to output file
    %writematrix(["Relative_concentration", "separration", "sites", "prop_formation_rate"]);
    current_data = [O.relative_concentration, p.separation_distance, p.sites, p.prop_formation_rate];
    fid = fopen(matrixname, 'a+'); 
    fprintf(fid, "%d,%d,%d,%d\n", current_data);
    fclose(fid); 
    fprintf("Data saved. Happy days! >w<\n")

end