function [] = neu_sim(options)

% Smarter version of the original CHTC sim.
% Reduces hard-coding and therefore hopefully mistakes.
% Also, put the whole thing in a function so every value is easily
% changeable.

% Ideally near-everything should be alterable just by changing the function
% inputs.

    arguments (Input)
        % Inputs are modified into lists to allow more flexibility
        options.dim = 0;
        options.t_max = 10000;
        options.default_concentration = 0;
        options.food_concentration = [500];
        options.seed_concentration = [25,25];
        options.site_concentration = [500];
        options.independence_disadvantage = [Inf];
        options.flow_rate = 0.1;
        options.in_rate = 0.1;
        options.out_rate = 0.1;
        options.adsorp_rate = 0.01;
        options.fierce_reaction_rate = 0.01;
        options.prop_form_rate = 0.01/5;
        options.food_set = [];
        options.site_set = [];
        options.seed_set = [];

    end
    
    rng("shuffle");     % Ensure stochastics sims are different every time
    
    % Simply convert the options into the famed parameter struct p

    p = options;

end