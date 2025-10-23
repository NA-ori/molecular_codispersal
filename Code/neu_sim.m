function [] = neu_sim(dim, t_max, default_concentration, food_concentration, seed_concentration, site_concentration, independence_disadvantage, flow_rate, in_rate, ...
    out_rate, adsorp_rate, fierce_reaction_rate, prop_form_rate, food_set, site_set, seed_set)

% Smarter version of the original CHTC sim.
% Reduces hard-coding and therefore hopefully mistakes.
% Also, put the whole thing in a function so every value is easily
% changeable.

% Ideally near-everything should be alterable just by changing the function
% inputs.

    arguments (Input)
        % Inputs are modified into lists to allow more flexibility
        dim = 0;
        t_max = 10000;
        default_concentration = 0;
        food_concentration = [500];
        seed_concentration = [25,25];
        site_concentration = [500];
        independence_disadvantage = [Inf];
        flow_rate = 0.1;
        in_rate = flow_rate;
        out_rate = flow_rate;
        adsorp_rate = 0.01;
        fierce_reaction_rate = 0.01;
        prop_form_rate = fierce_reaction_rate/5;
        food_set = [];
        site_set = [];
        seed_set = [];

    end
    
    rng("shuffle");     % Ensure stochastics sims are different every time
    
    % Simply stick everything into a structure

    p = struct;
    p.dim = dim;
    p.t_max = t_max;
    p.default_concentration = default_concentration;
    p.food_concentration = food_concentration;
    p.seed_concentration = seed_concentration;
    p.sites = site_concentration;
    p.independence_disadvantage = independence_disadvantage;
    p.flow_rate = flow_rate;
    p.in_rate = in_rate;
    p.out_rate = out_rate;
    p.adsorp_rate = adsorp_rate;
    p.fierce_reaction_rate = fierce_reaction_rate;
    p.prop_form_rate = prop_form_rate;
    p.food_set = food_set;
    p.site_set = site_set;
    p.seed_set = seed_set;

end