function [reactions, rates, tags, ring_list] = new_init_reactions(p, options)
    % Test a new reaction format.
    % This will be based entirely on cell arrays containing numbers. No
    % string manipulation necessary.
    % This will maybe make things faster and easier? I'll see as I go I
    % guess

    arguments (Input)
        p;
        options.rings {mustBeInteger} = 2;
        options.subcycles_per_ring = [3,3];
        options.prop_forms = [0,1];
        options.formation_type = ["cheater", "cheater"];
        options.reac_rate = [0.01,0.01];
        options.must_adsorb = [1,1];
        options.fac_rings = [0,0];
    end

    rings = options.rings; subcycles_per_ring = options.subcycles_per_ring; prop_forms = options.prop_forms; formation_type = options.formation_type;
    reac_rate = options.reac_rate; must_adsorb = options.must_adsorb; fac_rings = options.fac_rings

    % The basic structure:
    % Reactions are stored in 1 cell array of rows = reactions and cols =
    % species. Species in a row are assigned 1 or -1 etc depending on
    % whether they are products or reactants. Most reverse reactions can
    % easily be created by inverting these numbers.
    % Reaction rates and identifying tags are stored in other lists. Add a
    % function to easily bind these lists together for readability.

    reactions = {};
    rates = {};
    tags = {};

    % Keep track of the 'species' (ie columns) that are present in each
    % ring bc sorting the results out would otherwise be impossible

    ring_list = cell(1, options.rings);

    % Initialize food and waste in positions 1 and 2.

    % Make the rings.
    for current_ring = 1:rings
        

    end

end