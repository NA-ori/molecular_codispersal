function [I] = init_reactions(p, I)
    %rings, subcycles_per_ring, prop_forms, formation_type, reac_rate, must_adsorb, fac_rings)

% Initialize a basic set of reactions.

    arguments (Input)
        p;
        I;
    end

    % Ring generation method written by A. Khanov

    reactions = {};
    ring_list = cell(1,p.rings);
    base_species = [];
    ring_species = []; % will contain central ring species ONLY
    member_species_record = cell(1,p.rings);

    for current_ring = 1:p.rings
        member_species_list = [];
        if p.must_adsorb(current_ring) == 1
            adsorb_tag = "_ad";
        else
            adsorb_tag = "_diff";
        end

        for current_sub = 1:p.subcycles_per_ring(current_ring)
            curren_reac_rate = p.reac_rate(current_ring);

            % generate a member species and intermediates:
            member_species = "sp_r" + string(current_ring) + "_" + string(current_sub);
            member_species_intermediate = "sp_r" + string(current_ring) + "_" + string(current_sub) + "_i";
            member_species_waste = "sp_r" + string(current_ring) + "_" + string(current_sub) + "_w";

            % link subcycles through waste:
            mutualist_waste = "";
            if current_sub == p.subcycles_per_ring
                mutualist_waste = "sp_r" + string(current_ring) + "_" + string(1) + "_w";
            else
                mutualist_waste = "sp_r" + string(current_ring) + "_" + string(current_sub + 1) + "_w";
            end

            member_species_list = [member_species_list, member_species];

            ring_species = [ring_species, member_species, member_species_intermediate, member_species_waste];
            ring_list{1, current_ring} = [ring_list{1, current_ring}, member_species, member_species_intermediate, member_species_waste];
            member_species_record{1, current_ring} = member_species_list;

            % generate the list of reactions for the subcycle
            reactions{end + 1,1} = {curren_reac_rate, {append(member_species, adsorb_tag), "F", "site"}, {1, 1, 1}, {append(member_species_intermediate, adsorb_tag), append(member_species_waste, adsorb_tag)}, {1, 1}, "autocat"};
            reactions{end + 1,1} = {curren_reac_rate, {append(member_species_intermediate, adsorb_tag), append(member_species_waste, adsorb_tag)}, {1, 1}, {append(member_species, adsorb_tag), "F", "site"}, {1, 1, 1}, "autocat"};
            reactions{end + 1,1} = {curren_reac_rate, {append(member_species_intermediate, adsorb_tag), append(mutualist_waste, adsorb_tag)}, {1, 1}, {append(member_species, adsorb_tag)}, {2}, "autocat"};
            reactions{end + 1,1} = {curren_reac_rate, {append(member_species, adsorb_tag)}, {2}, {append(member_species_intermediate, adsorb_tag), append(mutualist_waste, adsorb_tag)}, {1, 1}, "autocat"};
            if p.fac_rings(current_ring) == 1
                reactions{end + 1,1} = {curren_reac_rate/p.independence_disadvantage, {append(member_species_intermediate, adsorb_tag), "F"}, {1, 1}, {append(member_species, adsorb_tag)}, {2}, "autocat"};
                reactions{end + 1,1} = {curren_reac_rate/p.independence_disadvantage, {append(member_species, adsorb_tag)}, {2}, {append(member_species_intermediate, adsorb_tag), "F"}, {1, 1}, "autocat"};
            end

        end


        if p.prop_forms(current_ring) == 1
            % Add propagule formation here if necessary
            % Two methods for doing this

            % 'cheater' makes a higher-order reaction if necessary and
            % then the propensity calculating code will base the propensity
            % on the counts of only two randomly chosen reactants
            % This is an egregious violation of the gillespie algorithm but it
            % should perhaps make prop formation rates comparable for the
            % purposes of isolating the role of propagules

            % 'split' splits everything into max second-order reactions.
            % Much more sensible at the cost of being a pain.

            reactants_list = {};
            for i = 1:length(member_species_list)
                reactants_list{end+1} = append(member_species_list(i), adsorb_tag);
            end

            free_reactants_list = {};
            for i = 1:length(member_species_list)
                free_reactants_list{end+1} = append(member_species_list(i), "_diff");
            end            

            reac_stoich = num2cell(ones(1, length(member_species_list)));

            if p.formation_type(current_ring) == "cheater"
                
                reactions{end + 1,1} = {p.prop_formation_rate, reactants_list, reac_stoich, {("prop_" + string(current_ring)), "site"}, {1, length(member_species_list)}, "prop_form"};
                reactions{end + 1,1} = {p.prop_formation_rate, {("prop_" + string(current_ring))}, {1}, free_reactants_list, reac_stoich, "prop_break"};

            elseif p.formation_type(current_ring) == "split"

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Cheating for now. Need to add auto-generation!%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if p.subcycles_per_ring(current_ring) == 2
                    reactions{end + 1,1} = {p.prop_formation_rate, reactants_list, reac_stoich, {("prop_" + string(current_ring)), "site"}, {1, 2}, "prop_form"};
                    reactions{end + 1,1} = {p.prop_break_rate, {("prop_" + string(current_ring))}, {1}, free_reactants_list, reac_stoich, "prop_break"};

                elseif p.subcycles_per_ring(current_ring) == 3
                    reactions{end + 1,1} = {p.prop_formation_rate, {"sp_r1_1_ad", "sp_r1_2_ad"}, {1,1}, {"sp_r1_1_2_prop"}, {1}, "init_prop"};
                    reactions{end + 1,1} = {p.prop_formation_rate, {"sp_r1_1_2_prop"}, {1}, {"sp_r1_1_ad", "sp_r1_2_ad"}, {1,1}, "init_prop"};
                    reactions{end + 1,1} = {p.prop_formation_rate, {"sp_r1_1_ad", "sp_r1_3_ad"}, {1,1}, {"sp_r1_1_3_prop"}, {1}, "init_prop"};
                    reactions{end + 1,1} = {p.prop_formation_rate, {"sp_r1_1_3_prop"}, {1}, {"sp_r1_1_ad", "sp_r1_3_ad"}, {1,1}, "init_prop"};
                    reactions{end + 1,1} = {p.prop_formation_rate, {"sp_r1_2_ad", "sp_r1_3_ad"}, {1,1}, {"sp_r1_2_3_prop"}, {1}, "init_prop"};
                    reactions{end + 1,1} = {p.prop_formation_rate, {"sp_r1_2_3_prop"}, {1}, {"sp_r1_2_ad", "sp_r1_3_ad"}, {1,1}, "init_prop"};
                    reactions{end + 1,1} = {p.prop_funnel_rate, {"sp_r1_1_2_prop", "sp_r1_3_ad"}, {1,1}, {("prop_" + string(current_ring)), "site"}, {1,3}, "end_prop"};
                    reactions{end + 1,1} = {p.prop_funnel_rate, {"sp_r1_1_3_prop", "sp_r1_2_ad"}, {1,1}, {("prop_" + string(current_ring)), "site"}, {1,3}, "end_prop"};
                    reactions{end + 1,1} = {p.prop_funnel_rate, {"sp_r1_2_3_prop", "sp_r1_1_ad"}, {1,1}, {("prop_" + string(current_ring)), "site"}, {1,3}, "end_prop"};
                    reactions{end + 1,1} = {p.prop_break_rate, {("prop_" + string(current_ring))}, {1}, {"sp_r1_1_diff", "sp_r1_2_diff", "sp_r1_3_diff"}, {1,1,1}, "prop_break"};

                elseif p.subcycles_per_ring(current_ring) == 4
                    reactions{end + 1,1} = {p.prop_formation_rate, {"sp_r1_1_ad", "sp_r1_2_ad"}, {1,1}, {"sp_r1_1_2_prop"}, {1}, "init_prop"};
                    reactions{end + 1,1} = {p.prop_formation_rate, {"sp_r1_1_2_prop"}, {1}, {"sp_r1_1_ad", "sp_r1_2_ad"}, {1,1}, "init_prop"};

                    reactions{end + 1,1} = {p.prop_formation_rate, {"sp_r1_1_ad", "sp_r1_3_ad"}, {1,1}, {"sp_r1_1_3_prop"}, {1}, "init_prop"};
                    reactions{end + 1,1} = {p.prop_formation_rate, {"sp_r1_1_3_prop"}, {1}, {"sp_r1_1_ad", "sp_r1_3_ad"}, {1,1}, "init_prop"};

                    reactions{end + 1,1} = {p.prop_formation_rate, {"sp_r1_1_ad", "sp_r1_4_ad"}, {1,1}, {"sp_r1_1_4_prop"}, {1}, "init_prop"};
                    reactions{end + 1,1} = {p.prop_formation_rate, {"sp_r1_1_4_prop"}, {1}, {"sp_r1_1_ad", "sp_r1_4_ad"}, {1,1}, "init_prop"};                    

                    reactions{end + 1,1} = {p.prop_formation_rate, {"sp_r1_2_ad", "sp_r1_3_ad"}, {1,1}, {"sp_r1_2_3_prop"}, {1}, "init_prop"};
                    reactions{end + 1,1} = {p.prop_formation_rate, {"sp_r1_2_3_prop"}, {1}, {"sp_r1_2_ad", "sp_r1_3_ad"}, {1,1}, "init_prop"};  

                    reactions{end + 1,1} = {p.prop_formation_rate, {"sp_r1_2_ad", "sp_r1_4_ad"}, {1,1}, {"sp_r1_2_4_prop"}, {1}, "init_prop"};
                    reactions{end + 1,1} = {p.prop_formation_rate, {"sp_r1_2_4_prop"}, {1}, {"sp_r1_2_ad", "sp_r1_4_ad"}, {1,1}, "init_prop"};

                    reactions{end + 1,1} = {p.prop_formation_rate, {"sp_r1_3_ad", "sp_r1_4_ad"}, {1,1}, {"sp_r1_3_4_prop"}, {1}, "init_prop"};
                    reactions{end + 1,1} = {p.prop_formation_rate, {"sp_r1_3_4_prop"}, {1}, {"sp_r1_3_ad", "sp_r1_4_ad"}, {1,1}, "init_prop"};

                    reactions{end + 1,1} = {p.prop_funnel_rate, {"sp_r1_1_2_prop", "sp_r1_3_4_prop"}, {1,1}, {("prop_" + string(current_ring)), "site"}, {1,4}, "end_prop"};
                    reactions{end + 1,1} = {p.prop_funnel_rate, {"sp_r1_1_3_prop", "sp_r1_2_4_prop"}, {1,1}, {("prop_" + string(current_ring)), "site"}, {1,4}, "end_prop"};
                    reactions{end + 1,1} = {p.prop_funnel_rate, {"sp_r1_1_4_prop", "sp_r1_2_3_prop"}, {1,1}, {("prop_" + string(current_ring)), "site"}, {1,4}, "end_prop"};
                    reactions{end + 1,1} = {p.prop_break_rate, {("prop_" + string(current_ring))}, {1}, {"sp_r1_1_diff", "sp_r1_2_diff", "sp_r1_3_diff", "sp_r1_4_diff"}, {1,1,1,1}, "prop_break"};

                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % End Cheating %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end

        end

    end

    % Now add adsorption/desorption reactions
    for spec = 1:length(ring_species)
        reactions{end + 1,1} = {p.adsorb_rate, {(ring_species(spec) + "_diff"), "site"}, {1,1}, {(ring_species(spec) + "_ad")}, {1}, "adsorb"};
        reactions{end + 1,1} = {p.adsorb_rate, {(ring_species(spec) + "_ad")}, {1}, {(ring_species(spec) + "_diff"), "site"}, {1,1}, "desorb"};
    end

    % Add inflow reaction
    if p.chemostat == false
        reactions{end + 1,1} = {p.in_rate, {"F"}, {0}, {"F"}, {1}, "inflow"};
    end

    % Retrieve a list of all species in the network
    for x = 1:size(reactions, 1)  % get a list of all the species from the basic reactions
        for y = 1:size(reactions{x}{2}, 2)
            base_species = [base_species, reactions{x}{2}{y}];
        end
        for z = 1:size(reactions{x}{4}, 2)
            base_species = [base_species, reactions{x}{4}{z}];
        end
    end
    base_species = unique(base_species);    % remove all non-unique species

    % Add everything to the internal struct

    I.reactions = reactions;
    I.ring_list = ring_list;
    I.member_species_record = member_species_record;
    I.base_species = base_species;

end