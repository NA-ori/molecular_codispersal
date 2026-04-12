function [current_chemical_counts] = disturb(p, I, current_chemical_counts, disturbed_coord)

    arguments (Input)
        p;
        I;
        current_chemical_counts;
        disturbed_coord;
    end
    
    site_mask = strcmp(I.base_species, "site");
    food_mask = strcmp(I.base_species, "F");

    current_chemical_counts(disturbed_coord,:) = zeros(1, length(I.base_species));
    current_chemical_counts(disturbed_coord,site_mask) = p.site_concentration;
    current_chemical_counts(disturbed_coord,food_mask) = p.food_concentration;

    % current_chemical_counts{1,disturbed_coord}(2,:) = {0};
    % current_chemical_counts{1,disturbed_coord}{2,site_mask} = p.site_concentration;
    % current_chemical_counts{1,disturbed_coord}{2,food_mask} = p.food_concentration;
end