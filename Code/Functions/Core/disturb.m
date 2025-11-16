function [current_chemical_counts] = disturb(p, current_chemical_counts, disturbed_coord, base_species)

    arguments (Input)
        p;
        current_chemical_counts;
        disturbed_coord;
        base_species;
    end
    
    site_mask = strcmp(base_species, "site");
    current_chemical_counts{1,disturbed_coord}(2,:) = {0};
    current_chemical_counts{1,disturbed_coord}{2,site_mask} = p.site_concentration;
end