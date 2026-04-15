function [O] = approx_unoccupied_pix(p, I, O)

    arguments (Input)
        p;
        I;
        O;
    end

    sites_index = strcmp(O.global_species_counts(1,:), "site");

    sites_count = O.global_species_counts{end, sites_index};

    O.approx_unoccupied_pixels = str2double(sites_count) / p.site_concentration;

end