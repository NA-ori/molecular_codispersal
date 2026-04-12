function [current_chemical_counts] = disperse(I, num_coords, diff_ind, current_chemical_counts, cloud)

    arguments (Input)
        I;
        num_coords;
        diff_ind;
        current_chemical_counts;
        cloud;
    end

    % Calculate where everything goes
    update_mask = zeros(num_coords,length(I.base_species));
    for coord = 1:num_coords
        for spec = 1:length(diff_ind)
            species = diff_ind(spec);
            if current_chemical_counts(coord,species) > 0

                sample_array = randsample(1:num_coords, current_chemical_counts(coord,species), true, cloud(coord,:));
                [counts,ordering] = groupcounts(sample_array');

                for d = 1:length(ordering)
                    update_mask(ordering(d),species) = update_mask(ordering(d),species) + counts(d);
                    update_mask(coord,species) = update_mask(coord,species) - counts(d);
                end

            end
        end
    end
    clear coord; clear species; clear spec; clear sample_array; clear d; clear counts; clear ordering;

    % Move things around
    current_chemical_counts = current_chemical_counts + update_mask;
    clear update_mask;

end








% function [current_chemical_counts] = disperse(p, update_mask, diff_mask, current_chemical_counts, prob_cloud, prop_prob_cloud)
% 
%     arguments (Input)
%         p;
%         update_mask;
%         diff_mask;
%         current_chemical_counts;
%         %coordinate_list;
%         prob_cloud;
%         prop_prob_cloud;
%     end
% 
%     % Calculate where everything goes
%     for coord = 1:size(current_chemical_counts,2)
%         for species = 1:size(current_chemical_counts{1,coord},2)
%             if (diff_mask(species) == 1) && (current_chemical_counts{1,coord}{2,species} > 0)
%                 % Must update prop_prob_cloud so it accounts for potential
%                 % multiple rings with different numbers of subcycles
% 
%                 particle_number = sum(randsample(0:1, current_chemical_counts{1,coord}{2,species}, true, [1-p.diffusion_divisor,p.diffusion_divisor]));
% 
%                 if contains(current_chemical_counts{1,coord}{1,species}, "prop")
%                      sample_array = randsample(1:size(current_chemical_counts,2), particle_number, true, cell2mat(prop_prob_cloud{coord}));
%                 else
%                     sample_array = randsample(1:size(current_chemical_counts,2), particle_number, true, cell2mat(prob_cloud{coord}));
%                 end
% 
%                 [counts, ordering] = groupcounts(sample_array');
% 
%                 for d = 1:length(ordering)
%                     update_mask{ordering(d)}{2,species} = update_mask{ordering(d)}{2,species} + counts(d);
%                     update_mask{coord}{2,species} = update_mask{coord}{2,species} - counts(d);
%                 end
%             end
%         end
%     end
% 
%     % Move things around
%     for coord = 1:size(current_chemical_counts,2)
%         for species = 1:size(current_chemical_counts{1,coord},2)
%             current_chemical_counts{1,coord}{2,species} = current_chemical_counts{1,coord}{2,species} + update_mask{1,coord}{2,species};
%         end
%     end
% 
% end