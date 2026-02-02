function [] = precalculate_occupation_p(p, I, options)

    % Precalculate and save probability distribution so this only has to be
    % calculated ONCE

arguments (Input)
    p;
    I;
    options.prop = false;
end

    R = I.R;
    prob_cloud = cell(1,size(I.all_coordinates,2));

    parfor nZero = 1:size(I.all_coordinates,2)
        missing_prob = 1;
        for nOne = 1:size(I.all_coordinates,2)

           if options.prop == false, t = R; else, t = floor(R/sqrt(max(p.subcycles_per_ring))); end  % Revise t to reflect diffusion rates (or make 2 versions of prob_cloud)
            n01 = I.all_coordinates{nZero}{1}; n02 = I.all_coordinates{nZero}{2};
            n1 = I.all_coordinates{nOne}{1}; n2 = I.all_coordinates{nOne}{2};

            P = single(round(occupation_p(p.shape, R, 1, t, n01, n02, n1, n2), 4));
            % ~~~ The following is a temporary approximation ~~~
            % ~~~ Figure out why probs are going negative ~~~~~~
            if P < 0, P = 0; end

            missing_prob = missing_prob - P;
            prob_cloud{nZero}{1,end+1} = P;
            
        end
        %prob_cloud{nZero}{nZero} = prob_cloud{nZero}{nZero} + missing_prob; % Make total probs = 1. Higher probability of a particle staying in original site.
    end
    
    if options.prop == false
        save(I.map_name, "prob_cloud", "-v7.3");
    else
        prop_prob_cloud = prob_cloud;
        save(I.prop_map_name, "prop_prob_cloud", "-v7.3");
    end

end