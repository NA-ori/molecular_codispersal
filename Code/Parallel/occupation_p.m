function [P] = occupation_p(shape, R, q, t, n01, n02, n1, n2)

    % Calculate the occupation probability of a site on a hexagonal lattice

    % Note that the parallelogram form is defective
    
    arguments (Input)
        shape;
        R;          % Circumradius/side length
        q;          % Probability of movement
        t;          % Number of movement steps
        n01; n02;   % Starting coords
        n1; n2;     % Ending coords
    end
    
    if shape == "hex"

        omega = 3*R*(R+1) + 1;
        summation = 0;
        for m = 0:(R-1)
            for s = 0:((3*m)+2)
                k1 = ((R+1)*s) - m + R;
                k2 = (((3*R)+1)*m) - (R*s) + (2*R) + 1;
                current_term = cos( (n1-n01)*((2*k1*pi)/omega) + (n2-n02)*((2*k2*pi)/omega) ) * ...
                    ( 1-q + (q/3)*cos(((k1-k2)*2*pi)/omega) + (q/3)*cos((2*k1*pi)/omega) + (q/3)*cos((2*k2*pi)/omega) )^t;
                summation = summation + current_term;
            end
        end
        P = (1/omega) + (2/omega)*summation;

    elseif shape == "parallelogram"

        omega = R^2;
        summation = 0;
        for k1 = 0:(R-1)
            for k2 = 0:(R-1)
                current_term = cos( (n1-n01)*((2*k1*pi)/omega) + (n2-n02)*((2*k2*pi)/omega) ) * ...
                    ( 1-q + (q/3)*cos(((k1-k2)*2*pi)/omega) + (q/3)*cos((2*k1*pi)/omega) + (q/3)*cos((2*k2*pi)/omega) )^t;
                summation = summation + current_term;
            end
        end
        P = (1/omega)*summation;

    end
end