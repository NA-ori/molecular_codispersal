function [friends] = neighbors(x,y,z,sep,I)

    % Calculate the coordinates of all neighboring pixels at a set distance

    arguments (Input)
        x;
        y;
        z;
        sep;
        I;
    end

    friends = cell(1,6);

    friends{1} = [x, y+sep, z-sep];
    friends{2} = [x+sep, y, z-sep];
    friends{3} = [x+sep, y-sep, z];
    friends{4} = [x, y-sep, z+sep];
    friends{5} = [x-sep, y, z+sep];
    friends{6} = [x-sep, y+sep, z];

    % Determine whether you need to wrap around to the other side
    for friend = 1:size(friends,2)
        if any(abs(friends{friend}) > I.R)
            % Find shortest distance to a reflected center
            subtracted = cell(1,6);
            distances = zeros(1,6);
            for centre = 1:size(I.centres,2)
                subtracted{centre} = friends{friend} - I.centres{centre};
                distances(centre) = max(abs(subtracted{centre}));
            end
            [~, dex] = min(distances);
            friends{friend} = subtracted{dex};
        end
    end


    % friends{1}{1} = x; friends{1}{2} = y+sep; friends{1}{3} = z-sep;
    % friends{2}{1} = x+sep; friends{2}{2} = y; friends{2}{3} = z-sep;
    % friends{3}{1} = x+sep; friends{3}{2} = y-sep; friends{3}{3} = z;
    % friends{4}{1} = x; friends{4}{2} = y-sep; friends{4}{3} = z+sep;
    % friends{5}{1} = x-sep; friends{5}{2} = y; friends{5}{3} = z+sep;
    % friends{6}{1} = x-sep; friends{6}{2} = y+sep; friends{6}{3} = z;

end