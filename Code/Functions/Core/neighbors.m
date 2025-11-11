function [friends] = neighbors(x,y,z,sep)

    % Calculate the coordinates of all neighboring pixels at a set distance

    arguments (Input)
        x;
        y;
        z;
        sep;
    end

    friends = cell(1,6);

    friends{1}{1} = x+0; friends{1}{2} = y+sep; friends{1}{3} = z-sep;
    friends{2}{1} = x+sep; friends{2}{2} = y+0; friends{2}{3} = z-sep;
    friends{3}{1} = x+sep; friends{3}{2} = y-sep; friends{3}{3} = z-0;
    friends{4}{1} = x+0; friends{4}{2} = y-sep; friends{4}{3} = z+sep;
    friends{5}{1} = x-sep; friends{5}{2} = y+0; friends{5}{3} = z+sep;
    friends{6}{1} = x-sep; friends{6}{2} = y+sep; friends{6}{3} = z-0;

end