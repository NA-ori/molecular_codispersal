function [] = print_reactions(reactions)

    % Printing to test function functionality :3
    for current_row = 1:size(reactions, 1)
        fprintf("{" + string(reactions{current_row}{1}) + " ")
        fprintf("{ ")
        for column = 1:size(reactions{current_row}{2}, 2)
            fprintf(string(reactions{current_row}{2}{column}) + " ")
        end
        fprintf("} ")
        fprintf("{ ")
        for column = 1:size(reactions{current_row}{3}, 2)
            fprintf(string(reactions{current_row}{3}{column}) + " ")
        end
        fprintf("} ")
        fprintf("{ ")
        for column = 1:size(reactions{current_row}{4}, 2)
            fprintf(string(reactions{current_row}{4}{column}) + " ")
        end
        fprintf("} ")
        fprintf("{ ")
        for column = 1:size(reactions{current_row}{5}, 2)
            fprintf(string(reactions{current_row}{5}{column}) + " ")
        end
        fprintf("} ")
        fprintf("{ ")
        for column = 1:size(reactions{current_row}{6}, 2)
            fprintf(reactions{current_row}{6}{column} + "}")
        end
        fprintf("} ")
        fprintf("\n")
    end
end