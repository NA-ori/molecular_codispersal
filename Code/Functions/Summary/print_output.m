function [] = print_output(p, O, options)

    % Print output of one simulation; output files are stitchable into a
    % csv and can then be graphed and such

    arguments (Input)
        p; O;
        options.file_prefix = "";
        options.file_suffix = "";
        options.random_id = true;
        options.output_names = ["End_time", "Runtime", "Max_cycles", "Facultative_modifier", "Relative_concentration", "Approx_unoccupied_pixels", "Sites_occupied"];
        options.output_values = [O.end_time, O.run_time, max(p.subcycles_per_ring), max(p.independence_disadvantage), O.relative_concentration, O.approx_unoccupied_pixels, O.site_occupation(1)];
    end

    % Name output file
    matrixname = options.file_prefix + "results" + options.file_suffix + "_" + O.sim_id;
    matrixname = matrixname + ".txt";

    % ~~~~~~~~~~~~~~~~~~~~~ Write data ~~~~~~~~~~~~~~~~~~~~~~~~~

    fields = fieldnames(p);
    fid = fopen(matrixname, 'a+');

    % Add labels
    for param = 1:length(fields)
        if isscalar(p.(fields{param})) && ~isstring(p.(fields{param}))
            fprintf(fid, "%s,", fields{param});
        end
    end

    for output = 1:length(options.output_names)
        if output == length(options.output_names)
            fprintf(fid, "%s\n", options.output_names(output));
        else
            fprintf(fid, "%s,", options.output_names(output));
        end
    end

    % Add values
    for param = 1:length(fields)
        if isscalar(p.(fields{param})) && ~isstring(p.(fields{param}))
            fprintf(fid, "%d,", p.(fields{param}));
        end
    end

    for output = 1:length(options.output_values)
        if output == length(options.output_values)
            fprintf(fid, "%d\n", options.output_values(output));
        else
            fprintf(fid, "%d,", options.output_values(output));
        end
    end

    fclose(fid);
    fprintf("Data saved. Happy days! >w<\n");

end