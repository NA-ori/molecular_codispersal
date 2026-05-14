function [] = print_output(p, I, O, options)

    % Print output of one simulation; output files are stitchable into a
    % csv and can then be graphed and such

    arguments (Input)
        p; I; O;
        options.file_prefix = "";
        options.file_suffix = "";
        options.random_id = true;
        options.output_names = ["Fastest_reac_rate", "End_time", "Runtime", "Max_cycles", "Facultative_modifier", "Relative_concentration", "Approx_unoccupied_pixels", "Sites_occupied", "ID", "ad_A_extinct", "ad_NA_extinct", "ad_none_extinct", "ad_all_extinct", "tot_A_extinct", "tot_NA_extinct", "tot_all_extinct", "tot_none_extinct"];
        options.output_values = [max(p.reac_rate), O.end_time, O.run_time, max(p.subcycles_per_ring), max(p.independence_disadvantage), O.relative_concentration, O.approx_unoccupied_pixels, O.site_occupation(1), O.sim_id, O.adsorbed_A_extinct, O.adsorbed_NA_extinct, O.adsorbed_none_extinct, O.adsorbed_all_extinct, O.total_A_extinct, O.total_NA_extinct, O.total_all_extinct, O.total_none_extinct];
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

    % Also print a csv of the concentration tracker
    csv_name = options.file_prefix + "curve" + options.file_suffix + "_" + O.sim_id + ".csv";
    writematrix(O.global_species_counts, csv_name);

    % Also save all the structs for later use if desired
    info_name = options.file_prefix + "info" + options.file_suffix + "_" + O.sim_id + ".mat";
    save(info_name, "p", "I", "O");

    fprintf("Data saved. Happy days! >w<\n");

end