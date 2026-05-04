function [O] = csvmaker(O)
    % Makes a more detailed csv output with aggregate data AND
    % concentration tracking throughout each individual pixel

    arguments (Input)
        O;
    end

    % Make labels for aggregate data
    csv_data = O.global_species_counts;
    labels = ["pixel"; zeros(length(O.time), 1)];
    csv_data = [csv_data, labels];

    % Get data for every individual pixel
    for pixel = 1:length(O.concentration_tracker)
        current_data = [O.time, O.concentration_tracker{pixel}];
        current_data = [current_data, ones(length(O.time), 1)*pixel];
        csv_data = [csv_data; current_data];
    end

    O.csv_data = csv_data;

end