function eject(what, data, filename)
    directory = fullfile('./data', what);

    if ~isfolder(directory)
        mkdir(directory)
    end

    filepath = fullfile(directory, filename);

    [fid, msg] = fopen(filepath, 'w+');

    if strcmp(what, 'receiver')
        spec = "%.16f%+.16fi";
    else
        spec = "%d";
    end

    % spec = strjoin(repelem(datatype, size(data, 2)));
    spec = strcat(spec, "\n");

    if fid > 0
        % Print the size of the data
        fprintf(fid, '%d %d\n', size(data));

        % Print the actual data
        if isreal(data)
            fprintf(fid, spec, data');
        else
            fprintf(fid, spec, [real(data(:)), imag(data(:))].');
        end

        fclose(fid);
    else
        error('Could not open the file: %s (%s)', filepath, msg);
    end
end
