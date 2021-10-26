function varargout = inject(what, filename)
    directory = fullfile('./data', what);

    if ~isfolder(directory)
        error('The directory does not exist: %s', directory);
    end

    filepath = fullfile(directory, filename);

    fid = fopen(filepath, 'r');

    if fid > 0
        header = strsplit(fgetl(fid), ' ');

        [lines, columns] = header{:};
        lines = str2double(lines);
        columns = str2double(columns);

        data = textscan(fid, '%f');

        data = data{:};

        if lines == size(data, 2) && columns == size(data, 1)
            data = transpose(data);
        end

        fclose(fid);
    else
        error('Could not open the file: %s', filepath);
    end

    if strcmp(what, 'window_alignment_stage')
        varargout = {data(1), data(2), data(3)};
    elseif strcmp(what, 'hdl/window_alignment_stage')
        varargout = {data(1), data(2), data(3)};
    else
        varargout = {data};
    end
end
