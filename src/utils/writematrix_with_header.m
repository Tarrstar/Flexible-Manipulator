function writematrix_with_header(filename, data, header)

    fid = fopen(filename, 'w');

    if fid < 0
        error('Cannot open file: %s', filename);
    end

    fprintf(fid, '%s\n', header);
    fclose(fid);

    dlmwrite(filename, data, '-append');
end