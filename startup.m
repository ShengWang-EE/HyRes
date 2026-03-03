function startup()
    root = fileparts(mfilename('fullpath'));

    % 1) Project code
    addpath(genpath(root), '-begin');

    % 2) External toolboxes (relative to project root)
    toolbox_root = fullfile(root, '..', 'Matlab repository', 'toolbox_user');
    addpath(genpath(fullfile(toolbox_root, 'YALMIP-master')), '-begin');
    addpath(genpath(fullfile(toolbox_root, 'matpower8.1')), '-begin');

    % 3) Gurobi (C:\gurobi1102)
    gurobi_home = 'C:\gurobi1102\win64';
    gurobi_matlab = fullfile(gurobi_home, 'matlab');
    gurobi_bin = fullfile(gurobi_home, 'bin');

    if isfolder(gurobi_matlab)
        addpath(gurobi_matlab, '-begin');
        setenv('GUROBI_HOME', gurobi_home);
        cur_path = getenv('PATH');
        if ~contains(lower(cur_path), lower(gurobi_bin))
            setenv('PATH', [gurobi_bin ';' cur_path]);
        end
    else
        warning('Gurobi MATLAB path not found: %s', gurobi_matlab);
    end

    rehash toolboxcache
end
