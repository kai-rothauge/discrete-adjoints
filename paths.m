function paths(set)

if strcmpi(set,'set')
    addpath(genpath([pwd '\bin']));
    addpath(genpath([pwd '\settings']));
%     % Add necessary paths for current test
%     addpath([pwd '\Time-Stepping']);
%     split = strsplit(time_stepping,' ');
%     addpath(genpath([pwd '\Time-Stepping\' split{1} '\' split{2}]));
%     addpath(genpath([pwd '\Misfit Functions\' misfit_function]));
%     addpath(genpath([pwd '\PDE Examples\' PDE_example]));
% 
%     cd([pwd '\Derivatives\Testing'])
else
    rmpath(genpath([pwd '\bin']));
    rmpath(genpath([pwd '\settings']));
%     % Remove paths from previous test
%     warning('off','MATLAB:rmpath:DirNotFound')
%     rmpath(genpath([pwd '\Time-Stepping']));
%     rmpath(genpath([pwd '\Misfit Functions']));
%     rmpath(genpath([pwd '\PDE Examples']));
%     warning('on','MATLAB:rmpath:DirNotFound')
end
