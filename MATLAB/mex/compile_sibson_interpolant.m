function compile_sibson_interpolant(compilationOptions)

% Validate input: Check if required fields exist in the struct
if (nargin < 1)
    compilationOptions = struct();
end

% Eigen directory (required), e.g. '/usr/local/include/eigen-3.4.0'
if ~isfield(compilationOptions, 'eigenDir') || isempty(compilationOptions.eigenDir)
    error('Please supply Eigen path in compilationOptions.eigenDir');
end

% LibIGL directory (required), e.g. '/usr/local/include/libigl/include'
if ~isfield(compilationOptions, 'iglDir') || isempty(compilationOptions.iglDir)
    error('Please supply LibIGL path in compilationOptions.iglDir');
end

% CGAL directory (required), e.g. '/usr/local/include/CGAL-6.0/include'
if ~isfield(compilationOptions, 'CGALDir') || isempty(compilationOptions.CGALDir)
    error('Please supply CGAL path in compilationOptions.CGALDir');
end

% Boost directory (optional, defaults to system locations)
if ~isfield(compilationOptions, 'boostDir') || isempty(compilationOptions.boostDir)
    if isfolder('/usr/include/boost')
        boostDir = '/usr/include/boost';
    elseif isfolder('/usr/local/include/boost')
        boostDir = '/usr/local/include/boost';
    else
        error('Boost library not found in default locations');
    end
else
    boostDir = compilationOptions.boostDir;
end

% TBB directory (optional, defaults to system locations)
if ~isfield(compilationOptions, 'tbbDir') || isempty(compilationOptions.tbbDir)
    if isfolder('/usr/include/tbb')
        tbbDir = '/usr/include/tbb';
    elseif isfolder('/usr/local/include/tbb')
        tbbDir = '/usr/local/include/tbb';
    else
        error('TBB library not found in default locations');
    end
else
    tbbDir = compilationOptions.tbbDir;
end

% Verbose flag (optional, default to false)
if ~isfield(compilationOptions, 'verbose') || isempty(compilationOptions.verbose)
    verbose = false;  % Default to no verbose output
else
    verbose = compilationOptions.verbose;
end

[projectDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(projectDir);

% Compilation options
CXXOPTIMFLAGS = '"-O3"';
CXXFLAGS = '"$CXXFLAGS -march=native -fopenmp -fPIC -std=c++17"';
LDFLAGS = '"$LDFLAGS -fopenmp"';
CGALCompileFlags = '-DCGAL_EIGEN3_ENABLED=true ';

% Include flags
includeFlags = [ ...
    '-I' compilationOptions.eigenDir ' ' ...
    '-I' boostDir ' ' ...
    '-I' compilationOptions.CGALDir ' ' ...
    '-I' compilationOptions.iglDir ' ' ...
    '-I' tbbDir ' ' ...
    '-I/usr/include:/usr/local/include ' ];

% Library flags
libFlags = '-L/usr/lib:/usr/lib/x86_64-linux-gnu:/usr/local/lib ';
libFlags = [libFlags ...
    '-lgmp -lmpfr -lboost_thread -lboost_system -ltbb '];

% MEX compilation string
if verbose

    mexString = [ 'mex -v sibsonInterpolant.cpp ' ...
        'CXXOPTIMFLAGS=' CXXOPTIMFLAGS ' ' ...
        'CXXFLAGS=' CXXFLAGS ' ' ...
        'LDFLAGS=' LDFLAGS ' ' ...
        CGALCompileFlags includeFlags libFlags ];

else

    mexString = [ 'mex sibsonInterpolant.cpp ' ...
        'CXXOPTIMFLAGS=' CXXOPTIMFLAGS ' ' ...
        'CXXFLAGS=' CXXFLAGS ' ' ...
        'LDFLAGS=' LDFLAGS ' ' ...
        CGALCompileFlags includeFlags libFlags ];

end

% Execute the MEX command
eval(mexString);

end
