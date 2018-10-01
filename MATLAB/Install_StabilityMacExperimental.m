function [] = Install_StabilityMacExperimental()
% Run this function from within the matlab folder to install stability
curr_folder = pwd;

cd ..;

disp('Compiling files...');
%setenv('PATH', [getenv('PATH') ':/usr/local/include/lemon:/usr/local/include']);
%mex -v GCC='/usr/local/bin/gcc-4.9' CXX='/usr/local/bin/g++-4.9' MATLAB/louvain_matlab_interface.cpp  -output stability_louvain CXXFLAGS="-std=gnu++0x -fPIC" -I./ -lemon

% Detailed compilation with expokit:
% !g++-4.6 -c matlab/louvain_matlab_interface.cpp -I./ -I/usr/local/MATLAB/R2010b/extern/include -I/usr/local/MATLAB/R2010b/simulink/include -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -pthread -std=gnu++0x  -DUSE_BOOST -DMX_COMPAT_32 -O3 -DNDEBUG -lemon -L/usr/local/include/lemon 
% %-Wall g++ is sufficient if version >=4.4
% 
% !gfortran -c -I./ -I/usr/local/MATLAB/R2010b/extern/include -I/usr/local/MATLAB/R2010b/simulink/include -fexceptions -fPIC -fno-omit-frame-pointer  -DUSE_BOOST -DMX_COMPAT_32 -O3  /home/mts09/repositories/group_repository/graph-codes/cliques/cxx/cliques/expokit.f 
% %-Wall
% 
% !g++-4.6 -O3 -pthread -shared -Wl,--version-script,/usr/local/MATLAB/R2010b/extern/lib/glnxa64/mexFunction.map -Wl,--no-undefined -o  stability_louvain.mexa64  louvain_matlab_interface.o expokit.o -Wl,-rpath-link,/usr/local/MATLAB/R2010b/bin/glnxa64 -L/usr/local/MATLAB/R2010b/bin/glnxa64 -lmx -lmex -lmat -lm -lgfortran -llapack -lemon -lblas 
% %-Wall 

% compilation mac
!/usr/local/bin/g++-5 -c -DMX_COMPAT_32   -DMATLAB_MEX_FILE  -I./  -I"/Applications/MATLAB_R2015a.app/extern/include" -I"/Applications/MATLAB_R2015a.app/simulink/include" -std=gnu++0x -fPIC -O2 -DNDEBUG matlab/louvain_matlab_interface.cpp -o louvain_matlab_interface.o
!/usr/local/bin/g++-5 -L/usr/local/include/lemon -Wl,-twolevel_namespace -undefined error -arch x86_64 -mmacosx-version-min=10.9 -Wl,-syslibroot,/usr/local/include -Wl,-syslibroot,/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk -framework Cocoa -bundle  -Wl,-exported_symbols_list,"/Applications/MATLAB_R2015a.app/extern/lib/maci64/mexFunction.map" -O -Wl,-exported_symbols_list,"/Applications/MATLAB_R2015a.app/extern/lib/maci64/mexFunction.map" louvain_matlab_interface.o   -lemon   -L"/Applications/MATLAB_R2015a.app/bin/maci64" -lmx -lmex -lmat -o stability_louvain.mexmaci64

disp('Moving build files to bin directory...');

movefile('*.mex*',curr_folder);
!rm louvain_matlab_interface.o
cd(curr_folder);

disp('Done');


