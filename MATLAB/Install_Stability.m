function [] = Install_Stability()
% Run this function from within the matlab folder to install stability
curr_folder = pwd;

cd ..;

disp('Compiling files...');
% setenv('PATH', [getenv('PATH') ':/usr/local/bin:/usr/local/lib'])
% !export PATH=$PATH:/usr/local/lib:/usr/local/bin:/usr/local/include
mex -v MATLAB/louvain_matlab_interface.cpp  -output louvain CXXFLAGS="-std=gnu++0x -fPIC" -I./CPP/cliques/ -I/usr/local/include -lemon
disp('Moving build files to bin directory...');

movefile('*.mex*',curr_folder);

% remove part below if you do not want the files ot be added to your path
cd(curr_folder);
disp('Adding bin directory to path...');
path(genpath(pwd),path);
savepath;

disp('Done');


