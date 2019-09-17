###############################################################################  
Copyright (C) 2018 A. Delmotte, M. Schaub, S. Yaliraki, M. Barahona

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

###############################################################################

-----------------------------------------------------------------------------
Generalized Louvain optimization (for graph partitioning problems)
-----------------------------------------------------------------------------

The code implements a generalized Louvain optimization algorithm which can be used to
optimize several objective functions, e.g., the ones discussed in the article:

Michael T. Schaub, Jean-Charles Delvenne, Renaud Lambiotte, Mauricio Barahona
"Multiscale dynamical embeddings of complex networks"
https://arxiv.org/abs/1804.03733

This code emerged from a previous repository that implemented the Louvain algorithm
for optimzation of Markov stability, see here
https://github.com/michaelschaub/PartitionStability
A legacy version of this code -- including the old C++ backend (no lemon library), with
an improved Matlab interface is included within this repository for convenience.
Please see the README file within the respective folder for further details.

This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie grant agreement No 702410.


***If you make use of any part of this toolbox, please cite our work.***

The C++ optimization toolbox (cliques) can be used independently or be called from Matlab.
If you want to use the code independently, you may also want to make use of the FORTRAN 
code implementing the computation of the matrix exponential function (see FORTRAN folder).

For detailed instructions on how to compile the code in MATLAB see below.
If you find a bug or have further comments, please send an email and if 
necessary the input file and the parameters that caused the error.

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Authors   : M. Schaub  
Email     : mschaub[at]mit.edu
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

###############################################################################

-----------------------------------------------------------------------------
Contributions to the code
-----------------------------------------------------------------------------

Please see CODE_HISTORY.txt for more information.

###############################################################################

-----------------------------------------------------------------------------
How to install the stability package
-----------------------------------------------------------------------------

Prerequisites:
a) Install Lemon Graph library -- a version is provided in the folder CPP/lemon-lib 
    for convenience. See https://lemon.cs.elte.hu/trac/lemon for further details

1. Open Matlab

2. Make sure you have a C++ compiler installed
  * For Linux, you can find one here: 
    http://www.gnu.org/software/gcc/
  * For Windows, you can use Visual C++ express: 
    http://www.microsoft.com/express/Windows/

3. Make sure mex is properly configured in Matlab:
  * Type "mex -setup" in Matlab, and choose your compiler.

4. In Matlab, go into the directory of the Stability toolbox.

5. Type "Install_Stability" in the Matlab command window.
  * If you get an error message concerning the libstdc++.so file, 
    you may want to try the following manipulation:

        cd "Matlab_root_directory"/sys/os/glnx86/
        sudo mv libgcc_s.so.1 libgcc_s.so.1.back
        sudo mv libstdc++.so.6 libstdc++.so.6.back

6. You will get a messge asking whether the stability toolbox should 
   be added to your Matlab path. Answering yes will allow you to use 
   the stability toolbox functions as standard Matlab functions.
            
7. Type "help stability" in Matlab to discover how to use the code.

8. Try this example to check that everything is working:
    
        cd('MATLAB/demo');   % go into the demo directory 
        load demo;    % load data and then run stability
        [S, N, VI, C] = partition_stability(Graph,Time,'plot','v');
        % for a more advanced example see also the example analysis 
        % of the Harari highland data in the demo folder

NOTES:

* The install script provides the option to add the bin folder to your 
Matlab path. This will enable you to use stability as a standard Matlab 
function from any directory. If you don't want this option any more,
just remove it from the path by going in File/Set Path.

* If you get a warning message concerning savepath, and you want the 
stability code to be in your path, go, after the installation, in 
File/Set Path, and choose "save". Then choose where you want pathdef.m
to be saved. If at the next matlab startup, you notice that stability is
not in your matlab path anymore, try editing/creating the "startup.m" file
from your matlab user folder (type userpath to know where it is located)
and add the following line: addpath(' path to bin folder of stability 
package '). Alternatively, if you are the only user on your machine, you
can start matlab as a superuser ("sudo matlab" in linux) and rerun the
"Install_Stability" script. This will permanently add the stability folder 
in the path for all users.

* To speed up the calculations, you might consider adding the
option 'noVI'. This disables the calculation of the variation of information, 
which is usually slow at small Markov times, when the number of 
communities found is big. 
Another option is to decrease the number of optimisations on which the variation 
of information is calculated. To do so, add the option 'M' and put a value
such that M < L (L is the number of louvain optimisations).  
Example:
 
    [S, N, VI, C] = partition_stability(Graph,time,'plot','v', 'L', 100, 'M', 10);

