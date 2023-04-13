%% This code will take a base feb file (with geometry and boundary conditions
%and partition the elements into individual materials (neo-Hookean). This code
%will also run febioStudio from command line via Matlab. Check the command line
%code to confirm FebioStudio will run on your computer
%Inputs:
% -base feb file with geometry and BC
% -material parameters for each element in FEBio mesh
%
%Output
% -new feb file with partitioned elements
% -xplt and log file
% -text file containing F for each element for each step
% -text file containing stress for each element for each step
%
%written by: Liz Gacek
%Date: 1/27/2022

%Input Structure variables for write_FEBio_HGO
Input.base_filename = 'cyl.feb';       %name of base feb file
Input.filename = 'cyl_multinet.feb';  %name of new feb file
Input.b = '/';                          %dependant on operation system
Input.noElem = 320;                    %number of elements to partition
Input.outFtensor = 'cyl_multinet_F.txt';           %name of F text file output
Input.outStress = 'cyl_multinet_stress.txt';       %name of stress text file output
Input.P = 5000;

%Material parameters for each element
Input.matrixP = 0.499;                  %neo-Hookean poisson's ratio
C1 = 0.01e6;                 %material properties for each element
C_input = C1.*ones(Input.noElem);
Input.Matfits.C = C_input;

dir_id = 'newNetworks';
mkdir(dir_id);

%write new Feb file
write_FEBio2(Input)

% network properties
% network scale
scale= 10e-6;
% fiber radii
radsE = 20e-9;
radsC = 150e-9;
radsA = 10e-9;

num_nA = 30; % 3D nodes
num_nE = 20; % planar nodes

plot_net=0;

% for net_id = 1:Input.noElem
%     % produce networks
%     createNets(num_nE, num_nA, radsE, radsC, radsA, scale, net_id, dir_id, plot_net)
% end

%run FebioStudio from command line
% executeCL(Input.filename);
%
% function executeCL(file) % execute commandline
%     run_string=['febio3 -i ',file]; %,' -silent'
%     system(run_string); %Run FEBio job
% end