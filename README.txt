Non-Affine Network Solver (NaNS) Plugin for FEBio

Ryan R Mahutga
Department of Biomedical Engineering
University of Minnesota, MN, USA

Updated: 04/11/2023

License: GNU General Public License Version 3

Disclaimer: 
===========
This work is an independent publication and is 
neither affiliated with, nor authorized, sponsored, or approved
by Microsoft Corporation, Eigen, or FEBio. 

This plugin is free software, which can be redistributed and 
modified under the terms of the GNU General Public License 
Version 3 as published by the Free Software Foundation. This plugin
is distributed in the hope that it will be useful and is provided 
without warranty of any kind. See the GNU General Public License 
for more details.

Interfacing the NaNS Plugin with FEBio: 
=======================================
In order to make this plugin work with FEBio, you need to include
 the .dll file (or .so file) in the location of the executable 
 FEBio4.exe. This is typically somewhere like 
 C:\Program Files\FEBioStudio2\bin. You also need to modify the 
 file febio.xml also in the same location by adding 
<import>NONAFFINE_NETWORK_SOLVER_PLUGIN.dll</import> between 
<febio_config version="3.0"> and </febio_config>. 

Setting up the NaNS Material in a simulation:
=============================================
The easiest way to set up a network material to be used in a 
simulation is to create the model, assign neo-Hookean materials to 
your parts that you want to be networks, export the FE model, then 
modify the corresponding .feb file using a text editor (I like
Notepad++ for this). You can update the material you want to be a
network by modifying the neo-Hookean material to be:

		<material id="1" name="M1" type="PeriodicNetwork">
			<netNum>1</netNum>
			<netSave lc="1">1</netSave>
			<netSolve>1</netSolve>
			<E>50000</E>
			<v>4.990000e-01</v>
		</material>

In this material definition, netNum is the numerical identifier of the 
network to be used (so PeriodicNetwork1.txt), netSave when it is 1 will
save pertinent network data to .txt files, netSolve when it is 1 solve 
the network using Newton Raphson iteration, when it is 0 if uses an affine
approximation, E is the neo-Hookean ground matrix Young's Modulus and v 
is the material Poisson's ratio. 

Structure of the Network File: 
===============================
The network file is structure like 

 719  90  1.200000e-05  85 -2 -2
 14  38  0  1  0  3  4.744273e-01  1.000000e-07
 14  17  0  1  0  3  5.181404e-01  1.000000e-07
 .
 .
 .
 2.522269e-01  1.298863e-01  1.299305e-01
-4.053110e-01  2.481414e-01  1.917741e-01
 3.624445e-01 -1.188753e-02  3.109951e-01
-1.460525e-01  7.657790e-03  2.038366e-01

where the first line gives the number of fibers, the number of nodes,
the scale of the network (i.e. m per computational unit length), and
the next three numbers are node numbers for stabilization. In this 
case node 85 if fixed in x,y,z to keep the periodic network from sliding
around in it's domain and the next two numbers (-2) are just placeholders.
The next, in this case 719, lines are fiber definitions. The first two 
values indicate which two nodes are connected. The next three values 
indicate how many times and in which direction the fiber crosses a boundary.

As an example, the first 5 values of 14  38  0  1  0 indicate that node 14 
is connected to 38 and to get from node 14 to node 38 we cross the x boundary
0 times, the y boundary 1 time, and the z boundary 0 times. Thus, the vector 
from node 14 to node 38 would be the difference between node 38 shifted one 
RVE size (basically moved over 1 RVE) in the y-direction and node 14 in the
original RVE. 

The 6th value is the fiber type which changes the constitutive behavior (in 
this work 1 is helical (collagen), 2 is linear (elastin), and 3 is active 
(actin). The 7th value is the fiber rest length in computational length units.
The final value is the fiber radius in real units [m].  
 
After the (719) fibers are defined, the file contains the location of the (90) 
nodes in computational space as x,y,z values. 

Generating Periodic Networks: 
=============================
The simplest way to generate periodic networks is using a Delaunay triangulation
algorith on nodes repeated in a 3x3x3 grid and extracting the center 1x1x1 RVE 
while tracking fiber crossings. Included with this plugin are several Matlab 
routines which create various types of networks and perform various manipulations
to them (i.e. pare them to lower connectivity).  

Generation of Advanced FE Simulations (i.e. every element gets its own network):
=================================================================================
Also included with this plugin is a Matlab code designed to take a simulation with 
one part and split it into one part for each element so one can easily use one
network per element. 

Working with the Source code:
=============================
To modify and recompile the plugin, the simplest avenue is to use VisualStudio2019 
or newer. This should open the .sln file with all the necessary compilation
information in tact. If not, the code should be compiled as a .dll (dynamic link library)
using Windows SDK 10.0 with platform toolset VisualStudio2019 v142. 

To compile you will need to update the Include directories and the Library directories. 
First ensure you have the FEBio SDK which is usually included as part of the FEBio download,
although you may have to redowload FEBio and select it as an optional part of the install. 
In include directories you should have: C:\Program Files\FEBioStudio2\sdk\include  
Under libraries you should have C:\Program Files\FEBioStudio2\sdk\lib\Release

The code also relies on the Eigen library available here https://eigen.tuxfamily.org/, which 
is already included with the Source code provided in the FEBIO_PLUGIN subfolder.  

Examples Files: 
===============
Included as a part of this plugin are many example files which were used in the 
creation of the paper on this Plugin (update citation on publication). 
