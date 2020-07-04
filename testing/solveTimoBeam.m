%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Nguyen Thanh Vien Minh
% Email: minh.nguyen@ikm.uni-hannover.de
% Working at: Institut fuer Kontinuumsmechanik, Hannover, Germany
% Website: https://www.ikm.uni-hannover.de/kontinuumsmechanik.html?&no_cache=1&L=1
% ------------------------------------------------------------------------------------------
% If you have any question, please do not hesitate to contact me
% immediately via my email.
%
% Please cite the paper if you would like to use my source code as a part of your
% project
% "A Virtual Element Method for 2D linear elastic fracture analysis" - V.M
% Nguyen-Thanh; X. Zhuang; H. Nguyen-Xuan; T. Rabczuk; P. Wriggers
%
% Thank you and have fun with my code, enjoy it !!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This programme aims to test Virtual Element Method programme by trying to solve the 
% TIMOSHENKO's BEAM PROBLEM
% Author: Minh T.V Nguyen
% Mesh type: Polygonal elements
%                                
% \----------------|--------------| |
% \|               |              |  |
% \|               |              |   |
% \|               |              |   |
% \|               |              |  |
% \-------------------------------| |   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------------------------------------------------------
%                   CLEAR ENVIRONMENT
%-------------------------------------------------------------------------------------------------------------------------
clear all;
close all;
projectdir = '/home/minh/PHD/dev/VEM_lightweight/';
addpath(projectdir)

%-------------------------------------------------------------------------------------------------------------------------
%                   LOAD POLYGONAL MESH
%-------------------------------------------------------------------------------------------------------------------------
% mesh = load('/home/minh/dev/VEMMinhIKM/meshes/TimoBeam_Poly');
% mesh = load('/home/minh/dev/VEMMinhIKM/meshes/TimoBeam_voronoi');
mesh = load('/home/minh/PHD/dev/VEM_lightweight/meshes/TimoBeam_voronoi');
% mesh = load('/home/minh/dev/VEMMinhIKM/meshes/TimoBeam_Q4');
vertices = mesh.vertices;
elements = mesh.elements;
%-------------------------------------------------------------------------------------------------------------------------
%                   SETTING BOUNDARY CONDITIONS
%-------------------------------------------------------------------------------------------------------------------------
bc_idx = 1;
loading_idx = 1;
loading_idx2 = 1;
loading_idx3 = 1;

for i = 1 : size(vertices,1)
    if (vertices(i,1) == 0.0)
        boundary(bc_idx,1) = i;
        bc_idx = bc_idx + 1;            
    elseif (vertices(i,1) == 16.0)
        loading(loading_idx,1) = i;
        loading_idx = loading_idx + 1;
    end
    if (vertices(i,2) == -2.0)
        loading2(loading_idx2,1) = i;
        loading_idx2 = loading_idx2 + 1;
    elseif (vertices(i,2) == 2.0) 
        loading3(loading_idx3,1) = i;
        loading_idx3 = loading_idx3 + 1;
    end
end

 dirichlet{1} = boundary;
 neumann{1} = loading;
 neumann{2} = loading2;
 neumann{3} = loading3;

for i = 1 : length(dirichlet)
    [Ab, Ib] = sort(vertices(dirichlet{i}(:,1),2),'descend');  %Ib contains the corresponding indices of Ab.
    boundary_temp = dirichlet{i};
    boundary_temp(:,1) = dirichlet{i}(Ib,1);
    dirichlet{i} = boundary_temp(:,:);
end

% Sorting along y-direction for traction
for i = 1 : length(neumann)
    if i == 1
        [At, It] = sort(vertices(neumann{i}(:,1),2));  %Ib contains the corresponding indices of Ab.
    else
        [At, It] = sort(vertices(neumann{i}(:,1),1));  %Ib contains the corresponding indices of Ab.
    end 
    loading_temp = neumann{i};
    loading_temp(:,1) = neumann{i}(It,1);
    neumann{i} = loading_temp(:,:);
end
dirichlet{1} = MeshOperation.convertNodes2Edges(dirichlet{1});
neumann{1}   = MeshOperation.convertNodes2Edges(neumann{1});
neumann{2}   = MeshOperation.convertNodes2Edges(neumann{2});
neumann{3}   = MeshOperation.convertNodes2Edges(neumann{3});

mesh.elements  = elements;
mesh.vertices  = vertices;
mesh.dirichlet = dirichlet;
mesh.neumann   = neumann;
mesh.elmType_  = 'Poly';
%-------------------------------------------------------------------------------------------------------------------------
%                   1. SETTING MATERIAL PROPERTIES
%                   2. SETTING VEM METHOD to SOLVE
%-------------------------------------------------------------------------------------------------------------------------
P = 1000;
D = 4;
I = D^3/12;
E = 1000000;
L = 16;
nu = 0.3;
k = 1; % the order of approximation scheme.
dim = 2; % Dimensional problem
thickness = 1;
typeF = 'traction';

% u Dirichlet
uD = @(x,y)[P*y/(6*E*I) .* ((6*L - 3*x).*x + (2+nu)*(y.*y - D^2/4)), ...
           -P/(6*E*I) * (3*nu*(y.*y).*(L-x) + (4+5*nu)*(D^2)*x/4 + (3*L-x).*(x.*x))];   
         
% t Neumann
tN = @(x,y)[y-y,-P/(2*I)*(D^2/4 - y.*y)];

uFunction{1} = uD;
tFunction{1} = tN;
tFunction{2} = @(x,y)[x-x,y-y];
tFunction{3} = @(x,y)[x-x,y-y];
functions.uFunction = uFunction;
functions.tFunction = tFunction;

global mat;
mat = Material(E, nu, 'plane_stress', thickness);
global method;
method = Method('vem', dim, k, typeF);

%-------------------------------------------------------------------------------------------------------------------------
%               VEM SOLVER
%-------------------------------------------------------------------------------------------------------------------------
VEM = VirtualElementMethod(mat, method, mesh, mesh.dirichlet, mesh.neumann, functions.uFunction, functions.tFunction);
[U, Pi_Star] = VEM.solveLinearElastic2D();
[strain_elms, stress_elms] = VEM.executePostProcessing(Pi_Star,U);
VEM.delete();
%-------------------------------------------------------------------------------------------------------------------------
%               PLOTTING
%-------------------------------------------------------------------------------------------------------------------------
Draw = Plotting(mesh, dim);
Draw.plotDeformation(U);
Draw.plotStress(stress_elms, 3);

