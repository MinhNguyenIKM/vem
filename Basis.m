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

function [M, GradM] = Basis(k, x, y, xCentroid, yCentroid, diameter)
%BASIS Summary of this function goes here
%   Detailed explanation goes here
% Input:
%         k: approximation order
%         x: dof at i-th (x direction) vertice
%         y: dof at i-th (y direction) vertice
    if k == 1
%         % BASIS OF G. H. PAULINO 1
%         M       =     {[1; 0];          [0; 1];         [y-yCentroid; -(x-xCentroid)];        [x-xCentroid; 0];         [0; y-yCentroid];         [y-yCentroid; x-xCentroid]};
%         GradM   =     {[0 0; 0 0],      [0 0; 0 0],     [0 1; -1 0],    [1 0; 0 0],     [0 0; 0 1],     [0 1; 1 0]};
%         
%         % BASIS OF G. H. PAULINO 2
%         M       =     {[1; 0];          [0; 1];         [y; -(x)];        [x; 0];         [0; y];         [y; x]};
%         GradM   =     {[0 0; 0 0],      [0 0; 0 0],     [0 1; -1 0],    [1 0; 0 0],     [0 0; 0 1],     [0 1; 1 0]};
% 
%         % BASIS OF WILHELM RUST
%         M       =     {[1; 0]; [0; 1]; [y; -(x)];  [y; x]; [x; -(y)]; [x; y]};
%         GradM   =     {[0 0; 0 0],      [0 0; 0 0],     [0 1; -1 0],    [0 1; 1 0],     [1 0; 0 -1],     [1 0; 0 1]};

        % MY BASIS 1
        M       =     {[1; 0]; [0; 1]; [y-yCentroid; -(x-xCentroid)];  [y-yCentroid; x-xCentroid]; [x-xCentroid; -(y-yCentroid)]; [x-xCentroid; y-yCentroid]};
        GradM   =     {[0 0; 0 0],      [0 0; 0 0],     [0 1; -1 0],    [0 1; 1 0],     [1 0; 0 -1],     [1 0; 0 1]};
        
        % MY BASIS 2
%         M       =     {[1; 0]; [0; 1]; [(y-yCentroid)/diameter; -(x-xCentroid)/diameter];  [(y-yCentroid)/diameter; (x-xCentroid)/diameter]; [(x-xCentroid)/diameter; -(y-yCentroid)/diameter]; [(x-xCentroid)/diameter; (y-yCentroid)/diameter]};
%         GradM   =     {[0 0; 0 0],      [0 0; 0 0],     [0 1; -1 0]/diameter,    [0 1; 1 0]/diameter,     [1 0; 0 -1]/diameter,     [1 0; 0 1]/diameter};
%           

    else 
       disp('We are implementing case k > 1');       
    end
end

