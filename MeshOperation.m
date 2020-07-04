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

classdef MeshOperation < handle
   methods(Static)
       function adjElm = get_AdjElem(nodeIdx, mesh)
            elements = mesh.elements;
            adjElm = zeros(length(elements),1);
            i = 0;
            for e = 1 : length(elements)
               connectivity = elements{e};
%                [yes, ~] = ismember(nodeIdx, connectivity);
               yes = any(connectivity == nodeIdx);
               if (yes)
                  i = i + 1;
                  adjElm(i) = e;
               end
            end
            adjElm = adjElm(1:i); % Trim extra zeros from the results
       end
       
       function adjElm = getAdjElem(nodeIdx, grid)
            elmMatrix = grid.elements.elmMatrix; % elements/connectivity is represented in matrix form
            indices = find(elmMatrix == nodeIdx);
            adjElm = rem(indices-1,size(elmMatrix,1))+1; % get the subscripts (m,n) of elmMatrix. Here m-th is the element index.
       end
     
       function bool = is_CornerNode(nodeIdx, grid)
           bool = ismember(nodeIdx,grid.nodes.cornerNodes);
       end
       
       function [nodeArray, num] = get_RecoveredNodes(adj_elm, mesh)
           nodeArray = [];
           elements = mesh.elements;
           for i = 1 : length(adj_elm)
                verts = elements{adj_elm(i)};
                for j = 1 : length(verts)
                    if ~ismember(verts(j), nodeArray)
                        nodeArray(end+1) = verts(j);
                    end
                end
           end
           num = length(nodeArray);
       end
       
       function elmMat = createElmMatrix(elmCell)
           [max_size, ~] = max(cellfun('size', elmCell, 1));
           elmMat = zeros(length(elmCell),max_size);
           for e = 1 : length(elmCell)
               elmMat(e,1:length(elmCell{e})) = elmCell{e};
           end
       end
       % Grid structure
       % Grid |------ nodes |--- num
       %      |             |--- coordinates
       %      |             |--- cornerNodes
       %      |
       %      |--- elements |--- num
       %                    |--- connectivity
       %                    |--- areas
       %                    |--- centroids
       %                    |--- diameters
       %                    |--- boundaryElms
       %                    |--- cornerElms
       %                    |--- elmMatrix
       %
       function grid = get_MeshInfo(mesh)
            nodes = mesh.vertices;
            elements = mesh.elements;
            grid.nodes.num = size(nodes,1);
            grid.nodes.coordinates = nodes(:,:);
            grid.nodes.cornerNodes = [];
            grid.elements.elmMatrix = MeshOperation.createElmMatrix(elements);
            grid.elements.num = length(elements);
            grid.elements.connectivity = elements;
            grid.elements.boundaryElms = [];
            grid.elements.cornerElms = [];
            grid.elements.centroids = [];
            grid.elements.diameters = [];
            grid.elements.areas = [];
            edges = [];
            for i = 1 : grid.elements.num
                con = grid.elements.connectivity{i};
                verts = nodes(con,:);
                [area, centroid, diameter] = MeshOperation.get_ElementInfo(verts);
                grid.elements.areas(i) = area;
                grid.elements.centroids(i,[1 2]) = centroid;
                grid.elements.diameters(i) = diameter;
                for j = 1 : length(con)
                   c = j;
                   n = j + 1;
                   if j == length(con)
                       n = 1;
                   end
                   larger  = max(con(c),con(n));
                   smaller = min(con(c),con(n));
                   edges(end+1,[1,2]) = [smaller,larger]; 
                end
            end
            grid.surface.edges = edges;
            % ------------------------------
            % get surface edges algorithm
            % ------------------------------
            [b, i1] = unique(grid.surface.edges, 'rows', 'first');
            [b, i2] = unique(grid.surface.edges, 'rows', 'last');
            grid.boundary.edges = b(i1==i2,:);
            %-------------------
            for e = 1 : grid.elements.num
               con = grid.elements.connectivity{e};
               count = 0;
               marker = [];
               for j = 1 : length(con)
                   c = j;
                   n = j + 1;
                   if j == length(con)
                       n = 1;
                   end
                   larger  = max(con(c),con(n));
                   smaller = min(con(c),con(n));
                   if ismember([smaller,larger],grid.boundary.edges,'rows')
                       count = count + 1;
                       marker(end+1) = smaller; % marker to get the corner node
                       marker(end+1) = larger;
                   end                   
               end
               if count == 1
                   grid.elements.boundaryElms(end+1) = e;
               end
               if count > 1
                   grid.elements.cornerElms(end+1) = e;
                   [m,i1] = unique(marker,'first');
                   [m,i2] = unique(marker,'last');
                   grid.nodes.cornerNodes(end+1:end+length(m(i1~=i2))) = m(i1~=i2);
               end
            end               
       end
       
		function [area, centroid, diameter] =  getElementInfo(verts)
			[area, centroid, diameter] = MeshOperation.get_ElementInfo(verts);
		end

        function [area, centroid, diameter] = get_ElementInfo(verts)
            Ne = size(verts,1); 
            area_components = (verts(:,1) .* verts([2:end,1],2)) - (verts([2:end,1],1) .* verts(:,2));
            area = 0.5 * abs(sum(area_components));
            %centroid = sum((verts + verts([2:end,1],:)) .* repmat(area_components,1,2)) / (6*area);
            centroid = MeshOperation.polygonCentroid(verts);
            diameter = 0;
            for i = 1 : (Ne-1)
               for j = (i+1) : Ne
                   diameter = max(diameter, norm(verts(i,:) - verts(j,:)));
               end
            end
			
        end
        

		function diameter = polygon_diameter ( n, v )

			%*****************************************************************************80
			%
			%% POLYGON_DIAMETER computes the diameter of a polygon.
			%
			%  Discussion:
			%
			%    The diameter of a polygon is the maximum distance between any
			%    two points on the polygon.  It is guaranteed that this maximum
			%    distance occurs between two vertices of the polygon.  It is
			%    sufficient to check the distance between all pairs of vertices.
			%    This is an N^2 algorithm.  There is an algorithm by Shamos which
			%    can compute this quantity in order N time instead.
			%
			%  Licensing:
			%
			%    This code is distributed under the GNU LGPL license.
			%
			%  Modified:
			%
			%    02 February 2005
			%
			%  Author:
			%
			%    John Burkardt
			%
			%  Parameters:
			%
			%    Input, integer N, the number of vertices of the polygon.
			%
			%    Input, real V(2,N), the vertices.
			%
			%    Output, real DIAMETER, the diameter of the polygon.
			%
			diameter = 0.0;
			for i = 1 : n
    			for j = i + 1 : n
      				diameter = max ( diameter, ...
	        		sqrt ( ( v(1,i) - v(1,j) ).^2 + ( v(2,i) - v(2,j) ).^2 ) );
		    	end

			end
  			return
		end

        %
        % Initial date: 03/08/2017
        % 
        function [centroid] = polygonCentroid(varargin)
              %POLYGONCENTROID Compute the centroid (center of mass) of a polygon
              %
              %   CENTROID = polygonCentroid(POLY)
              %   CENTROID = polygonCentroid(PTX, PTY)
              %   Computes center of mass of a polygon defined by POLY. POLY is a N-by-2
              %   array of double containing coordinates of vertices.
              %
              %   [CENTROID AREA] = polygonCentroid(POLY)
              %   Also returns the (signed) area of the polygon. 
              %
              %   Example
              %     % Draws the centroid of a paper hen
              %     x = [0 10 20  0 -10 -20 -10 -10  0];
              %     y = [0  0 10 10  20  10  10  0 -10];
              %     poly = [x' y'];
              %     centro = polygonCentroid(poly);
              %     drawPolygon(poly);
              %     hold on; axis equal;
              %     drawPoint(centro, 'bo');
              % 
              %   References
              %   algo adapted from P. Bourke web page
              %
              %   See also:
              %   polygons2d, polygonArea, drawPolygon
              %
              %   ---------
              %   author : David Legland 
              %   INRA - TPV URPOI - BIA IMASTE
              %   created the 05/05/2004.
              %

              % Algorithme P. Bourke, vectorized version

              % HISTORY
              % 2012.02.24 vectorize code


              % parse input arguments
              if nargin == 1
                  var = varargin{1};
                  px = var(:,1);
                  py = var(:,2);
              elseif nargin == 2
                  px = varargin{1};
                  py = varargin{2};
              end

              % vertex indices
              N = length(px);
              iNext = [2:N 1];

              % compute cross products
              common = px .* py(iNext) - px(iNext) .* py;
              sx = sum((px + px(iNext)) .* common);
              sy = sum((py + py(iNext)) .* common);

              % area and centroid
              area = sum(common) / 2;
              centroid = [sx sy] / 6 / area;
        end
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % minh.nguyen@ikm.uni-hannover.de
        % Initialize function convertEdges2Nodes
        % Date: 13/04/2017
        % Convert Edges array or cell to Nodes array or cell respectively
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function nodes = convertEdges2Nodes(edges)
            if iscell(edges)
                for k = 1 : length(edges)
                   edgesArray = edges{k};
                   [dummy, idx] = unique(edgesArray, 'first');
                   nodes{k} = edgesArray(sort(idx));
                end
            else
                [dummy, idx] = unique(edges, 'first');
                nodes = edges(sort(idx));
            end
        end
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % minh.nguyen@ikm.uni-hannover.de
        % Initialize function convertNodes2Edges
        % Date: 13/04/2017
        % Convert Nodes array or cell to Edges array or cell respectively
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
        function edges = convertNodes2Edges(nodes)
            if iscell(nodes)
                for k = 1 : length(nodes)
                    nodesArray = nodes{k};
                    idx = 1;
                    for i = 1 : length(nodesArray)-1
                        edges{k}(idx,1:2) = [nodesArray(i), nodesArray(i+1)];
                        idx = idx + 1;
                    end
%                     edges{k}(idx,1:2) = [nodesArray(i+1), nodesArray(1)];
                end
            else
                idx = 1;
                for i = 1 : length(nodes)-1
                    edges(idx,1:2) = [nodes(i), nodes(i+1)];
                    idx = idx + 1;
                end
%                 edges(idx,1:2) = [nodes(i+1), nodes(1)];
            end
        end
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
        % Minh T.V Nguyen minh.nguyen@ikm.uni-hannover.de
        % Renumber label or connectivity in a best order for newest vertex
        % bisection adaptivity
        % Initial date: 02/05/2017
        % Original author: Long Chen - Short Bisection Implementation in
        % Matlab
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
        function [node,elem] = label(node, elem)
            if size(elem{1},1) == 3
                elem = cell2mat(elem')';
            elseif size(elem{1},2) == 3
                elem = cell2mat(elem);
            end
            edgelength(:,1) = sqrt((node(elem(:,1),1)-node(elem(:,2),1)).^2 ...
                            + (node(elem(:,1),2)-node(elem(:,2),2)).^2);

            edgelength(:,2) = sqrt((node(elem(:,2),1)-node(elem(:,3),1)).^2 ...
                            + (node(elem(:,2),2)-node(elem(:,3),2)).^2);

            edgelength(:,3) = sqrt((node(elem(:,3),1)-node(elem(:,1),1)).^2 ...
                            + (node(elem(:,3),2)-node(elem(:,1),2)).^2);

            [~, I] = max(edgelength,[],2);
            elem((I==2),[1 2 3]) = elem((I==2),[2 3 1]);
            elem((I==3),[1 2 3]) = elem((I==3),[3 1 2]);
            elem  = mat2cell(elem,ones(1,size(elem,1)),size(elem,2));
        end
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
        % Minh T.V Nguyen minh.nguyen@ikm.uni-hannover.de
        % Get boundary edges of the domain
        % Initial date: 02/05/2017       
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
        function e = getBoundaryEdges(p,t)
        %BOUNDEDGES Find boundary edges from triangular mesh
        %   E=BOUNDEDGES(P,T)
        %   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.
        % Form all edges, non-duplicates are boundary edges
                edges=[t(:,[1,2]);
                       t(:,[1,3]);
                       t(:,[2,3])];
                node3=[t(:,3);t(:,2);t(:,1)];
                edges=sort(edges,2);
                [~,ix,jx]=unique(edges,'rows');
                vec=histc(jx,1:max(jx));
                qx=find(vec==1);
                e=edges(ix(qx),:);
                node3=node3(ix(qx));

                % Orientation
                v1=p(e(:,2),:)-p(e(:,1),:);
                v2=p(node3,:)-p(e(:,1),:);
                ix=find(v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1)>0);
                e(ix,[1,2])=e(ix,[2,1]);
        end
   end
end
