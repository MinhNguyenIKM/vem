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

%%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
% Initialize 2014
% Modified date: 13/11/1016
% Author: minh.nguyen@ikm.uni-hannover.de
% MINH T.V NGUYEN
% Utility file
%%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
classdef Utility < handle
    methods(Static)
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function arr = createSymArr(symb, m, n)
            for i = 1 : m
                for j = 1 : n
                    arr(i,j) = sym(symb);
                end
            end
        end
        
            
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Added on 22/03/2016
        % This function is to get the elements around the specific nodes
        % Output: elmRel = [Element <--- element index , NodeIndex <--- Node position index
        %                      ...                     ,  ...
        %                      ...                     ,  ... ]
        % Author: Minh T.V Nguyen minh.nguyen@ikm.uni-hannover.de
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function elmRelevant = getElmRelevantNode(nodes, elements)
            % elmRelevant = {[elmIndex1, nodeIndex1;
            %                 elmIndex2, nodeIndex2;
            %                 elmIndex3, nodeIndex3;],
            %                  [...], 
            %                  [...], ... 
            %               }
            elmRelevant = cell(length(nodes),1);
            for n = 1 : length(nodes)
                elmRel = [];
                i = 1;
                for e = 1 : length(elements)
                   elm = elements{e};
                   [yes, index] = ismember(n, elm);
                   if (yes)
                       elmRel(i,1) = e;
                       elmRel(i,2) = index;
                       i = i + 1;
                   end
                end
                elmRelevant{n} = elmRel;
            end
        end

        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Added on 17/08/2016
        % This function is to perform assigning displacement of dof to each
        % node.
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function U = assignUtoNode(Ue,dim)
            U = zeros(length(Ue)/dim,dim);
            for i = 1 : size(U,1)
               for j = 1 : dim 
                   U(i,j) = Ue(i*dim-(dim-j));
               end
            end
        end

        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Added on 31/01/2017
        % This function is to convert Displacement values in 2D (matrix
        % form) to 1D (vector form)
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function U_vector = convertDispIn2Dto1D(disp)
            U_vector = zeros(size(disp,1)*size(disp,2),1);
            for i = 1 : size(disp,1)
                for j = 1 : size(disp,2)
                    if j == 1
                       idx = i*(j+1)-1; 
                    elseif j == 2
                       idx = i*j;
                    end
                    U_vector(idx) = disp(i,j);
                end
            end
        end
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Convert 2D to 1D
        % Author: Minh T.V Nguyen minh.nguyen@ikm.uni-hannover.de
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function vec1D = convert2Dto1D(vec2D)
            vec1D = zeros(size(vec2D,1)*size(vec2D,2),1);
            vec1D(1:2:end) = vec2D(:,1);
            vec1D(2:2:end) = vec2D(:,2);    
        end

        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Initialize this function on 06/10/2016
        % This function is to convert polygonal tessellation in file tess generated
        % by NEPER into vtk file.
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++       
        function [vertices, elements] = getMeshInfoFromTessFile(infile, dim)
            if (isempty(infile))
                error('File does not exist'); 
            end

            fileID = fopen(infile, 'r');
            while (~feof(fileID))
               line = fgetl(fileID);
               if ~isempty(regexp(line,'\*\*vertex', 'once'))
                  line = fgetl(fileID);
                  num_vertex = textscan(line, '%d');
                  if (dim == 2)
                     vertices = zeros(num_vertex{1}, 2); % preallocating array for vertices
                  elseif (dim == 3)
                     vertices = zeros(num_vertex{1}, 3); % preallocating array for vertices
                  end                  
                  while 1
                     line = fgetl(fileID);
                     if ~isempty(regexp(line, '\*\*', 'once')) % If the current line begins with the *, it means that line in turning into another part
                         break 
                     end
                     c = textscan(line, '%d %f %f %f %d');
                     if (dim == 2)
                        vertices(c{1},:) = [c{2}, c{3}];
                     elseif (dim == 3)
                        vertices(c{1},:) = [c{2}, c{3}, c{4}];
                     end
                  end
               end
               
               if ~isempty(regexp(line,'\*\*face', 'once'))
                  line = fgetl(fileID);
                  num_element = textscan(line, '%d');
                  elements = cell(1, num_element{1});
                  while 1
                     line = fgetl(fileID);
                     if ~isempty(regexp(line, '\*\*', 'once')) % If the current line begins with the *, it means that line in turning into another part
                         break 
                     end
                     d = textscan(line, '%d');
                     elements{d{1}(1)} = zeros(1, d{1}(2));
                     for i = 1 : d{1}(2) % d{1}(1) is element id, d{1}(2) is the number of vertices.
                        elements{d{1}(1)}(i) = d{1}(d{1}(2)-i+1+2); % Reverse sequence of number ( 1 2 3 4 --> 4 3 2 1 )
                     end
                     for skip = 1 : 3
                         fgetl(fileID);
                     end
                  end
               end
            end
            fclose(fileID);
        end

        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Added on 17/08/2016
        % This function is to perform converting our format to vtk format
        % (classic) in order to view in Paraview
        %
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function mat2vtk_v1(xy, elms, u, stress, strain, output_filename, title)
            if ( isempty ( output_filename ) )
                output_filename = 'ns2d_fem.vtk';
            end
            nodenum = size(xy,1);
            xyz = zeros(nodenum,3);
            xyz(:,1:2) = xy;
            output_unit = fopen ( output_filename, 'w' );
            fprintf ( output_unit, '# vtk DataFile Version 2.0\n' );
            fprintf ( output_unit, '%s\n', title );
            fprintf ( output_unit, 'ASCII\n' );
            fprintf ( output_unit, 'DATASET UNSTRUCTURED_GRID\n' );
            fprintf ( output_unit, 'POINTS %d float\n', nodenum );
            for node = 1 : nodenum
                fprintf ( output_unit, '  %f  %f  %f\n', xyz(node,:) );
            end
            
            if iscell(elms)
		%
		% Add this condition on 06/10/2016 for polygonal mesh
		%
                cell_size = 0;
                element_num = length(elms);
                for e = 1 : element_num
                    cell_size = cell_size + (length(elms{e}) + 1);
                end
                fprintf ( output_unit, '\n' );
                fprintf ( output_unit, 'CELLS  %d  %d\n', element_num, cell_size );
                for e = 1 : element_num
                   fprintf ( output_unit, '  %d', length(elms{e}) );
                   for i = 1 : length(elms{e})
                      fprintf ( output_unit, '  %d', elms{e}(i)-1 ); 
                   end
                   fprintf ( output_unit, '\n' );
                end
                element_order = 0; %% it means we are in polygonal case
            else
                
                [element_num, element_order] = size(elms);
                cell_size   = element_num * (element_order+1);

                fprintf ( output_unit, '\n' );
                fprintf ( output_unit, 'CELLS  %d  %d\n', element_num, cell_size );
                for element = 1 : element_num
                    fprintf ( output_unit, '  %d', element_order );
                    for order = 1 : element_order
                        fprintf ( output_unit, '  %d', elms(element, order)-1 );
                    end
                    fprintf ( output_unit, '\n' );
                end
            end
            %
            %  VTK has a cell type 22 for quadratic triangles.  However, we
            %  are going to strip the data down to linear triangles for now,
            %  which is cell type 5.
            %
            fprintf ( output_unit, '\n' );
            fprintf ( output_unit, 'CELL_TYPES %d\n', element_num );

            if ( element_order == 3 )
                for element = 1 : element_num
                  fprintf ( output_unit, '5\n' );
                end
            elseif ( element_order == 4 )
                for element = 1 : element_num
                  fprintf ( output_unit, '9\n' );
                end
            elseif ( element_order == 6 )
                for element = 1 : element_num
                  fprintf ( output_unit, '22\n' );
                end
            elseif ( element_order == 0 )
                for element = 1 : element_num
                   fprintf ( output_unit, '7\n' ); 
                end
            end
            
            fprintf ( output_unit, '\n' );
            fprintf ( output_unit, 'POINT_DATA %d\n', nodenum );
            fprintf ( output_unit, 'SCALARS U1 float\n' );
            fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
            for node = 1 : nodenum
                fprintf ( output_unit, '  %.15f\n', u(node,1) );
            end
            
            fprintf ( output_unit, '\n' );
            fprintf ( output_unit, 'SCALARS U2 float\n' );
            fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
            for node = 1 : nodenum
                fprintf ( output_unit, '  %.15f\n', u(node,2) );
            end
            
            fprintf ( output_unit, '\n' );
            fprintf ( output_unit, 'CELL_DATA %d\n', element_num );
            fprintf ( output_unit, 'SCALARS S11 float\n' );
            fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
            for e = 1 : element_num
                fprintf ( output_unit, '  %.15f\n', stress{e}(1,1) ); %node 1, S11
            end
            
            fprintf ( output_unit, '\n' );
            fprintf ( output_unit, 'SCALARS S22 float\n' );
            fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
            for e = 1 : element_num
                fprintf ( output_unit, '  %.15f\n', stress{e}(2,1) ); %node 1, S22
            end
            
            fprintf ( output_unit, '\n' );
            fprintf ( output_unit, 'SCALARS S12 float\n' );
            fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
            for e = 1 : element_num
                fprintf ( output_unit, '  %.15f\n', stress{e}(3,1) ); %node 1, S12
            end
            
            fprintf ( output_unit, '\n' );
            fprintf ( output_unit, 'SCALARS E11 float\n' );
            fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
            for e = 1 : element_num
                fprintf ( output_unit, '  %.15f\n', strain{e}(1,1) ); %node 1, E11
            end
            
            fprintf ( output_unit, '\n' );
            fprintf ( output_unit, 'SCALARS E22 float\n' );
            fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
            for e = 1 : element_num
                fprintf ( output_unit, '  %.15f\n', strain{e}(2,1) ); %node 1, E22
            end
            
            fprintf ( output_unit, '\n' );
            fprintf ( output_unit, 'SCALARS E12 float\n' );
            fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
            for e = 1 : element_num
                fprintf ( output_unit, '  %.15f\n', strain{e}(3,1) ); %node 1, E12
            end
            fclose(output_unit);
        end
       
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Added on 02/12/2016 based on the old version mat2vtk_v1
        % This function is to get all information about a mesh which is
        % generated by ABAQUS. 
        % It will return information such as vertices, elementsm
        % DirichletBC nodes, NeumannBC nodes.
        %
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function mat2vtk(xy, elms, u, stress, strain, output_filename, title)
            if ( isempty ( output_filename ) )
                output_filename = 'ns2d_fem.vtk';
            end
            nodenum = size(xy,1);
            xyz = zeros(nodenum,3);
            xyz(:,1:2) = xy;
            output_unit = fopen ( output_filename, 'w' );
            fprintf ( output_unit, '# vtk DataFile Version 2.0\n' );
            fprintf ( output_unit, '%s\n', title );
            fprintf ( output_unit, 'ASCII\n' );
            fprintf ( output_unit, 'DATASET UNSTRUCTURED_GRID\n' );
            fprintf ( output_unit, 'POINTS %d float\n', nodenum );
            for node = 1 : nodenum
                fprintf ( output_unit, '  %f  %f  %f\n', xyz(node,:) );
            end
            
            if iscell(elms)
		%
		% Add this condition on 06/10/2016 for polygonal mesh
		%
                cell_size = 0;
                element_num = length(elms);
                for e = 1 : element_num
                    cell_size = cell_size + (length(elms{e}) + 1);
                end
                fprintf ( output_unit, '\n' );
                fprintf ( output_unit, 'CELLS  %d  %d\n', element_num, cell_size );
                for e = 1 : element_num
                   fprintf ( output_unit, '  %d', length(elms{e}) );
                   for i = 1 : length(elms{e})
                      fprintf ( output_unit, '  %d', elms{e}(i)-1 ); 
                   end
                   fprintf ( output_unit, '\n' );
                end
                element_order = 0; %% it means we are in polygonal case
            else
                
                [element_num, element_order] = size(elms);
                cell_size   = element_num * (element_order+1);

                fprintf ( output_unit, '\n' );
                fprintf ( output_unit, 'CELLS  %d  %d\n', element_num, cell_size );
                for element = 1 : element_num
                    fprintf ( output_unit, '  %d', element_order );
                    for order = 1 : element_order
                        fprintf ( output_unit, '  %d', elms(element, order)-1 );
                    end
                    fprintf ( output_unit, '\n' );
                end
            end
            %
            %  VTK has a cell type 22 for quadratic triangles.  However, we
            %  are going to strip the data down to linear triangles for now,
            %  which is cell type 5.
            %
            fprintf ( output_unit, '\n' );
            fprintf ( output_unit, 'CELL_TYPES %d\n', element_num );

            if ( element_order == 3 )
                for element = 1 : element_num
                  fprintf ( output_unit, '5\n' );
                end
            elseif ( element_order == 4 )
                for element = 1 : element_num
                  fprintf ( output_unit, '9\n' );
                end
            elseif ( element_order == 6 )
                for element = 1 : element_num
                  fprintf ( output_unit, '22\n' );
                end
            elseif ( element_order == 0 )
                for element = 1 : element_num
                   fprintf ( output_unit, '7\n' ); 
                end
            end
            
            fprintf ( output_unit, '\n' );
            fprintf ( output_unit, 'POINT_DATA %d\n', nodenum );
            fprintf ( output_unit, 'SCALARS U1 float\n' );
            fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
            for node = 1 : nodenum
                fprintf ( output_unit, '  %.15f\n', u(node,1) );
            end
            
            fprintf ( output_unit, '\n' );
            fprintf ( output_unit, 'SCALARS U2 float\n' );
            fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
            for node = 1 : nodenum
                fprintf ( output_unit, '  %.15f\n', u(node,2) );
            end
            
            fprintf ( output_unit, '\n' );
            fprintf ( output_unit, 'CELL_DATA %d\n', element_num );
            fprintf ( output_unit, 'SCALARS S11 float\n' );
            fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
            for e = 1 : element_num
                fprintf ( output_unit, '  %.15f\n', stress{e}(1,1) ); %node 1, S11
            end
            
            fprintf ( output_unit, '\n' );
            fprintf ( output_unit, 'SCALARS S22 float\n' );
            fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
            for e = 1 : element_num
                fprintf ( output_unit, '  %.15f\n', stress{e}(2,1) ); %node 1, S22
            end
            
            fprintf ( output_unit, '\n' );
            fprintf ( output_unit, 'SCALARS S12 float\n' );
            fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
            for e = 1 : element_num
                fprintf ( output_unit, '  %.15f\n', stress{e}(3,1) ); %node 1, S12
            end
            
            fprintf ( output_unit, '\n' );
            fprintf ( output_unit, 'SCALARS E11 float\n' );
            fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
            for e = 1 : element_num
                fprintf ( output_unit, '  %.15f\n', strain{e}(1,1) ); %node 1, E11
            end
            
            fprintf ( output_unit, '\n' );
            fprintf ( output_unit, 'SCALARS E22 float\n' );
            fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
            for e = 1 : element_num
                fprintf ( output_unit, '  %.15f\n', strain{e}(2,1) ); %node 1, E22
            end
            
            fprintf ( output_unit, '\n' );
            fprintf ( output_unit, 'SCALARS E12 float\n' );
            fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
            for e = 1 : element_num
                fprintf ( output_unit, '  %.15f\n', strain{e}(3,1) ); %node 1, E12
            end
            fclose(output_unit);
        end
        

        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function [vertices, elements] = getMeshInfoFromABAQUSFile_v1(filename, meshtype)
            fileID = fopen(filename, 'r');
            while (~feof(fileID))
               line = fgetl(fileID);
               if ~isempty(regexp(line,'^\*Node', 'once')) && isempty(regexp(line,'Output', 'once')) % This line marks the group of Nodes
                  while (isempty(regexp(line,'^\*Element', 'once'))) % If the current line begins with the *, it means that line in turning into another part
                      line = fgetl(fileID);
                      c = textscan(line, '%d, %f, %f'); % Coordinates of Nodes in Column 2 and Column 3
                      vertices(c{1},:) = [c{2} c{3}];
                  end
               end
               if ~isempty(regexp(line, '^\*Element', 'once')) && isempty(regexp(line, 'Output', 'Once')) % This line marks the group of Elements                  
                  while (isempty(regexp(line, '^\*Nset', 'once'))) % If the current line begins with the *, it means that line in turning into another part
                      if strcmpi(meshtype,'T3')
                        d = textscan(line, '%d, %d, %d, %d');
                        elements{d{1}} = [d{2} d{3} d{4}];
                      elseif strcmpi(meshtype, 'Q4')
                        d = textscan(line, '%d, %d, %d, %d, %d');
                        elements{d{1}} = [d{2} d{3} d{4} d{5}];
                      end
                      line = fgetl(fileID);
                  end
               end
            end
            fclose(fileID);
        end
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Modified on 02/12/2016 based on the old version getMeshInfoFromABAQUSFile_v1
        % This function modified because of the needs for collecting
        % elements and nodes in mixed elements mesh especially in crack
        % analysis.
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function [vertices, elements] = getMeshInfoFromABAQUSFile(filename, meshtype)
            fileID = fopen(filename, 'r');
            meshtype_flag = meshtype;
            while (~feof(fileID))
               line = fgetl(fileID);
               if ~isempty(regexp(line,'^\*Node', 'once')) && isempty(regexp(line,'Output', 'once')) % This line marks the group of Nodes
                  while (isempty(regexp(line,'^\*Element', 'once'))) % If the current line begins with the *, it means that line in turning into another part
                      line = fgetl(fileID);
                      c = textscan(line, '%d, %f, %f'); % Coordinates of Nodes in Column 2 and Column 3
                      vertices(c{1},:) = [c{2} c{3}];
                  end
               end
               if ~isempty(regexp(line, '^\*Element', 'once')) && isempty(regexp(line, 'Output', 'Once')) % This line marks the group of Elements                   
                  while (isempty(regexp(line, '^\*Nset|^\*End', 'once'))) % If the current line begins with the *, it means that line in turning into another part
                      if strcmpi(meshtype,'T3')
                        d = textscan(line, '%d, %d, %d, %d');
                        elements{d{1}} = [d{2} d{3} d{4}];
                      elseif strcmpi(meshtype, 'Q4')
                        d = textscan(line, '%d, %d, %d, %d, %d');
                        elements{d{1}} = [d{2} d{3} d{4} d{5}];
                      end
                      %--------- These couple of codes for mixed-element
                      % mesh which is used for crack analysis ----------
                      if strcmpi(meshtype_flag, 'MIXED') && ~isempty(regexp(line, 'CPS4|CPE4|CPE4R', 'once')) % CPS4 in ABAQUS means Q4 elements - linear case
                          meshtype = 'Q4';
                      elseif strcmpi(meshtype_flag, 'MIXED') && ~isempty(regexp(line, 'CPS3|CPE3|CPE3R', 'once')) % CPS3 in ABAQUS means T3 elements - linear case
                          meshtype = 'T3';
                      end
                      %----------------------------------
                      line = fgetl(fileID);
                      %----------------------------------
                      if strcmpi(meshtype_flag, 'MIXED') && ~isempty(regexp(line, 'CPS4|CPE4|CPE4R', 'once')) % CPS4 in ABAQUS means Q4 elements - linear case
                          meshtype = 'Q4';
                          line = fgetl(fileID);
                      elseif strcmpi(meshtype_flag, 'MIXED') && ~isempty(regexp(line, 'CPS3|CPE3|CPE3R', 'once')) % CPS3 in ABAQUS means T3 elements - linear case
                          meshtype = 'T3';
                          line = fgetl(fileID);
                      end
                      %------------------------------------------------
                  end
               end
            end
            fclose(fileID);
        end
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Initial date: 17/05/2017
        % Author: H. Nguyen-Xuan
        % a part of PolyTree code
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function [x,y] = getEllipse(gcoord_edge,e)
            %Let (x1,y1) and (x2,y2) be the coordinates of the two vertices of the ellipse's major axis, 
            %and let e be its eccentricity.
            x1=gcoord_edge(1,1);
            x2=gcoord_edge(2,1);
            y1=gcoord_edge(1,2);
            y2=gcoord_edge(2,2);

            a = 1/2*sqrt((x2-x1)^2+(y2-y1)^2);
            b = a*sqrt(1-e^2);
            t = linspace(0,2*pi,11);
            X = a*cos(t);
            Y = b*sin(t);
            w = atan2(y2-y1,x2-x1);
            x = (x1+x2)/2 + X*cos(w) - Y*sin(w);
            y = (y1+y2)/2 + X*sin(w) + Y*cos(w);
        end
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Initial date: 18/05/2017
        % Author: H. Nguyen-Xuan
        % a part of PolyTree code
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function [nodesConnectivity] = addHangingNode(nodesConnectivity, nodesEdge, nodeAdded)
            % Xen nut moi vao phan tu chung canh voi phan tu muon adaptive
            i = 1;
            while i <= length(nodesConnectivity)
                j = 1;
                while j <= 2
                    if nodesConnectivity(i) == nodesEdge(j)
                        if i == 1 && j == 1 && nodesConnectivity(2) ~= nodesEdge(2)
                            left = nodesConnectivity(1:end);
                            right = [];
                        elseif i == 1 && j == 2 && nodesConnectivity(2) ~= nodesEdge(1)
                            left = nodesConnectivity(1:end);
                            right = [];
                        else
                            left = nodesConnectivity(1:i);
                            right = nodesConnectivity(i+1:end);
                        end
                        j = 3;
                        i = length(nodesConnectivity)+1;
                    end
                    j = j + 1;
                end
                i = i + 1;
            end
            nodesConnectivity = [left, nodeAdded, right];  % add new nodes
        end
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function [idx,D] = knnsearch(varargin)
            % KNNSEARCH   Linear k-nearest neighbor (KNN) search
            % IDX = knnsearch(Q,R,K) searches the reference data set R (n x d array
            % representing n points in a d-dimensional space) to find the k-nearest
            % neighbors of each query point represented by eahc row of Q (m x d array).
            % The results are stored in the (m x K) index array, IDX. 
            %
            % IDX = knnsearch(Q,R) takes the default value K=1.
            %
            % IDX = knnsearch(Q) or IDX = knnsearch(Q,[],K) does the search for R = Q.
            %
            % Rationality
            % Linear KNN search is the simplest appraoch of KNN. The search is based on
            % calculation of all distances. Therefore, it is normally believed only
            % suitable for small data sets. However, other advanced approaches, such as
            % kd-tree and delaunary become inefficient when d is large comparing to the
            % number of data points. On the other hand, the linear search in MATLAB is
            % relatively insensitive to d due to the vectorization. In  this code, the 
            % efficiency of linear search is further improved by using the JIT
            % aceeleration of MATLAB. Numerical example shows that its performance is
            % comparable with kd-tree algorithm in mex.
            %
            % See also, kdtree, nnsearch, delaunary, dsearch

            % By Yi Cao at Cranfield University on 25 March 2008

            % Example 1: small data sets
            %{
            R=randn(100,2);
            Q=randn(3,2);
            idx=knnsearch(Q,R);
            plot(R(:,1),R(:,2),'b.',Q(:,1),Q(:,2),'ro',R(idx,1),R(idx,2),'gx');
            %}

            % Example 2: ten nearest points to [0 0]
            %{
            R=rand(100,2);
            Q=[0 0];
            K=10;
            idx=knnsearch(Q,R,10);
            r=max(sqrt(sum(R(idx,:).^2,2)));
            theta=0:0.01:pi/2;
            x=r*cos(theta);
            y=r*sin(theta);
            plot(R(:,1),R(:,2),'b.',Q(:,1),Q(:,2),'co',R(idx,1),R(idx,2),'gx',x,y,'r-','linewidth',2);
            %}

            % Example 3: cputime comparion with delaunay+dsearch I, a few to look up
            %{
            R=randn(10000,4);
            Q=randn(500,4);
            t0=cputime;
            idx=knnsearch(Q,R);
            t1=cputime;
            T=delaunayn(R);
            idx1=dsearchn(R,T,Q);
            t2=cputime;
            fprintf('Are both indices the same? %d\n',isequal(idx,idx1));
            fprintf('CPU time for knnsearch = %g\n',t1-t0);
            fprintf('CPU time for delaunay  = %g\n',t2-t1);
            %}
            % Example 4: cputime comparion with delaunay+dsearch II, lots to look up
            %{
            Q=randn(10000,4);
            R=randn(500,4);
            t0=cputime;
            idx=knnsearch(Q,R);
            t1=cputime;
            T=delaunayn(R);
            idx1=dsearchn(R,T,Q);
            t2=cputime;
            fprintf('Are both indices the same? %d\n',isequal(idx,idx1));
            fprintf('CPU time for knnsearch = %g\n',t1-t0);
            fprintf('CPU time for delaunay  = %g\n',t2-t1);
            %}
            % Example 5: cputime comparion with kd-tree by Steven Michael (mex file) 
            % <a href="http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=7030&objectType=file">kd-tree by Steven Michael</a> 
            %{
            Q=randn(10000,10);
            R=randn(500,10);
            t0=cputime;
            idx=knnsearch(Q,R);
            t1=cputime;
            tree=kdtree(R);
            idx1=kdtree_closestpoint(tree,Q);
            t2=cputime;
            fprintf('Are both indices the same? %d\n',isequal(idx,idx1));
            fprintf('CPU time for knnsearch = %g\n',t1-t0);
            fprintf('CPU time for delaunay  = %g\n',t2-t1);
            %}

            % Check inputs
            [Q,R,K,fident] = Utility.parseinputs(varargin{:});

            % Check outputs
            nargoutchk(0,2);

            % C2 = sum(C.*C,2)';
            [N,M] = size(Q);
            L=size(R,1);
            idx = zeros(N,K);
            D = idx;

            if K==1
                % Loop for each query point
                for k=1:N
                    d=zeros(L,1);
                    for t=1:M
                        d=d+(R(:,t)-Q(k,t)).^2;
                    end
                    if fident
                        d(k)=inf;
                    end
                    [D(k),idx(k)]=min(d);
                end
            else
                for k=1:N
                    d=zeros(L,1);
                    for t=1:M
                        d=d+(R(:,t)-Q(k,t)).^2;
                    end
                    if fident
                        d(k)=inf;
                    end
                    [s,t]=sort(d);
                    idx(k,:)=t(1:K);
                    D(k,:)=s(1:K);
                end
            end
            if nargout>1
                D=sqrt(D);
            end
        end

        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function [Q,R,K,fident] = parseinputs(varargin)
            % Check input and output
            narginchk(1,3);

            Q=varargin{1};

            if nargin<2
                R=Q;
                fident = true;
            else
                fident = false;
                R=varargin{2};
            end

            if isempty(R)
                fident = true;
                R=Q;
            end

            if ~fident
                fident = isequal(Q,R);
            end

            if nargin<3
                K=1;
            else
                K=varargin{3};
            end
        end

        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Added on 19/04/2015
        % Modified function getCoordinateFromFile function with 2 input
        % values
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function nodeCoor = getCoordinatesFromFile(filename, nodeCoordOld)
            if (nargin == 2)
                nodeCoor = nodeCoordOld;
            elseif (nargin == 1)
                nodeCoor = [];
            end
            fileID = fopen(filename, 'r');
            while (~feof(fileID))
                line = fgetl(fileID);
                if isempty(line) || ~isempty(regexp(line,'^\%*'))
                    continue;
                end
                c = textscan(line,'%d %f %f');  % Coordinates of Points in Column-2 Column-3
                disp(c{1});
                nodeCoor(c{1}, :) = [c{2} c{3}];
            end
            fclose(fileID);
        end

        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Added on 19/04/2015
        % This function will renumber node along the interface.
        % An array will be given through a file that map the old numbers
        % to the new numbers
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function [reNodeCoor, rePlanes] = reNumberNode(filename, nodeCoor, planes)
            reNodeCoor = nodeCoor;
            rePlanes   = planes;
            fileID = fopen(filename, 'r');
            map = [];
            while (~feof(fileID));
                line = fgetl(fileID);
                if isempty(line) || ~isempty(regexp(line, '^\%*'))
                    continue;
                end
                c = textscan(line, '%d %d');  %Old number in Column-1, new number in Column-2
                map(end+1,:) = [c{1} c{2}];
            end
            fclose(fileID);
            for i = 1 : size(map,1)
                % Rearrange the order of the node coordinate
%                 tempCoor = reNodeCoor(map(i,1),:);
%                 reNodeCoor(map(i,1),:) = reNodeCoor(map(i,2),:);
%                 reNodeCoor(map(i,2),:) = tempCoor;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % As a result, the enumerate number of each plane must be
                % changed
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for m = 1 : size(planes,1)
                    for n = 1 : size(planes,2)
                        if rePlanes(m,n) == map(i,1)
                            rePlanes(m,n) = map(i,2);
                        elseif rePlanes(m,n) == map(i,2)
                            rePlanes(m,n) = map(i,1);
                        end
                    end
                end
            end           
        end

        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function elements = getElementsFromFile(filename)
            % Example structure of elements
            % elements = [1      2      3
            %             4      5      6
            %             7      8      9
            %             10     11     12]
            elements = [];
            fileID = fopen(filename, 'r');
            while(~feof(fileID))
                line = fgetl(fileID);
                if isempty(line)
                    continue;
                end
                c = textscan(line, '%d %d %d %d');
                elements(c{1},:) = [c{2} c{3} c{4}];
            end
            fclose(fileID);
        end

%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function Inodes = getInterfaceNodesFromFile(filename)
            fileID = fopen(filename, 'r');
            start = 0;
            Inodes = [];
            while(~feof(fileID))
                line = fgetl(fileID);
                if isempty(line)
                   continue;
                end
                c = textscan(line, '%d');
                for i = 1 : length(c{1})
                    Inodes(i+start) = c{1}(i);
                end
                start = length(Inodes);
            end
        end
        
%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
           
%         function planes = getPlanesFromFile(filename)
%             % Example structure of planes
%             % planes = [1	4	5	2
%             %           2   5   6	3
%             %           4	7	8	5  
%             %           5	8	9	6
%             %           10  8   1   7]
%             planes = [];
%             fileID = fopen(filename, 'r');
%             while(~feof(fileID))
%                 line = fgetl(fileID);
%                 if isempty(line) || ~isempty(regexp(line,'^\%*'))
%                     continue;
%                 end
%                 c = textscan(line, '%d %d %d %d %d');
%                 planes(end+1,:) = [c{2} c{3} c{4} c{5}]; % PointID in Column-2 Column-3 Column-4 Column-5  
%             end
%             fclose(fileID);
%         end

%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Added on 19/04/2015
        % Modified function getPlanesFromFile2 function with 2 input
        % values
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function planes = getPlanesFromFile(filename, planes1) 
            % Example structure of planes
            % planes = [1	4	5	2
            %           2   5   6	3
            %           4	7	8	5  
            %           5	8	9	6
            %           10  8   1   7]
            if (nargin == 2)
                planes = planes1;
            elseif (nargin == 1)
                planes = [];
            end
            fileID = fopen(filename, 'r');
            while(~feof(fileID))
                line = fgetl(fileID);
                if isempty(line) || ~isempty(regexp(line,'^\%*'))
                    continue;
                end
                c = textscan(line, '%d %d %d %d %d');
                planes(c{1},:) = [c{2} c{3} c{4} c{5}]; % PointID in Column-2 Column-3 Column-4 Column-5  
            end
            fclose(fileID);
        end


        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function [nodes, planes] = generateQ4Mesh(filenode, fileplane)
            nodes = Utility.getCoordinatesFromFile(filenode);
            planes = Utility.getNodesOfPlane(fileplane);
        end

        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function bcc = applyBCC(filename)
            fileID = fopen(filename, 'r');
            start = 0;
            bcc = [];
            while(~feof(fileID))
                line = fgetl(fileID);
                if isempty(line) || ~isempty(regexp(line,'^\%*'))
                   continue;
                end
                c = textscan(line, '%d');
                for i = 1 : length(c{1})
                    bcc(i+start) = c{1}(i);
                end
                start = length(bcc);
            end
        end

          
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function [scaledU, scaledT, scaledCoordinates, E] = scaling(nodeCoordinates, YoungModulus)
            % Scale the system coordinates (Catersian coordinate system) to
            % unity scaled. That means the largest distance between 2 nodes
            % in the system is 1.
            maximum(1) = max(nodeCoordinates(:,1));
            maximum(2) = max(nodeCoordinates(:,2));
            minimum(1) = min(nodeCoordinates(:,1));
            minimum(2) = min(nodeCoordinates(:,2));
            scaledU = max(maximum - minimum);
            scaledT = 1/YoungModulus;
            scaledCoordinates = 1/scaledU * nodeCoordinates;
            E = 1;
        end
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %
        %
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function t = removeScaledT(st, scaledT)
            t = st / scaledT;
        end

        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function u = removeScaledU(su, scaledU)
            u = su * scaledU;
        end
 
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function [uF, tB, uB, uBI] = removingScales(u_F, t_B, uB_boundary, uB_inside, scaledT, scaledU)    
            uF = Utility.removeScaledU(u_F, scaledU);
            tB = Utility.removeScaledT(t_B, scaledT);
            uB = Utility.removeScaledU(uB_boundary, scaledU);
            uBI = Utility.removeScaledU(uB_inside, scaledU);
        end
 
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %
        %
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function new = copy(this)
            save('temp.mat', 'this');
            Foo = load('temp.mat');
            new = Foo.this;
            delete('temp.mat');
        end
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Initial date: 05.09.2018
        % Author: minh.nguyen@ikm.uni-hannover.de
        % create polygon counter-clockwise
        % Ben Voigt's algorithm. https://stackoverflow.com/questions/13935324/sorting-clockwise-polygon-points-in-matlab
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function [xo,yo] = createPolygonCW(x,y)
%                 Step 1: Find the centroid:
                cx = mean(x);
                cy = mean(y);
%                 Step 2: Find the angles:
                a = atan2(y - cy, x - cx);
%                 Step 3: Find the correct sorted order:
                [~, order] = sort(a);
%                 Step 4: Reorder the coordinates:
                xo = x(order);
                yo = y(order);
        end
    end
end
