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

classdef VirtualElementMethod
    %VIRTUALELEMENTMETHOD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private, GetAccess = public)
        material_;
        method_; % method -- name_
                 %        -- dim_
                 %        -- k_     (1 is Linear, 2 is Quadratic Apprx
                 %        -- TypeF_ ('traction' or 'concentrated')
        mesh_;
        grid_;
        dirichlet_;  % dirichlet is either an arraycell or an array
        neumann_; % and neumann is either an arraycell or an array
        uFunc_;
        tFunc_;
    end
    
    methods (Access = public)
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Constructor of Virtual Element Method
        % Minh T.V Nguyen minh.nguyen@ikm.uni-hannover.de
        % Initial date: 13/04/2017
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function this = VirtualElementMethod(material, method, mesh, dirichlet, neumann, uD, tN)
            disp('---> Constructor of VEM'); 
            this.material_ = material;
            this.method_   = method;
            this.mesh_     = mesh;
            this.dirichlet_ = dirichlet;
            this.neumann_   = neumann;
            this.uFunc_    = uD; % @hanle function for displacements
            this.tFunc_    = tN; % @handle function for tractions
%             this.grid_     = MeshOperation.get_MeshInfo(mesh);
        end
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Solving structure in 2D Linear Elasticity
        % Minh T.V Nguyen minh.nguyen@ikm.uni-hannover.de
        % Date: 11/04/2017
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function [U, PiNstar] = solveLinearElastic2D(this)
            dim     = this.method_.dim_;  
            elements = this.mesh_.elements;
            nodes  = this.mesh_.vertices;
            nodeNum = size(nodes,1);
            K = sparse(dim*nodeNum,dim*nodeNum); % LHS - Stiffness matrix
            F = zeros(dim*nodeNum,1); % RHS - Force summation
            U = zeros(dim*nodeNum,1); % Unknown field - Displacement
            PiNstar = cell(1, size(elements,1)); % Coefficient of Projection of Ansatz funtions
                    
            [K, F, U, PiNstar, inDOFs, bDOFs] = this.solveEL2DGlobal(K, F, U, PiNstar);
            %------------- Solving Equations ------------------------------
            if ~isempty(bDOFs)
                F_boundary = K(inDOFs,bDOFs)*U(bDOFs);
            else
                F_boundary = zeros(length(F(inDOFs)),1);
            end
            if ~isempty(inDOFs)
                U(inDOFs) = K(inDOFs,inDOFs) \ (F(inDOFs) - F_boundary);
            end
            %--------------------------------------------------------------
        end
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Assembly whole structure
        % Minh T.V Nguyen minh.nguyen@ikm.uni-hannover.de
        % Date: 11/04/2017
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function [K, F, U, PiNstar, inDOFs, bDOFs] = solveEL2DGlobal(this, K, F, U, PiNstar)
            elements  = this.mesh_.elements;
            vertices  = this.mesh_.vertices; 
            dim       = this.method_.dim_;
            nodeNum   = length(vertices);

            % CALCULATE STIFFNESS MATRIX OF EACH ELEMENT & THEN ASSEMBLY TO GLOBAL SYSTEM
            for e = 1 : length(elements)
               [Ke, Fe, PiNstar0] = this.solveEL2DLocal(elements{e}, vertices(elements{e},:)); 
               dofK = this.node2dof(elements{e}, 'element'); 
               K(dofK, dofK) = K(dofK, dofK) + Ke;
               F(dofK, 1) = F(dofK, 1) + Fe;
               PiNstar{e} = PiNstar0;
            end
            % APPLY BOUNDARY CONDITIONS %

%             [dirichletBC, neumannBC] = this.assignBoundaryConditions();
%             bDOFs2     = this.node2dof_old(dirichletBC, 'boundary'); % Boundary Degrees of freedom
            bDOFs     = this.node2dof(this.dirichlet_, 'boundary'); % Boundary Degrees of freedom
            % SET SUBSTRACTION: B\A which B is a big set and A is a small set.
            inDOFs = setdiff(1:dim*nodeNum, bDOFs);
%             F      = this.applyNeumann(F, neumannBC);
%             U      = this.applyDirichlet(U, dirichletBC);
            F = this.applyNeumannBC(F);
            U = this.applyDirichletBC(U);
        end
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Solving for each discretization, each element
        % Minh T.V Nguyen minh.nguyen@ikm.uni-hannover.de
        % Date: 11/04/2017
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function [Ke, Fe, PiNstar0] = solveEL2DLocal(this, element, coordinates)
            vert = coordinates;
            k = this.method_.k_;
            dim = this.method_.dim_;
            thickness = this.material_.getThickness();
            E = this.material_.getMaterial();
            % GET AREA & CENTROID of an element
            [area, centroid, diameter] = MeshOperation.get_ElementInfo(vert);

            Ndof = length(element);
            nk = k*3*dim;
            B = cell(nk, Ndof);       
            D = cell(Ndof, nk);
            
            for i = 1:Ndof % Loop over the NDOF. 
                % get Basis of Polynomials
                % get Gradient of Basis
                [M, DM] = Basis(k, vert(i,1), vert(i,2), centroid(1), centroid(2), diameter);
                for alpha = 1 : nk
                    if alpha == 1
                        D{i,alpha} = [1 0]';
                        B{alpha,i} = 1/length(element) * [1 0];
                    elseif alpha == 2
                        D{i,alpha} = [0 1]';
                        B{alpha,i} = 1/length(element) * [0 1];
                    else
                        D{i,alpha} = this.getD(M{alpha});
                        B{alpha,i} = this.getB(DM{alpha}, vert, Ndof, E, i, alpha, area);
                    end
                end
            end
            
            % B_matrix is a matrix 2*(k*3 x Ndof)
            B_matrix = cell2mat(B);
            % D_matrix is a matrix 2*(Ndof x k*3)
            D_matrix = cell2mat(D);
            % G_matrix is a matrix 2*(k*3 x k*3)
            G_matrix = B_matrix * D_matrix ;
            G_matrixtilde = G_matrix;
            % G_tilde(1:2,:) = 0;
            G_matrixtilde(1:3,:) = 0;
            
            % a is a coefficient matrix  2*(3k x Ndof). a is moreover coefficent
            % of Projection in Polynomial space Pk
            a = linsolve(G_matrix,B_matrix);
            
            % Pi (Projection operator in Element space) 2*(Ndof x Ndof)
            Pi = D_matrix * a;
            
            % Consistency Matrix dimension is 2*(Ndof x Ndof)
            Ke_c = a' * G_matrixtilde * a;
            
            % Stabilization Matrix
            % Stabilization coefficient is chosen from "On the use of the
            % virtual element method for geomechanics on reservoir grids"
            alpha = 1/10*abs(area)*thickness*trace(E)*trace(inv(D_matrix(:,3:6)'*D_matrix(:,3:6)));
            Ke_s = alpha * ((eye(dim*Ndof) - Pi)' * (eye(dim*Ndof) - Pi));
            
            % Stiffness matrix element
            Ke = Ke_c + Ke_s;
            
            % Force: We assume no body force
            Fe = zeros(size(Ke,1),1);
            
            % Projection of Shape Functions
            PiNstar0 = a;
        end
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Minh T.V Nguyen minh.nguyen@ikm.uni-hannover.de
        % 24/04/2017
        % move function postProcessing to a part of the
        % VirtualElementMethod class
        % It will perform the post processing process.
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function [Strain_elms, Stress_elms] = executePostProcessing(this, a_PiN, U)
            mesh = this.mesh_;
            k = this.method_.k_;
            material = this.material_.getMaterial();
            vertices = mesh.vertices;
            elements = mesh.elements;
            for e = 1 : length(elements)
               elm = elements{e};
               Ndof = length(elm);
               vert = vertices(elm,:);
               Epsilon = zeros(3,Ndof);
               Sigma = zeros(3,Ndof);
               [area, centroid, diameter] = MeshOperation.get_ElementInfo(vert);
               % need a loop over vertices here instead of i:=1 always
               i = 1; % i = 1 : the number of node (But in this case, DM is always unchangable because we are using linear approximation)
               [M, DM] = Basis(k, vert(i,1), vert(i,2), centroid(1), centroid(2), diameter);
               Ue = this.getDisplacementOfElement(U, elm);
               GradDisp = 0;
               for alpha = 1 : length(DM)
                   GradDisp = GradDisp + dot(a_PiN{e}(alpha,:), Ue) .* DM{alpha};
               end
               StrainVector = [GradDisp(1,1); GradDisp(2,2); GradDisp(1,2) + GradDisp(2,1)];
               for point = 1 : Ndof
                   Epsilon(:,point) = StrainVector;
                   Sigma(:,point) = material * StrainVector;
               end
               Strain_elms{e} = Epsilon;
               Stress_elms{e} = Sigma;
            end
        end
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % minh.nguyen@ikm.uni-hannover.de
        % 14/04/2017 (Eastern Day)
        % Destructor of VEM Class
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function delete(this)
            disp('---> Destructor of VEM'); 
        end
    end
    
    methods (Access = private)
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % minh.nguyen@ikm.uni-hannover.de
        % 13/04/2017
        % Initialize this function assignBoundaryConditions
        % Apply described u function for dirichlet nodes and described t
        % function for neumann nodes
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function [dirichletBC, neumannBC] = assignBoundaryConditions(this)
            uNodes = this.dirichlet_; 
            tNodes   = this.neumann_;   
            uFunction = this.uFunc_;
            tFunction = this.tFunc_;
            vertices  = this.mesh_.vertices;
            if iscell(uNodes) && iscell(tNodes)
                for i = 1 : length(uFunction)
                    for n = 1 : length(uNodes{i})
                        dirichletBC{i}(n,1)   = uNodes{i}(n);
                        dirichletBC{i}(n,2:3) = uFunction{i}(vertices(uNodes{i}(n),1),vertices(uNodes{i}(n),2)); % Substitute coordinates of prescribed values to the displacement function
                    end
                end
                for i = 1 : length(tFunction)
                    for n = 1 : length(tNodes{i})
                        neumannBC{i}(n,1)   = tNodes{i}(n);
                        neumannBC{i}(n,2:3) = tFunction{i}(vertices(tNodes{i}(n),1),vertices(tNodes{i}(n),2)); % Substitute coordinates of prescribed values to the traction function
                    end
                end
            elseif ~iscell(uNodes) && ~iscell(tNodes)
                for n = 1 : length(uNodes)
                    dirichletBC(n,1)   = uNodes(n);
                    dirichletBC(n,2:3) = uFunction(vertices(uNodes(n),1),vertices(uNodes(n),2)); % Substitute coordinates of prescribed values to the displacement function
                end
                for n = 1 : length(tNodes)
                    neumannBC(n,1)   = tNodes(n);
                    neumannBC(n,2:3) = tFunction(vertices(tNodes(n),1),vertices(tNodes(n),2)); % Substitute coordinates of prescribed values to the traction function
                end
            else
               error('Should be aware of this case');
            end
        end
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % minh.nguyen@ikm.uni-hannover.de
        % Modified date: 13/04/2017
        % Update function getD to a method this.getD of class
        % VirtualElementMethod
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function D_alphai = getD(this, m)
           %m: element alpha-th m of BASIS M 
            D_alphai = m;
        end
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % minh.nguyen@ikm.uni-hannover.de
        % Modified date: 13/04/2017
        % Update the old function getD to a method this.getB of class
        % VirtualElementMethod
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function B_alphai = getB(this, dm, vert, Ne, mat, i, alpha, area)
        % calculate each element of matrix B
        %       m: ith element of polynomial basis
        %       vert: vertices of this element
        %       mat: material
            thickness = this.material_.thickness_;
            k = this.method_.k_;
            strain = [dm(1,1) dm(2,2) dm(1,2)+dm(2,1)]';
            stressvec = mat * strain; % matrix multiplication
            stress = [stressvec(1) stressvec(3); stressvec(3) stressvec(2)];
            Phi = [1 0; 0 1]; % Phi(i) at the specific vertice
            if i == 1
                P_prev = vert(Ne,:);
                P_next = vert(i+1,:);
            elseif i == Ne
                P_prev = vert(i-1,:);
                P_next = vert(1,:);
            else
                P_prev = vert(i-1,:);
                P_next = vert(i+1,:);
            end
            P = vert(i,:);
            if alpha ~= 3 % khong phai la rotation (rigid body motion)
                % Calculation using Gauss-Lobatto
                % normalvector = [(y2 - y1), -(x2-x1)];
                B_ai_boundary = 1/2*[P_next(2)-P_prev(2), -P_next(1)+P_prev(1)] * Phi * thickness * stress;

                if k == 1  % Linear (k == 1) approximation with 2 points on an edge.
                    B_ai_volume = zeros(1,2);
                else
                   % Quadratic or higher approximation
                   % B_volume = ...
                   error(message('Case k > 1 has not been implemented yet'));
                end
                % Vector 1x2
                B_alphai = B_ai_boundary + B_ai_volume;
            else
                % Calculating mean value of vertices
                B_alphai = 1/area*thickness*1/2*([P_next(2)-P(2), -(P_next(1)-P(1))] - [P(2)-P_prev(2), -(P(1)-P_prev(1))] * Phi);
            end
        end
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Updating node2dof --> node2dof2
        % minh.nguyen@ikm.uni-hannover.de
        % 3/2/2017
        % Update node2dof2 --> this.node2dof
        % 11/04/2017
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       function dof = node2dof(this, setnodes, type)
            % There are 2 types: BOUNDARY type and ELEMENT type for 2 different
            % functions
            % ELEMENT type - nodes is an array containing node number i.e [2, 4, 6, 9, 10] is case
            % BOUNDARY type - nodes is the boundary array
            dim = this.method_.dim_;
            if strcmpi(type,'boundary')
                for i = 1 : length(setnodes)
                    dirichlet = setnodes{i};
                    dof = [];
                    if ~isempty(dirichlet)
                        if isrow(dirichlet)
                            dirichlet = dirichlet';
                        end
                        dirnode = unique(dirichlet);
                        %uXY = this.uFunc_{i}(this.mesh_.vertices(dirnode,1),this.mesh_.vertices(dirnode,2));
                        uXY = feval(this.uFunc_{i},this.mesh_.vertices(dirnode,1),this.mesh_.vertices(dirnode,2));
                        idxX = find(uXY(:,1) ~= -999999); % -999999 means these nodes in (x-direction in this case) are free movement
                        idxY = setdiff(find(uXY ~= 999999),idxX); % -999999 means these nodes in (y-direction in this case) are free movement
                        repDirnode = repmat(dirnode,1,dim); % Use for generating the D.O.F corresponding to node number
                        dof = [repDirnode(idxX)*dim-1;repDirnode(idxY)*dim];
                    end
                end            
            elseif strcmpi(type, 'element')
                dof = zeros(length(setnodes)*dim,1);
                dof(1:2:end) = setnodes*dim-1;
                dof(2:2:end) = setnodes*dim;
            end
       end
        
        function dof = node2dof_old(this, nodes, type)
            % There are 2 types: BOUNDARY type and ELEMENT type for 2 different
            % functions
            % ELEMENT type - nodes is an array containing node number i.e [2, 4, 6, 9, 10] is case
            % BOUNDARY type - nodes is the boundary array
            dof = [];
            dim = this.method_.dim_;
            if strcmpi(type,'boundary')                
                if ~iscell(nodes)
                    % -999999 la node do se free (Khong bi rang buoc)
                    if size(nodes,1) ~= 0
                        for i = 1 : size(nodes,1)
                            node = nodes(i,1);
                            if nodes(i,2) ~= -999999 % x-direction of the node is NOT free in INHOMOGENEOUS BOUNDARY CASE
                                dof(end+1) = node*dim - 1;
                            end
                            if nodes(i,3) ~= -999999 % x-direction of the node is NOT free in INHOMOGENEOUS BOUNDARY CASE
                                dof(end+1) = node*dim;
                            end
                        end
                    end
                else
                    for k = 1 : length(nodes)
                        nods = nodes{k};
                        % -999999 la node do se free (Khong bi rang buoc)
                        if size(nods,1) ~= 0
                            for i = 1 : size(nods,1)
                                node = nods(i,1);
                                if nods(i,2) ~= -999999 % x-direction of the node is NOT free in INHOMOGENEOUS BOUNDARY CASE
                                    dof(end+1) = node*dim - 1;
                                end
                                if nods(i,3) ~= -999999 % x-direction of the node is NOT free in INHOMOGENEOUS BOUNDARY CASE
                                    dof(end+1) = node*dim;
                                end
                            end
                        end
                    end
                end
            elseif strcmpi(type, 'element')
                for i = 1 : length(nodes)
                    dof(end+1) = nodes(i)*dim - 1;
                    dof(end+1) = nodes(i)*dim;
                end
            end
        end
        
        
        function U = applyDirichletBC(this, U)
           dim  = this.method_.dim_;
           vertices  = this.mesh_.vertices;
           uFunction = this.uFunc_;
           dirichlet   = this.dirichlet_;
           for i = 1 : length(dirichlet)
               pairNodes = dirichlet{i};
               if ~isempty(pairNodes)
                   if isrow(pairNodes) 
                       pairNodes = pairNodes'; 
                   end
                   nodes = unique(pairNodes);
                   %uXY = uFunction{i}(vertices(nodes,1),vertices(nodes,2));
                   uXY = feval(uFunction{i},vertices(nodes,1),vertices(nodes,2));
                   idxX = find(uXY(:,1) ~= -999999); % values -999999 means this node is free to move. We don't take account to the free nodes
                   idxY = setdiff(find(uXY ~= -999999),idxX);
                   uX   = uXY(idxX);
                   uY   = uXY(idxY);
                   nodeXY = repmat(nodes,1,dim); % Use for generating the D.O.F corresponding to node number
                   U([nodeXY(idxX)*dim-1;nodeXY(idxY)*dim]) = U([nodeXY(idxX)*dim-1;nodeXY(idxY)*dim]) + [uX;uY]; 
               end
           end
        end
        
        function F = applyNeumannBC(this, F)
           type = this.method_.typeF_;
           dim  = this.method_.dim_;
           thickness = this.material_.thickness_;
           vertices  = this.mesh_.vertices;
           tFunction = this.tFunc_;
           neumann   = this.neumann_;
           if strcmp(type,'traction')
               for i = 1 : length(neumann)
                   pairNodes = neumann{i};
                   Node1 = vertices(pairNodes(:,1),:);
                   Node2 = vertices(pairNodes(:,2),:);
                   area = (sqrt((Node2(:,1)-Node1(:,1)).^2 + (Node2(:,2)-Node1(:,2)).^2))/2 * thickness;
                   for j = 1 : 2
%                        txy = feval(tFunction{i}(vertices(pairNodes(:,j),1),vertices(pairNodes(:,j),2));
                       txy = feval(tFunction{i},vertices(pairNodes(:,j),1),vertices(pairNodes(:,j),2));
                       F([pairNodes(:,j)*dim-1;pairNodes(:,j)*dim]) = F([pairNodes(:,j)*dim-1;pairNodes(:,j)*dim]) + [area.*txy(:,1);area.*txy(:,2)]; 
                   end
               end
           end
        end
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % apply force which is either concentrated force or traction force for the node prescribed
        % Author: MINH T.V NGUYEN 
        % Email: minh.nguyen@ikm.uni-hannover.de
        % Change this function to a method of the class Virtual Element
        % From applyLoading --> applyNeumann
        % 11/04/2017
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function F = applyNeumann(this, F, loading)
            type = this.method_.typeF_;
            dim  = this.method_.dim_;
            thickness = this.material_.thickness_;
            vertices = this.mesh_.vertices;
            if ~iscell(loading)
                if strcmpi(type, 'concentrated')
                    for i = 1 : size(loading, 1)
                        F(loading(i,1)*dim-1) = loading(i,2);
                        F(loading(i,1)*dim)   = loading(i,3);
                    end
                elseif strcmpi(type, 'traction')
                   for i = 1 : size(loading,1)-1
                       Node1 = vertices(loading(i,1),:);
                       Node2 = vertices(loading(i+1,1),:);
                       w1 = (sqrt((Node2(1)-Node1(1))^2 + (Node2(2)-Node1(2))^2))/2 * thickness;
                       w2 = (sqrt((Node1(1)-Node2(1))^2 + (Node1(2)-Node2(2))^2))/2 * thickness;
                       T1 = w1 * loading(i,2:3)';
                       T2 = w2 * loading(i+1,2:3)';
                       F(loading(i,1)*dim-1:loading(i,1)*dim) = F(loading(i,1)*dim-1:loading(i,1)*dim) + T1;
                       F(loading(i+1,1)*dim-1:loading(i+1,1)*dim) = F(loading(i+1,1)*dim-1:loading(i+1,1)*dim) + T2;
                   end
                end
            else
                if strcmpi(type, 'concentrated')
                    for k = 1 : length(loading)
                        load = loading{k};
                        for i = 1 : size(load, 1)
                            F(load(i,1)*dim-1) = load(i,2);
                            F(load(i,1)*dim)   = load(i,3);
                        end
                    end
                elseif strcmpi(type, 'traction')
                    for k = 1 : length(loading)
                        load = loading{k};
                        for i = 1 : size(load,1)-1
                           Node1 = vertices(load(i,1),:);
                           Node2 = vertices(load(i+1,1),:);
                           w1 = (sqrt((Node2(1)-Node1(1))^2 + (Node2(2)-Node1(2))^2))/2 * thickness;
                           w2 = (sqrt((Node1(1)-Node2(1))^2 + (Node1(2)-Node2(2))^2))/2 * thickness;
                           T1 = w1 * load(i,2:3)';
                           T2 = w2 * load(i+1,2:3)';
                           F(load(i,1)*dim-1:load(i,1)*dim) = F(load(i,1)*dim-1:load(i,1)*dim) + T1;
                           F(load(i+1,1)*dim-1:load(i+1,1)*dim) = F(load(i+1,1)*dim-1:load(i+1,1)*dim) + T2;
                        end
                    end
                end
            end
        end
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Change this function to a method of the class Virtual Element
        % From applyBoundaryConditions --> applyDirichlet
        % 11/04/2017
        % apply boundary conditions for the node prescribed
        % Author: MINH T.V NGUYEN 
        % Email: minh.nguyen@ikm.uni-hannover.de
        % The code currently is applying for homogenous boundary which means the
        % boundary is either fixed or not fixed
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function U = applyDirichlet(this, U, boundary)
            dim  = this.method_.dim_;
            if ~iscell(boundary)
                for i = 1 : size(boundary, 1)
                    if boundary(i,2) ~= -999999 % truc x khong chuyen dong tu do (bi fixed o truc x)
                        U(boundary(i,1)*dim-1) = boundary(i,2);
                    end
                    if boundary(i,3) ~= -999999 % truc y khong chuyen dong tu do (bi fixed o truc y)
                        U(boundary(i,1)*dim) = boundary(i,3);
                    end
                end
            else
                for k = 1 : length(boundary)
                    bound = boundary{k};
                    for i = 1 : size(bound, 1)
                        if bound(i,2) ~= -999999 % truc x khong chuyen dong tu do (bi fixed o truc x)
                            U(bound(i,1)*dim-1) = bound(i,2);
                        end
                        if bound(i,3) ~= -999999 % truc y khong chuyen dong tu do (bi fixed o truc y)
                            U(bound(i,1)*dim) = bound(i,3);
                        end
                    end
                end
            end
        end
        
       %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       % Minh T.V Nguyen minh.nguyen@ikm.uni-hannover.de
       % Initial date: 24/04/2017
       % Integrate function getDisplacementOfElement to the
       % VirtualElementClass
       %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       function Ue = getDisplacementOfElement(this, U, element)
            dim = this.method_.dim_;
            Ue = zeros(dim*length(element),1);
            for i = 1 : length(element)
                Ue(i*dim-1:i*dim) = U(element(i)*dim-1:element(i)*dim);
            end
       end
    end 
end

