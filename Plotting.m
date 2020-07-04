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

classdef Plotting < handle
    %PLOTTING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(GetAccess = public, SetAccess = private)
        mesh_;
        dim_;
    end
    
    methods(Access = public)
        function this = Plotting(mesh, dim)
           disp('!!! PLOTTING !!!');
           this.mesh_ = mesh;
           this.dim_  = dim;
        end
        
        %
        % minh.nguyen@ikm.uni-hannover.de
        % Initial date: 03/09/2017
        %
        function plotQNode(this, crack, radius, JDomain, atTipNumber, qnode)
            nodes = this.mesh_.vertices;
            elements = this.mesh_.elements;
            figure;
            for i = 1 : length(JDomain)
               elementsJDomain{i} = elements{JDomain(i)};
            end
            for i = 1 : length(elements)
                elm = elements{i};
                X = zeros(1,length(elm));
                Y = zeros(1,length(elm));
                for n = 1 : length(elm)
                    X(n) = nodes(elm(n),1);
                    Y(n) = nodes(elm(n),2); 
                end
                patch(X, Y, double(qnode{i}),'LineStyle', 'None'); % convert logical array to numeric array
%                 hold on; 
%                 X(end+1) = X(1);
%                 Y(end+1) = Y(1);
%                 plot(X,Y,'b');
                colormap jet;
                hold on; 
            end
            hold on
            % plot the circle
            theta = -pi:0.1:pi;
            if exist('atTipNumber','var')
                xTip = crack.crack_tip_(atTipNumber,:); 
            else
                xTip = crack.crack_tip_(1,:);
            end
            xCr = crack.seam_;
            xo = xTip(1) + radius*cos(theta) ;
            yo = xTip(2) + radius*sin(theta) ;
             plot([xo xo(1)],[yo yo(1)],'k-');
             this.plot_mesh(nodes,elements,'Poly','b-');
             this.plot_mesh(nodes,elementsJDomain,'Poly','r-');
             cr = plot(xCr(:,1),xCr(:,2),'k-');
            set(cr,'LineWidth',2.0);
        end
        
        % minh.nguyen@ikm.uni-hannover.de
        % Added on 20/Jan/2017
        function plotDomainIntegrationAtCrackTip(this, crack, radius, JDomain, atTipNumber, typeElm)
            nodes = this.mesh_.vertices;
            elements = this.mesh_.elements;
            if exist('typeElm','var')
                elmType = typeElm; 
            else
                elmType = 'T3';
            end
            for i = 1 : length(JDomain)
               elementsJDomain{i} = elements{JDomain(i)};
            end
           % plot
            figure;
            hold on
            % plot the circle
            theta = -pi:0.1:pi;
            if exist('atTipNumber','var')
                xTip = crack.crack_tip_(atTipNumber,:); 
            else
                xTip = crack.crack_tip_(1,:);
            end
            xCr = crack.seam_;
            xo = xTip(1) + radius*cos(theta) ;
            yo = xTip(2) + radius*sin(theta) ;
            plot([xo xo(1)],[yo yo(1)],'k-');
            this.plot_mesh(nodes,elements,elmType,'b-');
            this.plot_mesh(nodes,elementsJDomain,elmType,'r-')
            cr = plot(xCr(:,1),xCr(:,2),'k-');
            set(cr,'LineWidth',2.0);
            axis off;
            % ----
        end
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % PLOT FIELD
        % minh.nguyen@ikm.uni-hannover.de
        % Initial date: 04.08.2018
        % Original code: Nguyen Vinh Phu's code
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function plot_field(this, field)
                X = this.mesh_.vertices;
                connect = this.mesh_.elements;
                holdState=ishold;
                if (size(field) == length(connect))
                    elementalField = 1;
                else
                    elementalField = 0;
                end
                % fill X if needed
                if (size(X,2) < 3)
                    for c = size(X,2)+1:3
                        X(:,c) = zeros(size(X,1),1);
                    end
                end
                for e=1:length(connect)
                    xpt = zeros(length(connect{e})+1,1);
                    ypt = zeros(length(connect{e})+1,1);
                    zpt = zeros(length(connect{e})+1,1);
                    for n = 1 : length(connect{e})
                        xpt(n) = X(connect{e}(n),1);
                        ypt(n) = X(connect{e}(n),2);
                        zpt(n) = X(connect{e}(n),3);
                    end
                    xpt(end) = xpt(1);
                    ypt(end) = ypt(1);
                    zpt(end) = zpt(1);
                    if (elementalField)
                        fpt = [field{e},field{e}(1)];
                    else
                        fpt = field([connect{e},connect{e}(1)]);
                    end
                    fill3(xpt,ypt,zpt,fpt);
                    msh = plot3(xpt,ypt,zpt,'w-');
                    set(msh,'LineWidth',1.0);
                end
                shading interp
                axis equal
                if (~holdState )
                    hold off
                end
        end

        % minh.nguyen@ikm.uni-hannover.de
        % Added on 20/Jan/2017
           function plot_mesh(this, X,connect,elem_type,se)

            % function plot_mesh(X,connect,elem_type,linespec)
            % 
            % plots a nodal mesh and an associated connectivity.  X is
            % teh nodal coordinates, connect is the connectivity, and
            % elem_type is either 'L2', 'L3', 'T3', 'T6', 'Q4', or 'Q9' 
            % depending on the element topology.

            if ( nargin < 5 )
               se='w-';
            end

            holdState=ishold;
            hold on

            % fill X if needed
            if (size(X,2) < 3)
               for c=size(X,2)+1:3
                  X(:,c)=[zeros(size(X,1),1)];
               end
            end

            for e=1:length(connect)
%                  ord = 1:length(connect{e});
%                  ord(end+1) = 1;
%                if ( strcmp(elem_type,'Q9') )       % 9-node quad element
%                   ord=[1,5,2,6,3,7,4,8,1];
%                elseif ( strcmp(elem_type,'Q8') )  % 8-node quad element
%                   ord=[1,5,2,6,3,7,4,8,1];
%                elseif ( strcmp(elem_type,'T3') )  % 3-node triangle element
%                   ord=[1,2,3,1];
%                elseif ( strcmp(elem_type,'T6') )  % 6-node triangle element
%                   ord=[1,4,2,5,3,6,1];
%                elseif ( strcmp(elem_type,'Q4') )  % 4-node quadrilateral element
%                   ord=[1,2,3,4,1];
%                elseif ( strcmp(elem_type,'L2') )  % 2-node line element
%                   ord=[1,2];   
%                elseif ( strcmp(elem_type,'L3') )  % 3-node line element
%                   ord=[1,3,2];   
%                elseif ( strcmp(elem_type,'H4') )  % 4-node tet element
%                   ord=[1,2,4,1,3,4,2,3];   
%                elseif ( strcmp(elem_type,'B8') )  % 8-node brick element
%                   ord=[1,5,6,2,3,7,8,4,1,2,3,4,8,5,6,7];   
%                end

               xpt = zeros(length(connect{e})+1,1);
               ypt = zeros(length(connect{e})+1,1);
               zpt = zeros(length(connect{e})+1,1);
               for n=1:length(connect{e})
                    xpt(n)=X(connect{e}(n),1);
                    ypt(n)=X(connect{e}(n),2);
                    zpt(n)=X(connect{e}(n),3);
               end
               xpt(end) = xpt(1);
               ypt(end) = ypt(1);
               zpt(end) = zpt(1);
               msh = plot3(xpt,ypt,zpt,se);
               set(msh,'LineWidth',1.0);
            end

            rotate3d on
            axis equal

            if ( ~holdState )
              hold off
            end
        end
            
        function plotUndeformation(this)
            vertices = this.mesh_.vertices;
            elements = this.mesh_.elements;
            for e = 1 : length(elements)
                elm = elements{e};
                vert = vertices(elm,:);
                Ndof = length(elm);
                [area, centroid, diameter] = getElementInfo(vert, Ndof);
                
                X = zeros(1, length(elm));
                Y = zeros(1, length(elm));
                for n = 1 : length(elm)
                    X(n) = vertices(elm(n),1);
                    Y(n) = vertices(elm(n),2);
                end
                X(n+1) = X(1);
                Y(n+1) = Y(1);
                plot(X,Y,'b');
                text(centroid(1),centroid(2),num2str(e),'FontWeight','bold');
                hold on;
            end
        end
        
        function plotDeformation(this, displacements)
%  
%         X, Y:   Coordinates of un-deform structure.
%         XD, YD: Coordinates of deform structure.
%             
            width = 13;
            height = 5;
            alw = 0.75;
            fsz = 11;
            lw = 1.0;
            msz = 8;
            
            figure('Name', 'Plotting Deformation');
            pos = get(gcf,'Position');
            set(gcf, 'Position', [pos(1) pos(2) width*100 height*100]);
%             set(gca,'position',[pos(1) pos(2) width*100 height*100],'units','normalized')
            set(gca, 'FontSize', fsz, 'LineWidth', alw);
            box off;
            vertices = this.mesh_.vertices;
            elements = this.mesh_.elements;
            for e = 1 : length(elements)
% Array         elm = elements(e,:);
                elm = elements{e};
                % Preallocate the arrays of X and Y, XD and YD
                X  = zeros(1, length(elm));
                Y  = zeros(1, length(elm));
                XD = zeros(1, length(elm));
                YD = zeros(1, length(elm));
                fac = 1; % de cho nhin ro hon su bien dang
                   for n = 1 : length(elm)
                      X(n)  = vertices(elm(n),1);
                      XD(n) = X(n) + fac*displacements(elm(n)*this.dim_-1);
                      Y(n)  = vertices(elm(n),2);
                      YD(n) = Y(n) + fac*displacements(elm(n)*this.dim_);
                   end
                   X(n+1)  = X(1);
                   Y(n+1)  = Y(1);
                   XD(n+1) = XD(1);
                   YD(n+1) = YD(1);
%                    plot(X,Y,'b');
                    
                    plot(XD,YD,'b', 'LineWidth',lw,'MarkerSize',msz);
                   axis off;
                   hold on;
                   
            end
        end
        
        function plotStress(this, Stress_elms, direction)
           figure('Name', 'Stress visualization');
           axis equal; colorbar('vert');
           vertices = this.mesh_.vertices;
           elements = this.mesh_.elements;
           for e = 1 : length(elements)
               elm = elements{e};
               X = zeros(1, length(elm));
               Y = zeros(1, length(elm));
               for n = 1 : length(elm)
                  X(n) = vertices(elm(n),1);
                  Y(n) = vertices(elm(n),2);
                  stress = Stress_elms{e}(:,n);
               end
               patch(X, Y, stress(direction), 'LineStyle', 'None');
               colormap jet;
               hold on;
           end
        end
        
        %
        % Initial date: 25.08.2017
        % Author: minh.nguyen@ikm.uni-hannover.de
        %
        function plotStressInDeformation(this, Stress_elms, direction, u2D)
           figure('Name', 'Stress visualization');
           width = 7;
           height = 16;
%            axis equal;
           pos = get(gcf,'Position');
%            set(hc,'xaxisloc','top');
%            colorbar('vert');
           set(gcf, 'Position', [pos(1) pos(2) width*100 height*100]);
%            set(gca,'position',[0 0 0.87 1],'units','normalized');
           box off;
           vertices = this.mesh_.vertices;
           elements = this.mesh_.elements;
           for e = 1 : length(elements)
               elm = elements{e};
               X = zeros(1, length(elm));
               Y = zeros(1, length(elm));
               for n = 1 : length(elm)
                  X(n) = vertices(elm(n),1) + u2D(elm(n),1);
                  Y(n) = vertices(elm(n),2) + u2D(elm(n),2);
                  stress = Stress_elms{e}(:,n);
               end
               patch(X, Y, stress(direction), 'LineStyle', 'None');
               colormap jet;
               axis off;
               hold on;
           end
        end
        
    end
    
end
