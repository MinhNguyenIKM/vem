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

classdef Material < handle
    %MATERIAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = private)
        E_;  % Young Modulus
        mu_; % Shear Modulus
        nu_; % Poisson ratio
        kappa_;
        material_;
        type_;
        thickness_; % Thickness
    end
    
    methods
        function this = Material(E, nu, type, h)
            % type : 'plane_strain' or 'plane_stress'
            this.E_  = E;
            this.nu_ = nu;
            this.mu_ = E/(2*(1+nu));
            this.thickness_ = h;
            this.type_ = type;
            if strcmp(type, 'plane_strain')
                this.material_ = E/((1+nu)*(1-2*nu))*[1-nu    nu      0;   ...
                                                        nu   1-nu    0;  ...
                                                         0   0   (1-2*nu)/2]; 
                this.kappa_ = 3 - 4*nu;
            elseif strcmp(type, 'plane_stress')
                this.material_ = E /(1 - nu^2) * [1       nu      0; ... 
                                                  nu      1       0;  ...
                                                  0       0      (1-nu)/2]; ...
                this.kappa_ = (3-nu)/(1+nu);
            else
                error('Plane Strain or Plane Stress only');
            end
        end
        
        function M = getMaterial(this)
            M = this.material_;
        end
        
        function E = getYoung(this)
            E = this.E_; 
        end
        
        function mu = getShear(this)
            mu = this.mu_;
        end
        
        function nu = getPoisson(this)
        	nu = this.nu_;
        end
        
        function kappa = getKappa(this)
            kappa = this.kappa_; 
        end
        
        function h = getThickness(this)
            h = this.thickness_;
        end
        
        function E_star = getEStar(this)
            if strcmp(this.type_, 'plane_strain')
                E_star = this.E_ / (1 - this.nu_*this.nu_); 
            elseif strcmp(this.type_, 'plane_stress')
                E_star = this.E_;
            end
        end
    end
    
end

