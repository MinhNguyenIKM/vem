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

classdef Method < handle
    %METHOD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private, GetAccess = public)
        name_;
        dim_;
        k_;
        typeF_;
    end
    
    methods
        function this = Method(name, dimension, degree, forceType)
            this.name_ = name;
            this.dim_ = dimension;
            this.k_ = degree;
            if exist('forceType','var')
                this.typeF_ = forceType;
            end
        end
    end
    
end

