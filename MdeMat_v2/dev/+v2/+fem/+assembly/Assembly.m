classdef Assembly < handle
    % ASSEMBLY 
    
    properties
        feSp % Finite element space
        feOp % Finite element operator
    end
    
    methods
        function obj = Assembly(feSp, feOp)
            % ASSEMBLY Constructor for Assembly.
            %
            %   obj = Assembly(feSp, feOp) creates an assembly object with
            %   the specified finite element space and finite element
            %   operator.
            %
            % Inputs:
            %   feSp - FiniteElementSpace object
            %
            % Outputs:
            %   obj - Constructed Assembly object

            obj.feSp = feSp;
            obj.feOp = feOp;
        end
    end
end

