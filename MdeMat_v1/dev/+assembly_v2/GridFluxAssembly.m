classdef GridFluxAssembly < fem.assembly.GridAssembly
    
    properties
        interior
        boundary
    end

    methods        
        function T = assembleLeftFlux(obj, d, c)
            % ASSEMBLELEFTFLUX Assemble left-biased flux terms.
            %
            %   T = assembleLeftFlux(obj, d, c) combines interior and
            %   boundary contributions for left-biased flux terms along the
            %   specified dimension.
            %
            % Inputs:
            %   obj - The EulerianGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Combined triplet matrix for left flux terms

            T1 = obj.interior.assembleInnerLeftFlux(d, c);
            T2 = obj.interior.assembleOuterLeftFlux(d, c);
            T3 = obj.boundary.assembleInnerLeftFlux(d, c);
            T4 = obj.boundary.assembleOuterLeftFlux(d, c);
            T = [T1; T2; T3; T4];
        end

        function T = assembleRightFlux(obj, d, c)
            % ASSEMBLYIGHTFLUX Assemble right-biased flux terms.
            %
            %   T = assembleRightFlux(obj, d, c) combines interior and
            %   boundary contributions for right-biased flux terms along
            %   the specified dimension.
            %
            % Inputs:
            %   obj - The EulerianGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Combined triplet matrix for right flux terms

            T1 = obj.interior.assembleInnerRightFlux(d, c);
            T2 = obj.interior.assembleOuterRightFlux(d, c);
            T3 = obj.boundary.assembleInnerRightFlux(d, c);
            T4 = obj.boundary.assembleOuterRightFlux(d, c);
            T = [T1; T2; T3; T4];
        end

        function T = assembleCentralFlux(obj, d, c)
            % ASSEMBLECENTRALFLUX Assemble central flux terms.
            %
            %   T = assembleCentralFlux(obj, d, c) combines left and
            %   right flux contributions with equal weighting to create
            %   central difference-type operators.
            %
            % Inputs:
            %   obj - The EulerianGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Combined triplet matrix for central flux terms

            T1 = obj.assembleLeftFlux(d, c / 2);
            T2 = obj.assembleRightFlux(d, c / 2);
            T = [T1; T2];
        end
    end
end

