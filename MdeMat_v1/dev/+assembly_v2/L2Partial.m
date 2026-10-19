classdef L2Partial < fem.assembly.Operator
    % L2PARTIAL 

    properties
        volume % Volume assembly
        flux % Flux assembly on interior elements
    end

    methods
        function obj = Partial(volume, interior, boundary)
            obj.volume = volume;
            obj.interior = interior;
            obj.boundary = boundary;
        end
    
        function A = assemble(obj, d, c, fluxType)
            c = obj.volume.assembleConstant(c);
            T1 = obj.volume.assembleVolumePartial(d, c);
            switch fluxType
                case 1
                    T2 = obj.flux.assembleLeftFlux(d, c);
                case 2
                    T2 = obj.flux.assembleRightFlux(d, c);
                case 3
                    T2 = obj.flux.assembleCentralFlux(d, c);
            end
            n = obj.volume.nGlobalDofs;
            T = [T1; T2];
            A = obj.sparseFromTriplets(T, n, n);
        end

        function b = assembleBc(obj, d, c, a, f)
            c = obj.volume.assembleConstant(c);
            if a > 0
                T = obj.assembleLeftTrace(d, c, f);
            else
                T = obj.assembleRightTrace(d, c, f);
            end
            n = obj.volume.nGlobalDofs;
            b = obj.sparseFromTriplets(T, n, 1);
        end
    end
end