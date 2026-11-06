classdef OrthotopeRule < mdemat.common.quadrature.IntegrationRule
    % ORTHOTOPERULE Integration rule over the orthotope.
    %
    %   OrthotopeRule creates multidimensional integration rules by taking
    %   tensor products of one-dimensional rules. This approach is
    %   efficient for integrating functions that can be separated into
    %   products of univariate functions or for rectangular integration
    %   domains.
    %
    %   The orthotope rule constructs multidimensional quadrature by
    %   taking tensor products of univariate rules. This is most efficient
    %   for functions that exhibit separability or low-rank structure.
    %
    % See also:
    %   mdemat.common.quadrature.IntegrationRule
    
    methods
        function [X, w] = generateOnBoundary(obj, geometry, np)
            % GENERATEONBOUNDARY Integrate on one axis-aligned face.
            %   geometry has exactly one collapsed coordinate. That fixed
            %   coordinate carries unit point measure; the remaining axes
            %   carry their physical length measures. This differs from
            %   generate(), whose ambient-dimensional volume is zero.
            arguments
                obj mdemat.common.quadrature.OrthotopeRule
                geometry mdemat.data.geometry.Orthotope
                np {mustBeVector, mustBePositive, mustBeInteger}
            end
            fixed = geometry.Lower == geometry.Upper;
            assert(sum(fixed) == 1, 'MdeMat:Quadrature:InvalidFace', ...
                'A face must have exactly one collapsed coordinate.');
            if isscalar(np)
                np = repmat(np, 1, geometry.NDims);
            end
            assert(numel(np) == geometry.NDims, 'MdeMat:Quadrature:InvalidPointCount', ...
                'Specify one point count per ambient coordinate.');
            np(fixed) = 1;
            lower = geometry.Lower;
            upper = geometry.Upper;
            lower(fixed) = -0.5;
            upper(fixed) = 0.5;
            bounds = [lower; upper];
            auxiliary = mdemat.data.geometry.Orthotope(bounds(:));
            [X, w] = obj.generate(auxiliary, np);
            X(fixed, :) = geometry.Lower(fixed);
        end
    end

    methods (Access = protected)
        function [X, w] = generateImpl(obj, geometry, np)
            % GENERATEIMPL Generate integration nodes and weights.
            %
            %   [X, w] = generateImpl(obj, geometry, np) generates
            %   @a np integration points @a X and weights @a
            %   w for the geometric domain @a geometry.
            
            arguments
                obj mdemat.common.quadrature.OrthotopeRule
                geometry mdemat.data.geometry.Orthotope
                np {mustBePositive, mustBeInteger}
            end
            
            nd = geometry.NDims;
            
            if isscalar(np), np = repmat(np, 1, nd); end
            
            assert(all(length(np) == nd, 'all'), 'mdemat:common:quadrature:OrthotopeRule:InvalidInput', ...
                '%s', 'Length of np must equal the number of dimensions.');
            
            X = cell(1, nd);
            w = cell(1, nd);
            for i = 1:nd
                [X{i}, w{i}] = obj.generate1d(np(i), ...
                    lower = geometry.Lower(i), upper = geometry.Upper(i));
            end
            [X{1:nd}] = ndgrid(X{:});
            X = reshape(cat(nd+1, X{:}), [], nd).';
            
            [w{1:nd}] = ndgrid(w{:});
            w = reshape(prod(cat(nd+1, w{:}), nd+1), 1, []);
        end
    end
    
    methods (Abstract, Access = protected)
        % GENERATE1D Generate 1D integration nodes and weights.
        [x, w] = generate1d(obj, np, options)
    end
end
