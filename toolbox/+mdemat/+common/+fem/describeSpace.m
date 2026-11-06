function description = describeSpace(object)
%DESCRIBESPACE Describe built-in representations using a fixed data grammar.
%   D = describeSpace(S) captures the current mesh, geometry, basis, element,
%   quadrature and composite-space parameters. D contains only values and
%   explicit versioned tags, never objects, callbacks or class constructors.
%   Mesh/SL/observer caches are omitted. Physical quadrature Nodes/Weights,
%   basis normalization/order, mass/design matrices and sum/scaling weights
%   are retained. Unsupported user classes fail before a record is written.
%
%   See also mdemat.common.fem.restoreSpace.
type = class(object);
data = struct();
switch type
    case 'mdemat.data.mesh.UniformGrid'
        tag = 'mesh.uniform';
        data.bounds = reshape(cell2mat(cellfun(@(x) [x(1); x(end)], ...
            object.Vertices, UniformOutput=false)), 1, []);
        data.resolution = object.Resolution;
        data.boundary = object.Boundary;
    case 'mdemat.data.mesh.Grid'
        tag = 'mesh.grid'; data.vertices = object.Vertices; data.boundary = object.Boundary;
    case 'mdemat.data.mesh.Triangulation'
        tag = 'mesh.simplex';
        data.vertices = object.Vertices; data.elements = object.Elements;
        data.boundary = object.Boundary;
        data.vertexPeriodicity = object.VertexToVertexTables;
        data.facePeriodicity = object.FaceToFaceTable;
    case 'mdemat.data.mesh.Graph'
        tag = 'graph'; data.vertices = object.Vertices; data.edges = object.Edges;
    case 'mdemat.data.mesh.SeparableGraph'
        tag = 'graph.separable'; data.vertices = object.Vertices; data.edges = object.Edges;
    case 'mdemat.data.geometry.Orthotope'
        tag = 'geometry.orthotope'; data.bounds = object.BBox;
    case {'mdemat.data.geometry.Sphere', 'mdemat.data.geometry.Ball'}
        if endsWith(type, '.Sphere'), tag = 'geometry.sphere'; else, tag = 'geometry.ball'; end
        data.center = object.Center; data.radius = object.Radius;
    case 'mdemat.data.geometry.Simplex'
        tag = 'geometry.simplex'; data.vertices = object.Vertices;
    case 'mdemat.data.geometry.Polytope'
        tag = 'geometry.polytope'; data.vertices = object.Vertices; data.faces = object.Faces;
    case {'mdemat.common.fem.MeshSpace', 'mdemat.common.fem.FiniteElementSpace'}
        data.mesh = mdemat.common.fem.describeSpace(object.Mesh);
        if endsWith(type, '.FiniteElementSpace')
            tag = 'space.fem'; data.element = mdemat.common.fem.describeSpace(object.Element);
        else
            tag = 'space.mesh';
        end
    case {'mdemat.common.fem.SpectralSpace', 'mdemat.common.fem.AffineSpace', ...
            'mdemat.common.fem.ScaledAffineSpace'}
        data.element = mdemat.common.fem.describeSpace(object.Element);
        switch type
            case 'mdemat.common.fem.SpectralSpace', tag = 'space.spectral';
            case 'mdemat.common.fem.AffineSpace', tag = 'space.affine';
            otherwise, tag = 'space.scaled'; data.scales = object.scales;
        end
    case 'mdemat.common.fem.SumSpace'
        tag = 'space.sum';
        data.lhs = describePart(object.Lhs); data.rhs = describePart(object.Rhs);
        data.lhsWeight = object.LhsWeight; data.rhsWeight = object.RhsWeight;
    case {'mdemat.common.fem.Element', 'mdemat.common.fem.L2Element', ...
            'mdemat.common.fem.H1Element', 'mdemat.common.fem.BH1Element', ...
            'mdemat.common.fem.L2OrthotopeElement', 'mdemat.common.fem.L2SphereElement', ...
            'mdemat.common.fem.BH1OrthotopeElement', 'mdemat.common.fem.BH1SLOrthotopeElement'}
        names = {'Element', 'L2Element', 'H1Element', 'BH1Element', 'L2OrthotopeElement', ...
            'L2SphereElement', 'BH1OrthotopeElement', 'BH1SLOrthotopeElement'};
        tags = {'base', 'l2', 'h1', 'broken', 'l2-orthotope', 'l2-sphere', ...
            'broken-orthotope', 'sl-orthotope'};
        tag = ['element.', tags{strcmp(type, strcat('mdemat.common.fem.', names))}];
        data.geometry = mdemat.common.fem.describeSpace(object.Geometry);
        data.approximator = mdemat.common.fem.describeSpace(object.Approximator);
        if ~strcmp(tag, 'element.base')
            data.volume = mdemat.common.fem.describeSpace(object.Volume);
        end
        if isa(object, 'mdemat.common.fem.BH1Element'), data.flux = describeArray(object.Flux); end
        if strcmp(tag, 'element.l2-sphere'), data.reduction = object.DimReduction; end
    case {'mdemat.common.approximation.LinearApproximator', ...
            'mdemat.common.approximation.ModalApproximator', ...
            'mdemat.common.approximation.NodalApproximator'}
        if endsWith(type, '.ModalApproximator')
            tag = 'approximation.modal';
        elseif endsWith(type, '.NodalApproximator')
            tag = 'approximation.nodal';
        else
            tag = 'approximation.linear';
        end
        data.basis = mdemat.common.fem.describeSpace(object.Basis);
        data.mass = object.Mass; data.design = object.Design;
    case 'mdemat.common.quadrature.Integrator'
        tag = 'quadrature.integrator';
        data.rule = mdemat.common.fem.describeSpace(object.Rule);
        data.nodes = object.Nodes; data.weights = object.Weights;
    case 'mdemat.common.quadrature.MultiLevelIntegrator'
        tag = 'quadrature.multilevel'; data.integrators = describeArray(object.Integrators);
    case {'mdemat.common.quadrature.GaussLegendreRule', 'mdemat.common.quadrature.GaussLobattoRule', ...
            'mdemat.common.quadrature.ClosedNewtonCotesRule', ...
            'mdemat.common.quadrature.PeriodicTrapezoidalRule', ...
            'mdemat.common.quadrature.GaussHermiteRule', 'mdemat.common.quadrature.LebedevRule', ...
            'mdemat.common.quadrature.GaussTrapezoidalRule', ...
            'mdemat.common.quadrature.MonteCarloOrthotopeRule', ...
            'mdemat.common.quadrature.MonteCarloSphereRule'}
        tag = 'quadrature.rule';
        data.kind = erase(type, 'mdemat.common.quadrature.');
        data.dimensions = object.NDims;
        data.coordinates = ''; data.seed = [];
        if isa(object, 'mdemat.common.quadrature.SphereRule'), data.coordinates = object.Coord; end
        if isprop(object, 'Seed'), data.seed = object.Seed; end
    case {'mdemat.common.approximation.CompiledFunction', ...
            'mdemat.common.spectral.CompiledOrthogonalBasis', ...
            'mdemat.common.interpolation.CompiledInterpolationBasis'}
        if endsWith(type, '.CompiledOrthogonalBasis')
            tag = 'basis.compiled-orthogonal';
        elseif endsWith(type, '.CompiledInterpolationBasis')
            tag = 'basis.compiled-interpolation';
        else
            tag = 'basis.compiled';
        end
        data.dimensions = object.NDims; data.components = object.NCodims;
        data.resource = object.Metadata.basisName;
        data.definition = struct('schemaVersion', object.Metadata.schemaVersion, ...
            'derivativeOrders', object.Metadata.derivativeOrders, ...
            'functionText', {object.Metadata.functionText});
    case {'mdemat.common.spectral.LegendreBasis', 'mdemat.common.spectral.HermiteBasis', ...
            'mdemat.common.spectral.ChebyshevBasis', 'mdemat.common.spectral.FourierBasis', ...
            'mdemat.common.spectral.SphericalHarmonicBasis'}
        tag = 'basis.recurrence'; data.kind = erase(type, 'mdemat.common.spectral.');
        data.degree = object.MaxDegree; data.normalized = object.IsNormalized;
        data.lower = []; data.upper = [];
        if strcmp(data.kind, 'FourierBasis'), data.lower = object.Lower; data.upper = object.Upper; end
    case 'mdemat.common.interpolation.CardinalKernelBasis'
        tag = 'basis.cardinal'; data.nodes = object.Nodes; data.weights = object.Weights;
    case 'mdemat.common.interpolation.RbfKernelBasis'
        tag = 'basis.rbf'; data.nodes = object.Nodes;
        data.kernel = object.Kernel; data.epsilon = object.Epsilon;
    case 'mdemat.common.approximation.SeparableFunction'
        tag = 'function.separable'; data.factors = describeArray(object.Factors);
        data.pattern = '';
        if ~isempty(object.Pattern)
            if isa(object.Pattern.indexer, 'mdemat.data.index.L1MultiIndexer')
                data.pattern = 'P';
            elseif isa(object.Pattern.indexer, 'mdemat.data.index.LxMultiIndexer')
                data.pattern = 'Q';
            else
                unsupported(object);
            end
        end
    case 'mdemat.data.function.Constant'
        tag = 'function.constant'; data.dimensions = object.NDims; data.value = object.c;
    case 'mdemat.data.function.Linear'
        tag = 'function.linear'; data.matrix = object.A; data.shift = object.b;
    case 'mdemat.data.function.Maxwellian'
        tag = 'function.maxwellian'; data.dimensions = object.NDims;
        data.temperature = object.Temperature; data.center = object.Center;
    case {'mdemat.data.function.Product', 'mdemat.data.function.Composition'}
        if endsWith(type, '.Product'), tag = 'function.product'; else, tag = 'function.composition'; end
        data.lhs = mdemat.common.fem.describeSpace(object.Lhs);
        data.rhs = mdemat.common.fem.describeSpace(object.Rhs);
    otherwise
        unsupported(object);
end
description = struct('schemaVersion', 1, 'type', tag, 'data', data);
end

function descriptions = describeArray(objects)
descriptions = cell(1, numel(objects));
for index = 1:numel(objects)
    if iscell(objects), object = objects{index}; else, object = objects(index); end
    descriptions{index} = mdemat.common.fem.describeSpace(object);
end
end

function description = describePart(object)
if isnumeric(object) && isempty(object)
    description = struct('schemaVersion', 1, 'type', 'space.empty', 'data', struct());
else
    description = mdemat.common.fem.describeSpace(object);
end
end

function unsupported(object)
error('mdemat:persistence:UnsupportedRepresentation', ...
    'No data-only representation is registered for %s.', class(object));
end
