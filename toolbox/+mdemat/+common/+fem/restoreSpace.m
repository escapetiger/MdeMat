function object = restoreSpace(description)
%RESTORESPACE Reconstruct a built-in representation from a fixed value schema.
%   Object constructors are selected by explicit registered tags. Unknown
%   tags/fields/versions fail; saved class names, handles and expressions are
%   never executed. Quadrature sample coordinates/physical weights are kept;
%   basis values/derivatives and mesh transformation caches are rebuilt.
arguments
    description (1, 1) struct
end
try
    if ~valueTree(description), invalid('Representation descriptions must contain plain values.'); end
    fields(description, {'schemaVersion', 'type', 'data'});
    if ~isequal(description.schemaVersion, 1) || ...
            ~isstruct(description.data) || ~isscalar(description.data) || ...
            ~(ischar(description.type) && isrow(description.type))
        invalid('Invalid representation schema.');
    end
    data = description.data;
    switch description.type
        case 'mesh.uniform'
            fields(data, {'bounds', 'resolution', 'boundary'});
            object = mdemat.data.mesh.UniformGrid(data.bounds, data.resolution);
            object.Boundary = data.boundary;
        case 'mesh.grid'
            fields(data, {'vertices', 'boundary'});
            object = mdemat.data.mesh.Grid(data.vertices); object.Boundary = data.boundary;
        case 'mesh.simplex'
            fields(data, {'vertices', 'elements', 'boundary', 'vertexPeriodicity', 'facePeriodicity'});
            object = mdemat.data.mesh.Triangulation(data.vertices, data.elements);
            object.Boundary = data.boundary;
            object.VertexToVertexTables = data.vertexPeriodicity;
            object.FaceToFaceTable = data.facePeriodicity;
        case {'graph', 'graph.separable'}
            fields(data, {'vertices', 'edges'});
            if strcmp(description.type, 'graph')
                object = mdemat.data.mesh.Graph(data.vertices, data.edges);
            else
                object = mdemat.data.mesh.SeparableGraph(data.vertices, data.edges);
            end
        case 'geometry.orthotope'
            fields(data, {'bounds'}); object = mdemat.data.geometry.Orthotope(data.bounds);
        case {'geometry.sphere', 'geometry.ball'}
            fields(data, {'center', 'radius'});
            if strcmp(description.type, 'geometry.sphere')
                object = mdemat.data.geometry.Sphere(data.center, data.radius);
            else
                object = mdemat.data.geometry.Ball(data.center, data.radius);
            end
        case 'geometry.simplex'
            fields(data, {'vertices'}); object = mdemat.data.geometry.Simplex(data.vertices);
        case 'geometry.polytope'
            fields(data, {'vertices', 'faces'});
            object = mdemat.data.geometry.Polytope(data.vertices, data.faces);
        case 'space.empty'
            fields(data, {}); object = [];
        case 'space.mesh'
            fields(data, {'mesh'});
            object = mdemat.common.fem.MeshSpace(restore(data.mesh));
        case 'space.fem'
            fields(data, {'mesh', 'element'});
            object = mdemat.common.fem.FiniteElementSpace(restore(data.mesh), restore(data.element));
        case {'space.spectral', 'space.affine', 'space.scaled'}
            expected = {'element'};
            if strcmp(description.type, 'space.scaled'), expected{end+1} = 'scales'; end
            fields(data, expected);
            element = restore(data.element);
            switch description.type
                case 'space.spectral', object = mdemat.common.fem.SpectralSpace(element);
                case 'space.affine', object = mdemat.common.fem.AffineSpace(element);
                otherwise, object = mdemat.common.fem.ScaledAffineSpace(element, data.scales);
            end
        case 'space.sum'
            fields(data, {'lhs', 'rhs', 'lhsWeight', 'rhsWeight'});
            object = mdemat.common.fem.SumSpace(restore(data.lhs), restore(data.rhs), ...
                lhsWeight=data.lhsWeight, rhsWeight=data.rhsWeight);
        case {'element.base', 'element.l2', 'element.h1', 'element.broken', ...
                'element.l2-orthotope', 'element.l2-sphere', ...
                'element.broken-orthotope', 'element.sl-orthotope'}
            object = restoreElement(description.type, data);
        case {'approximation.linear', 'approximation.modal', 'approximation.nodal'}
            fields(data, {'basis', 'mass', 'design'});
            basis = restore(data.basis);
            switch description.type
                case 'approximation.linear', object = mdemat.common.approximation.LinearApproximator(basis);
                case 'approximation.modal', object = mdemat.common.approximation.ModalApproximator(basis);
                otherwise, object = mdemat.common.approximation.NodalApproximator(basis);
            end
            object.Mass = data.mass; object.Design = data.design;
        case 'quadrature.integrator'
            fields(data, {'rule', 'nodes', 'weights'});
            if ~isnumeric(data.nodes) || ~isnumeric(data.weights) || ...
                    ~ismatrix(data.nodes) || numel(data.weights) ~= size(data.nodes, 2)
                invalid('Quadrature coordinates and weights have inconsistent axes.');
            end
            object = mdemat.common.quadrature.Integrator(restore(data.rule));
            object.Nodes = data.nodes; object.Weights = data.weights;
        case 'quadrature.multilevel'
            fields(data, {'integrators'});
            parts = restoreArray(data.integrators);
            object = mdemat.common.quadrature.MultiLevelIntegrator(parts{1}.Rule, nLevels=numel(parts));
            object.Integrators = [parts{:}];
        case 'quadrature.rule'
            fields(data, {'kind', 'dimensions', 'coordinates', 'seed'});
            object = restoreRule(data);
        case {'basis.compiled', 'basis.compiled-orthogonal', 'basis.compiled-interpolation'}
            fields(data, {'dimensions', 'components', 'resource', 'definition'});
            switch description.type
                case 'basis.compiled'
                    object = mdemat.common.approximation.CompiledFunction( ...
                        data.dimensions, data.components, data.resource);
                case 'basis.compiled-orthogonal'
                    object = mdemat.common.spectral.CompiledOrthogonalBasis(data.components, data.resource);
                otherwise
                    family = regexprep(data.resource, '_[0-9]+$', '');
                    object = mdemat.common.interpolation.CompiledInterpolationBasis(data.components, family);
            end
            if object.NDims ~= data.dimensions, invalid('Compiled basis dimension changed.'); end
            definition = struct('schemaVersion', object.Metadata.schemaVersion, ...
                'derivativeOrders', object.Metadata.derivativeOrders, ...
                'functionText', {object.Metadata.functionText});
            if ~isequal(data.definition, definition)
                error('mdemat:persistence:ResourceChanged', ...
                    'Trusted numeric resource %s differs from the saved basis definition.', data.resource);
            end
        case 'basis.recurrence'
            fields(data, {'kind', 'degree', 'normalized', 'lower', 'upper'});
            object = restoreRecurrence(data);
        case 'basis.cardinal'
            fields(data, {'nodes', 'weights'});
            object = mdemat.common.interpolation.CardinalKernelBasis(data.nodes);
            if ~isequal(size(object.Weights), size(data.weights)), invalid('Invalid cardinal weights.'); end
            object.Weights = data.weights;
        case 'basis.rbf'
            fields(data, {'nodes', 'kernel', 'epsilon'});
            object = mdemat.common.interpolation.RbfKernelBasis(data.nodes, ...
                kernel=data.kernel, epsilon=data.epsilon);
        case 'function.separable'
            fields(data, {'factors', 'pattern'});
            object = mdemat.common.approximation.SeparableFunction( ...
                factors=restoreArray(data.factors), pattern=data.pattern);
        case 'function.constant'
            fields(data, {'dimensions', 'value'});
            object = mdemat.data.function.Constant(data.dimensions, data.value);
        case 'function.linear'
            fields(data, {'matrix', 'shift'}); object = mdemat.data.function.Linear(data.matrix, b=data.shift);
        case 'function.maxwellian'
            fields(data, {'dimensions', 'temperature', 'center'});
            object = mdemat.data.function.Maxwellian(data.dimensions, ...
                temperature=data.temperature, center=data.center);
        case {'function.product', 'function.composition'}
            fields(data, {'lhs', 'rhs'});
            if strcmp(description.type, 'function.product')
                object = mdemat.data.function.Product(restore(data.lhs), restore(data.rhs));
            else
                object = mdemat.data.function.Composition(restore(data.lhs), restore(data.rhs));
            end
        otherwise
            error('mdemat:persistence:UnsupportedTag', 'Unknown representation tag: %s.', description.type);
    end
catch exception
    if startsWith(exception.identifier, 'mdemat:persistence:'), rethrow(exception); end
    failure = MException('mdemat:persistence:InvalidDescriptor', 'Representation reconstruction failed.');
    throwAsCaller(addCause(failure, exception));
end
end

function object = restore(description)
object = mdemat.common.fem.restoreSpace(description);
end

function objects = restoreArray(descriptions)
if ~iscell(descriptions), invalid('Representation children must be a cell array.'); end
objects = cell(1, numel(descriptions));
for index = 1:numel(objects), objects{index} = restore(descriptions{index}); end
end

function object = restoreElement(tag, data)
expected = {'geometry', 'approximator'};
if ~strcmp(tag, 'element.base'), expected{end+1} = 'volume'; end
hasFlux = ismember(tag, {'element.broken', 'element.broken-orthotope', 'element.sl-orthotope'});
if hasFlux, expected{end+1} = 'flux'; end
if strcmp(tag, 'element.l2-sphere'), expected{end+1} = 'reduction'; end
fields(data, expected);
geometry = restore(data.geometry); approximator = restore(data.approximator);
if strcmp(tag, 'element.base')
    object = mdemat.common.fem.Element(geometry, approximator); return;
end
volume = restore(data.volume);
if hasFlux, flux = restoreArray(data.flux); flux = [flux{:}]; end
switch tag
    case 'element.l2', object = mdemat.common.fem.L2Element(geometry, approximator, volume);
    case 'element.h1', object = mdemat.common.fem.H1Element(geometry, approximator, volume);
    case 'element.broken', object = mdemat.common.fem.BH1Element(geometry, approximator, volume, flux);
    case 'element.l2-orthotope', object = mdemat.common.fem.L2OrthotopeElement(geometry, approximator, volume);
    case 'element.l2-sphere', object = mdemat.common.fem.L2SphereElement(geometry, approximator, volume, data.reduction);
    case 'element.broken-orthotope'
        object = mdemat.common.fem.BH1OrthotopeElement(geometry, approximator, volume, flux);
    case 'element.sl-orthotope'
        object = mdemat.common.fem.BH1SLOrthotopeElement(geometry, approximator, volume, flux);
end
end

function object = restoreRule(data)
switch data.kind
    case 'GaussLegendreRule', object = mdemat.common.quadrature.GaussLegendreRule(data.dimensions);
    case 'GaussLobattoRule', object = mdemat.common.quadrature.GaussLobattoRule(data.dimensions);
    case 'ClosedNewtonCotesRule', object = mdemat.common.quadrature.ClosedNewtonCotesRule(data.dimensions);
    case 'PeriodicTrapezoidalRule', object = mdemat.common.quadrature.PeriodicTrapezoidalRule(data.dimensions);
    case 'GaussHermiteRule', object = mdemat.common.quadrature.GaussHermiteRule(data.dimensions);
    case 'LebedevRule', object = mdemat.common.quadrature.LebedevRule(coord=data.coordinates);
    case 'GaussTrapezoidalRule'
        object = mdemat.common.quadrature.GaussTrapezoidalRule(data.dimensions, coord=data.coordinates);
    case 'MonteCarloOrthotopeRule'
        object = mdemat.common.quadrature.MonteCarloOrthotopeRule(data.dimensions, seed=data.seed);
    case 'MonteCarloSphereRule'
        object = mdemat.common.quadrature.MonteCarloSphereRule( ...
            data.dimensions, seed=data.seed, coord=data.coordinates);
    otherwise, error('mdemat:persistence:UnsupportedTag', 'Unknown quadrature rule tag.');
end
if object.NDims ~= data.dimensions, invalid('Quadrature dimension changed.'); end
end

function object = restoreRecurrence(data)
switch data.kind
    case 'LegendreBasis'
        object = mdemat.common.spectral.LegendreBasis(data.degree, isNormalized=data.normalized);
    case 'HermiteBasis'
        object = mdemat.common.spectral.HermiteBasis(data.degree, isNormalized=data.normalized);
    case 'ChebyshevBasis'
        object = mdemat.common.spectral.ChebyshevBasis(data.degree, isNormalized=data.normalized);
    case 'FourierBasis'
        object = mdemat.common.spectral.FourierBasis(data.degree, isNormalized=data.normalized, ...
            lower=data.lower, upper=data.upper);
    case 'SphericalHarmonicBasis'
        object = mdemat.common.spectral.SphericalHarmonicBasis(data.degree, isNormalized=data.normalized);
    otherwise, error('mdemat:persistence:UnsupportedTag', 'Unknown recurrence basis tag.');
end
end

function fields(data, expected)
if ~isstruct(data) || ~isscalar(data) || ~isequal(sort(fieldnames(data)), sort(expected(:)))
    invalid('Descriptor fields do not match the registered schema.');
end
end

function invalid(message)
error('mdemat:persistence:InvalidDescriptor', '%s', message);
end

function valid = valueTree(value)
allowed = {'double', 'single', 'logical', 'char', 'string', 'int8', 'uint8', ...
    'int16', 'uint16', 'int32', 'uint32', 'int64', 'uint64'};
if ismember(class(value), allowed)
    valid = true;
elseif iscell(value)
    valid = all(cellfun(@valueTree, value), 'all');
elseif isstruct(value)
    valid = all(cellfun(@valueTree, struct2cell(value)), 'all');
else
    valid = false;
end
end
