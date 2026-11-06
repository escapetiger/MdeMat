classdef RadiationScheme < mdemat.internal.execution.MolScheme
    % RADIATIONSCHEME Base class for all schemes solving radiation
    % transport equations.
    %
    %   RadiationScheme implements a numerical scheme for solving radiative
    %   transfer equation. This base class defines the basic validation and
    %   parse logic of configuration struct.
    %
    % See also:
    %   mdemat.internal.execution.MolScheme

    methods (Access = protected)
        function obj = setup(obj)
            % SETUP Setup default configuration options for radiation
            % schemes.

            %< Call parent setup
            obj = setup@mdemat.internal.execution.MolScheme(obj);
            mdemat.internal.execution.defineStudyInputs(obj,'radiation');

            %< Initial and boundary conditions
            obj.addConfig('ic', default=[], ...
                validator=@(x) isempty(x) || isa(x, 'function_handle') || ...
                    (isnumeric(x) && isscalar(x)));
            obj.addConfig('bc', default=[], ...
                validator=@(x) true);

            %< Time stepping
            obj.addConfig('cfl', default=0.5, ...
                validator=@(x) isempty(x) || (isscalar(x) && x > 0));
            obj.addConfig('dt', default=[], ...
                validator=@(x) isempty(x) || (isscalar(x) && x > 0));

            %< Scaling parameter
            obj.addConfig('epsilon', default=1, ...
                validator=@(x) isnumeric(x) && isscalar(x) && x > 0);

            %< Physical coefficients
            obj.addConfig('scattering', default=[], ...
                validator=@(x) isnumeric(x) || isa(x, 'function_handle') || isempty(x));
            obj.addConfig('absorption', default=[], ...
                validator=@(x) isnumeric(x) || isa(x, 'function_handle') || isempty(x));
            obj.addConfig('source', default=[], ...
                validator=@(x) isnumeric(x) || isa(x, 'function_handle') || isempty(x));

            %< Discretization options
            obj.addConfig('vDimReduction', default='symmetry', ...
                validator=@(x) ismember(x, {'', 'symmetry', 'topology'}));

            %< Analysis and visualization
            obj.addConfig('eId', default='Radiation', ...
                validator=@(x) (ischar(x) || isstring(x)));
            obj.addConfig('errorTypes', default={'L1', 'L2', 'max'}, ...
                validator=@(x) iscell(x) && all(cellfun(@(y) ischar(y) || isstring(y), x)));
            obj.addConfig('exact', default=[], ...
                validator=@(x) isstruct(x) || isa(x, 'function_handle') || isempty(x));
            obj.addConfig('components', default=struct('U', [], 'G', [], 'F', 1), ...
                validator=@(x) isstruct(x) && all(structfun(@(y) isempty(y) || isnumeric(y) && isscalar(y) && (y > 0), x)));
        end

        function obj = configure(obj)
            % CONFIGURE Generate radiation-specific configuration
            % dependencies.

            %< Call parent configure
            obj = configure@mdemat.internal.execution.MolScheme(obj);

            %< Determine boundary condition type
            if isempty(obj.Config.bc)
                obj.setDerived('bcType', 'periodic');
            else
                obj.setDerived('bcType', 'dirichlet');
            end

        end

    end
end
