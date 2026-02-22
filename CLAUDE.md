# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MdeMat is a modular and scalable MATLAB library for numerical algorithms research on differential equations. The project follows "The Elements of MATLAB Style" by Richard K. Johnson.

## Architecture

The codebase is organized in a strict layered architecture with dependency constraints:

* Level 1 (Core): +core/
* Level 2 (Approximation): +approx/
* Level 3 (Physics): +physics/

**Critical Rule**: Each level can only depend on lower levels. Higher levels cannot depend on same or higher levels.

### Core Package Structure (`+core/`)

The `+core/` package contains fundamental components that higher levels depend on:

* **`+function/`**: Function hierarchy and abstract function definitions
* **`+geometry/`**: Geometric primitives (Hyperball, Hypersphere, Orthotope, etc.)
* **`+linalg/`**: Linear algebra utilities and matrix operations
* **`+optim/`**: Optimization algorithms and solvers
* **`+symbolic/`**: Symbolic computation helpers
* **`+chrono/`**: Performance measurement and timing utilities
* **`+except/`**: Exception handling and assertion functions

### Approximation Package Structure (`+approx/`)

The `+approx/` package implements numerical approximation methods and depends only on `+core/`:

* **`+assembly/`**: Assembly operations for finite element methods. Base class is `FiniteElementAssembly`; subclasses drop the prefix: `AdjointAssembly`, `EllipticAssembly`, `SourceAssembly`, `UpwindAssembly`, `MultiplierAssembly`, `SemiLagrangianAssembly`, `WeightedAdjointAssembly`
* **`+element/`**: Finite element definitions. Base class `Element` with hierarchy: `H1Element`, `L2Element`, `BH1Element`, `BH1OrthotopeElement`, `BH1SLOrthotopeElement`, `L2OrthotopeElement`, `L2SphereElement`
* **`+integrate/`**: Numerical integration schemes (GaussLegendreRule, GaussLobattoRule, etc.)
* **`+linear/`**: Approximation solvers: `LinearApproximator` (base), `ModalApproximator`, `NodalApproximator`
* **`+mesh/`**: Mesh generation and manipulation tools
* **`+odeint/`**: ODE integration methods. Base class `OdeIntegrator` with hierarchy: `DirkIntegrator` (SDIRK2/3/4, ESDIRK3), `BdfIntegrator` (BDF2/3), `BeIntegrator`, `FeIntegrator`, `ExrkIntegrator` (EXRK4, SSPRK2/3, Heun), `ImexrkIntegrator` (ARS111/222/443), `SeparableOdeIntegrator` (LieTrotter, Strang), `EmpIntegrator`, `ImpIntegrator`. Uses positional arguments `step(L, S, M)` for performance. Has `+scripts/` with 6 example scripts (Prothero-Robinson, heat equation, convergence analysis, Van der Pol, stability regions, custom step functions)
* **`+space/`**: Function spaces and basis function definitions. Key classes:
  * `LinearApproximator`: Base approximator. `eval(X, C)` computes `Y = C.' * B` where `C` is `(NDofs × NData)`, `B = basis.eval(X)` is `(NDofs × NPoints)`, and `Y` is `(NData × NPoints)`.
  * `SumSpace`: Container for two `SpectralSpace` objects (`Lhs` + `Rhs`). `evalLhs(X, C)` multiplies `LhsWeight .* C` before evaluating Lhs; `evalRhs` same for Rhs. `invRhs(Y, X, CLhs)` removes macro contribution `Lhs.eval(X, CLhs)` from `Y`, then divides by `RhsWeight` to recover micro coefficients.

### Physics Package Structure (`+physics/`)

The `+physics/` package implements specific physics simulations and depends on both `+core/` and `+approx/`:

* **`+state/`**: State classes (`State` base, `SpatialState`, `KineticState`)
* **`+visual/`**: Visualization infrastructure (`Visualizer`, strategy classes, `FigureManager`, `Database`, `Dataset`)
* **`+analysis/`**: Analysis and post-processing (`Analyzer`, `FiniteElementAbsoluteMetric`, `FiniteElementRichardsonMetric`)
* **`+scheme/`**: Scheme infrastructure (`Scheme` base → `OdeScheme`/`SteadyScheme` → `MolScheme`, plus `ConfigParser`)
* **`+advection/`**: Advection equation solvers (`AdvectionScheme` base, `UwdgScheme`, `SldgScheme`)
* **`+diffusion/`**: Diffusion equation solvers (`DiffusionScheme` base, `LdgScheme`)
* **`+poisson/`**: Poisson equation solvers (`PoissonScheme` base, `LdgScheme`, `PdgScheme`)
* **`+radiation/`**: Radiation transport methods (`RadiationScheme`, `MmDgScheme`, `MacroMicroState`)
* **`+vlasov/`**: Vlasov equation implementations

Each physics subpackage contains `+scripts/` directories with example implementations and `config.txt` files for simulation parameters.

### Radiation Transport Architecture (`+physics/+radiation/`)

The macro-micro DG scheme decomposes the kinetic distribution `f = u + g` where `u` is the macroscopic (modal spherical harmonics) component and `g` is the microscopic (nodal RBF) component. Understanding the data flow across these three classes is essential:

**Velocity discretization (`SumSpace`)**: Built by the example scripts as `SumSpace(vMacroDisc, vMicroDisc)` where `vMacroDisc = SpectralSpace(L2SphereElement.modal(..., nu))` (Lhs) and `vMicroDisc = SpectralSpace(L2SphereElement.nodal(..., nv))` (Rhs). The Lhs and Rhs can have **independent** quadrature nodes — macro nodes chosen for exactness of macro-macro transport, micro nodes for RBF resolution. The LhsWeight/RhsWeight encode the decomposition scaling (ε powers).

**`MacroMicroState` lazy cache pattern**: Boundary evaluations are split across separate assembly phases (`macroTraceEvaluate` for `bU`, then `microTraceEvaluate` for `bG`). A bitmask `Status` tracks which components have been computed:

* `Cache.C` — macro modal coefficients `(NGlobalDofs × NMacroDofs)`
* `Cache.F` — function values at micro (Rhs) velocity nodes `(NTraceDofs × NMicroDofs)`
* `IsMacroComputed = 0b01`, `IsMicroComputed = 0b10`, `IsAllComputed = 0b11`
* After both bits are set, `reset()` is called automatically to clear stale cache
* **Critical**: `MmDgScheme.addBc` must call `state.reset()` between the `bU` and `bG` assembly loops to prevent stale cache from the last macro face contaminating the first micro face evaluation

**`lazyTraceEvaluate` dispatch**: When both macro (`m > 0`) and micro (`n > 0`) DOFs are present and no explicit velocity is supplied, the method evaluates the BC function `f` at **both** macro and micro velocity nodes separately — `F_M` for macro projection, `F_u` stored in `Cache.F` for micro extraction. This avoids the dimension mismatch that occurs when macro and micro node counts differ.

## Development Commands

### Running Tests

Tests use MATLAB's unittest framework and are located in `+tests/` subdirectories within each package. Test classes follow the naming pattern `Test<ClassName>`:

```matlab
% Run all tests for a specific package
runtests('+core/+geometry/+tests')

% Run a specific test class
runtests('+core/+geometry/+tests/TestHyperball')

% Run all tests in core package (recursive)
runtests('+core', 'IncludeSubfolders', true)

% Run all tests in approx package
runtests('+approx', 'IncludeSubfolders', true)

% Run all tests in physics package
runtests('+physics', 'IncludeSubfolders', true)
```

**Test Coverage**: Tests are available for core geometry (`+core/+geometry/+tests`), chrono (`+core/+chrono/+tests`), linear algebra (`+core/+linalg/+tests`), optimization (`+core/+optim/+tests`), mesh generation (`+approx/+mesh/+tests`), scheme infrastructure (`+physics/+scheme/+tests`), and visualization (`+physics/+visual/+tests`).

### Running Physics Examples

Physics packages include example scripts that demonstrate usage. Each physics package contains `+scripts/` directories with examples and a `config.txt` configuration file:

```matlab
% Run advection examples
run('+physics/+advection/+scripts/ex01_periodic.m')
run('+physics/+advection/+scripts/ex02_dirichlet.m')
run('+physics/+advection/+scripts/ex03_discontinuous.m')

% Run diffusion examples
run('+physics/+diffusion/+scripts/ex01_periodic.m')
run('+physics/+diffusion/+scripts/ex02_dirichlet.m')
run('+physics/+diffusion/+scripts/ex03_inhomogenous.m')

% Run Poisson examples
run('+physics/+poisson/+scripts/ex01_periodic.m')
run('+physics/+poisson/+scripts/ex02_dirichlet.m')

% Configuration files are located in +scripts/config.txt for each physics package
% These files contain simulation parameters like grid resolution, time steps, etc.
```

**Configuration File Format**: Each `config.txt` uses `OPTION=key=value` pairs. Required keys for the radiation scheme:

* `schemeName` — `'mmdg'`
* `tOdeIntName` — `'be'`, `'sdirk2'`, `'sdirk3'`, `'bdf2'`, etc.
* `tSplitName` — `''`, `'lie_trotter'`, or `'strang'`
* `xBasisOrder` — polynomial degree (integer)
* `xBasisType` — `'modal'` or `'nodal'`
* `xBasisPattern` — `'P'` (simplex) or `'Q'` (tensor)
* `xPenaltyType` — 2×2 cell array `{'boundary','right';'left',''}` controlling upwind flux directions
* `cfl` — CFL number for automatic time step selection

**Physics example script structure**: All function handles (`scattering`, `absorption`, `source`, `ic`, `bc`, `exact`) are assigned **inside the `run(config)` function**, not at the top level. The top-level loop only sets scalar parameters (`epsilon`, `order`, checkpoint paths) and calls `run`. Helper functions (`fInit`, `fBc`, `fScattering`, etc.) are defined at the end of the script and capture `config` via closure.

**Plotting `plt` struct**: Passed to `plotDensity1d`/`plotDensity2d` with fields: `figIdx` (figure number), `strategy` (`Strategy1d()` or `Strategy2d()`), `style` (`'line'` or `'scatter'`), `colorIdx`, `markerIdx`, `legend` (cell array), `time` (scalar or empty).

### MATLAB to C++ Conversion

Convert MATLAB files to C++ with Doxygen documentation:

```bash
perl m2cpp.pl path/to/file.m
```

This script converts MATLAB docstrings to Doxygen format following the project's documentation standards.

### Documentation Generation

Generate API documentation using Doxygen:

```bash
doxygen Doxyfile
```

This creates HTML and LaTeX documentation in the `doc/` directory from MATLAB docstrings and converted C++ files.

### Project Configuration

The `config/` directory contains important project standards:

* **`std.md`**: Complete MATLAB coding standards and documentation templates (required reading for all development)
* **`The Elements of MATLAB Style.pdf`**: Reference book for MATLAB best practices

## Coding Standards

### Documentation Requirements

All code must follow the templates in `config/std.md`:

**Function Documentation:**

```matlab
function [output1, output2] = functionName(input1, input2, varargin)
% FUNCTIONNAME One-line description of function purpose.
%
%   [output1, output2] = functionName(input1, input2) performs a 
%   specific computation or operation.
%
%   [output1, output2] = functionName(input1, input2, option1) performs 
%   another computation or operation.
%
%   Any special considerations, limitations, or implementation details.
%
%   Detailed description should explain the behavior, the steps
%   performed, and how it uses @a input1 and @a input2 to produce @a output1
%   and @a output2.
%
%   Example math: \f$ area = pi * r^2 \f$ or display math:
%
%   \f[
%       E = mc^2
%   \f]
%
% See also:
%   otherFunction1, otherFunction2

end
```

**Class Documentation:**

```matlab
classdef ClassName < ParentClass
    % CLASSNAME One-line description of class purpose.
    %
    %   ClassName represents [short description of functionality].
    %   The detailed description should explain major capabilities,
    %   usage patterns, and key behaviors.
    %
    % Notes: (optional)
    %   Any special considerations, limitations, or implementation details.
    %
    % See also:
    %   RelatedClass1, relatedFunction1

    properties
        Property1 % Regular property, use PascalCase
    end
    
    properties (Dependent)
        Property2 % Dependent property, computed from others
    end

    methods
        function obj = ClassName(arg1, arg2)
            % CLASSNAME Construct an instance of ClassName with a complete description.
            %
            %   obj = ClassName(arg1) creates an object with @a arg1. 
        end

        function output = publicMethod(obj, input)
            % PUBLICMETHOD Performs a specific operation and returns the result.
            %
            %   [output] = publicMethod(input) performs a specific computation or operation.
        end
    end

    methods
        function obj = set.Property1(obj, value)
            % SET.PROPERTY1 Sets the value of the property 'Property1'.
        end

        function value = get.Property2(obj)
            % GET.PROPERTY2 Returns the value of the dependent property 'Property2'.
        end
    end

    methods (Access = protected)
        function output = protectedMethod(obj, input)
            % PROTECTEDMETHOD One-line description only.
        end
    end

    methods (Access = private)
        function output = privateMethod(obj, input)
            % PRIVATEMETHOD One-line description only.
        end
    end
end
```

### Naming Conventions

* **Functions**: camelCase verbs (`computeArea`, `normalizeData`)
* **Classes**: PascalCase (`GeometryManager`, `BasisFunction`)  
* **Properties**: PascalCase (`Radius`, `MaxIterations`) - follows MATLAB toolbox conventions
* **Variables**: camelCase (`radius`, `maxIterations`)
* **Constants**: Prefer PascalCase (`MaxEpochs`), UPPERCASE for math constants (`PI`, `MAX_ITERATIONS`)

### Code Style Rules

* Use 4 spaces for indentation (no tabs)
* Maximum 80 characters per line
* Spaces around operators (`=`, `+`, `-`)
* Math formulas: Use LaTeX format `\f$ formula \f$` for inline, `\f[\n formula \n\f]` for display
* No in-function comments (explanations go in docstrings)
* One-line docstrings for private/protected methods only

### Vectorization

* Prefer vectorized operations over loops for performance
* Avoid growing arrays inside loops; preallocate memory when necessary
* Use built-in MATLAB functions for computations when possible

### Argument Validation

* **All public methods and functions** (except property getters) must use an `arguments` block for input validation and default values
* Do not use manual `nargin` checks or `inputParser` unless unavoidable
* Always use `name=value` syntax for optional arguments when calling functions

Example:

```matlab
function result = computeNorm(x, p)
% COMPUTENORM Compute the p-norm of vector x.

    arguments
        x double {mustBeVector}
        p double {mustBePositive} = 2
    end

    result = sum(abs(x).^p)^(1/p);
end
```

### Performance: Named Arguments vs Positional Arguments

**Critical**: MATLAB's named arguments (`options.X`) have significant overhead when called in tight loops, especially when parameters are large sparse matrices. The overhead comes from MATLAB copying and validating the options struct on every call.

**Rule**: In compute-intensive modules (`+odeint/`, `+linalg/`, inner assembly loops), use **positional arguments** instead of named arguments for methods called repeatedly during time stepping. Named arguments are fine for setup/configuration methods called once.

```matlab
% GOOD for hot path: positional arguments
function U = step(obj, L, S, M)
function x = solve(obj, A, b)

% OK for setup: named arguments
function obj = setTimeStep(obj, options)
```

The `+approx/+odeint/` module and `core.linalg.LinearSolver.solve` use positional arguments for this reason. The reference implementation is in `MdeMat_v2/+approx/+odeint/`.

## Working with the Codebase

### Adding New Components

1. **Identify the correct architectural level** for your component
2. **Respect dependency rules**: Only depend on lower architectural levels
3. **Use the documentation templates** from `config/std.md` exactly
4. **Add comprehensive tests** using MATLAB's unittest framework

### Package Organization

Each package should contain:

* Implementation files (`.m`)
* Test files in `+tests/` subdirectory
* Private utilities in `private/` subdirectory when needed

### Reference Implementation (MdeMat_v2)

The `MdeMat_v2/` directory contains the previous version of the library with the same package structure. Use it as a reference for algorithm design and interface patterns, particularly for `+odeint/` integrators.

### Development Directory

The `dev/` directory contains experimental and development components that may eventually be promoted to the core architecture:

* **`+approx/`**: Experimental approximation methods
* **`+radiation/`**: Radiation transport research implementations
* **`+vlasov/`**: Advanced Vlasov equation implementations

**Development Scripts**: The `dev/` directory also contains development analysis scripts (dev01_*.m through dev06_*.m) for numerical experiments with RBF methods, interpolation analysis, discrete ordinates, partition of unity methods, encoder-decoder approaches, and orthogonal decomposition.

### Legacy Directory

The `legacy/` directory at the project root contains older physics implementations (e.g., `RadiationScheme.m`, `MacroMicroState.m`, `MmDgScheme.m`). Some packages also have `legacy/` subdirectories (e.g., `+approx/+assembly/legacy/`) with deprecated implementations.

### DistMesh Integration

The `distmesh/` directory contains mesh generation algorithms based on the DistMesh library. It includes:

* **MATLAB implementations** (`.m` files): Core mesh generation functions like `distmesh2d.m`, `distmeshnd.m`
* **Optimized C++ implementations** (`.cpp` files with compiled MEX files): Performance-critical functions including:
  * `dellipse.cpp` - Ellipse distance function
  * `dellipsoid.cpp` - Ellipsoid distance function  
  * `dsegment.cpp` - Line segment distance function
  * `trisurfupd.cpp` - Surface triangulation updates
* **Demo scripts**: `meshdemo2d.m`, `meshdemond.m` for testing mesh generation

## Important Development Guidelines

### File Management

* **NEVER create files** unless absolutely necessary for achieving your goal
* **ALWAYS prefer editing** an existing file to creating a new one
* **NEVER proactively create documentation files** (*.md) or README files unless explicitly requested

## IDE Configuration

### VSCode Settings

The project includes VSCode configuration in `.vscode/settings.json` with Claude Code permissions properly configured. This ensures safe development practices with appropriate guardrails.
