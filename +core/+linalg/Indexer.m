classdef Indexer < handle
    % INDEXER Base class for all indexers.
    %
    %   Indexer provides foundation for classes that generate indices
    %   based on specific schemes. Defines basic properties and validation
    %   methods common to all indexers.
    %
    % See Also:
    %   core.linalg.MultiIndexer, core.linalg.CachedMultiIndexer

    properties
        Style % Storage ordering style: 'F' (column-major) or 'C' (row-major)
    end

    methods
        function obj = Indexer(options)
            % INDEXER Construct an instance of Indexer.
            %
            %   obj = Indexer() creates a Indexer with @a style = 'F'.
            %
            %   obj = Indexer(Name=Value) creates a Indexer with the
            %   specified options.

            arguments
                options.style{mustBeTextScalar} = 'F'
            end

            obj.Style = options.style;
        end
    end
end