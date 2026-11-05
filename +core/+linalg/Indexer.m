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
        style % Storage ordering style: 'F' (column-major) or 'C' (row-major)
    end

    methods
        function obj = Indexer(style)
            % INDEXER Constructor for indexer.
            %
            %   obj = Indexer() creates a Indexer with @a style = 'F'.
            %
            %   obj = Indexer(style) creates a Indexer with the specified
            %   @a style.
            %
            % Inputs:
            %   style - Storage ordering style (optional, default: 'F')
            %
            % Outputs:
            %   obj - Constructed Indexer object

            if nargin < 1, style = 'F'; end
            obj.style = style;
        end
    end
end