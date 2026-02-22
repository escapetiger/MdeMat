function S = block(As, Ms, Ns)
% BLOCK Constructs a block matrix from a 2D cell array of matrices.
%
%   S = block(As) infers the block dimensions and returns a block matrix
%   from a cell array of matrices.
%
%   S = block(As, Ms, Ns) returns a block matrix from a cell array of
%   matrices with specified row and column dimensions.
%
%   Each nonempty cell in @a As corresponds to a block in @a S. The
%   function automatically detects whether to create sparse or dense
%   output based on input matrix types. When @a Ms and @a Ns are not
%   provided, dimensions are inferred from the first nonempty block in
%   each row and column.

[nBlockRows, nBlockCols] = size(As);

if nargin < 2 || isempty(Ms) || isempty(Ns)
    Ms = zeros(nBlockRows, 1);
    Ns = zeros(nBlockCols, 1);
    
    for iRow = 1:nBlockRows
        found = false;
        for jCol = 1:nBlockCols
            if ~isempty(As{iRow, jCol})
                Ms(iRow) = size(As{iRow, jCol}, 1);
                found = true;
                break;
            end
        end
        
        core.except.assert(found, 'EmptyBlockRow', ...
            ['Block row %d has no nonempty blocks.' ...
            ' Provide Ms explicitly.'], iRow);
    end
    
    for jCol = 1:nBlockCols
        found = false;
        for iRow = 1:nBlockRows
            if ~isempty(As{iRow, jCol})
                Ns(jCol) = size(As{iRow, jCol}, 2);
                found = true;
                break;
            end
        end
        
        core.except.assert(found, 'EmptyBlockColumn', ...
            ['Block column %d has no nonempty blocks.' ...
            ' Provide Ns explicitly.'], jCol);
    end
else
    core.except.assert(numel(Ms) == nBlockRows, ...
        'InvalidRowSizes', ...
        'Length of Ms must equal the number of block rows (%d).', ...
        nBlockRows);
    
    core.except.assert(numel(Ns) == nBlockCols, ...
        'InvalidColumnSizes', ...
        'Length of Ns must equal the number of block columns (%d).', ...
        nBlockCols);
    
    Ms = Ms(:);
    Ns = Ns(:);
end

for iRow = 1:nBlockRows
    for jCol = 1:nBlockCols
        if ~isempty(As{iRow, jCol})
            [rowCount, colCount] = size(As{iRow, jCol});
            
            core.except.assert(rowCount == Ms(iRow), ...
                'RowSizeMismatch', ...
                'Block (%d,%d) has %d rows; expected %d.', ...
                iRow, jCol, rowCount, Ms(iRow));
            
            core.except.assert(colCount == Ns(jCol), ...
                'ColumnSizeMismatch', ...
                'Block (%d,%d) has %d columns; expected %d.', ...
                iRow, jCol, colCount, Ns(jCol));
        end
    end
end

totalRows = sum(Ms);
totalCols = sum(Ns);
rowOffsets = [0; cumsum(Ms)];
colOffsets = [0; cumsum(Ns)];

useSparse = true;
for iRow = 1:nBlockRows
    for jCol = 1:nBlockCols
        if ~isempty(As{iRow, jCol}) && ~issparse(As{iRow, jCol})
            useSparse = false;
            break;
        end
    end
    if ~useSparse
        break;
    end
end

if useSparse
    totalNnz = 0;
    for iRow = 1:nBlockRows
        for jCol = 1:nBlockCols
            if ~isempty(As{iRow, jCol})
                totalNnz = totalNnz + nnz(As{iRow, jCol});
            end
        end
    end

    rowIndices = zeros(totalNnz, 1);
    colIndices = zeros(totalNnz, 1);
    values = zeros(totalNnz, 1);
    nnzPos = 1;

    for iRow = 1:nBlockRows
        for jCol = 1:nBlockCols
            if ~isempty(As{iRow, jCol})
                block = As{iRow, jCol};
                [blockRows, blockCols, blockVals] = find(block);
                nValues = numel(blockVals);
                if nValues > 0
                    indexRange = nnzPos:(nnzPos + nValues - 1);
                    rowIndices(indexRange) = blockRows + rowOffsets(iRow);
                    colIndices(indexRange) = blockCols + colOffsets(jCol);
                    values(indexRange) = blockVals;
                    nnzPos = nnzPos + nValues;
                end
            end
        end
    end

    S = sparse(rowIndices, colIndices, values, totalRows, totalCols);
else
    S = zeros(totalRows, totalCols);
    for iRow = 1:nBlockRows
        for jCol = 1:nBlockCols
            if ~isempty(As{iRow, jCol})
                block = As{iRow, jCol};
                if issparse(block)
                    block = full(block);
                end
                rowRange = (rowOffsets(iRow) + 1):rowOffsets(iRow + 1);
                colRange = (colOffsets(jCol) + 1):colOffsets(jCol + 1);
                S(rowRange, colRange) = block;
            end
        end
    end
end
end