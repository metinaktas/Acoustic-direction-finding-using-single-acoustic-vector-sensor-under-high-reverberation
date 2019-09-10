function matrix = Vector2Matrix(vector, row, column, convertionType)

% Check the parameters validity
lenVect = length(vector);
if (lenVect ~= row * column)
    matrix = [];
else
    matrix = zeros(row, column);
    switch (convertionType)
        case 'column'
            for rowIndx = 1:row
                indx = (rowIndx-1)*column + [1:column];
                matrix(rowIndx,:) = vector(indx);
            end
        case 'row'
            for colIndx = 1:column
                indx = (colIndx-1)*row + [1:row];
                matrix(:,colIndx) = vector(indx);
            end
        otherwise
            matrix = [];
    end
end