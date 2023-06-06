function [A, E, D, b] = utility_read_matrix(filename, seed)

    rng(seed);

    graph = readDimacsFile(filename);

    [E, b] = createAdjacencyMatrix(graph);
    % Dimension of diagonal block
    m = size(E,2); 
    n = length(E)-1;
    A = zeros(m+n, m+n);
    
    % Generate random elements
    diagonal = rand(m, 1);
    
    % Create diagonal matrix
    D = diag(diagonal);
    
    E = E(1:n, 1:end);
    
    A(1:m,1:m) = D;
    A(1:m, m+1:end) = E';
    A(m+1:end, 1:m) = E;
    
    %b = rand(m+n,1);
    
end


function graph = readDimacsFile(filename)
    fid = fopen(filename, 'r');
    if fid == -1
        error('Unable to open file.');
    end

    graph = struct();
    graph.edgeList = [];
    
    edgeCount = 0;
    cCount = 0;

    line = fgetl(fid);
    
    while ~(line == -1)
        tokens = strsplit(line);
        if tokens{1} == 'c'
            cCount = cCount + 1;
            if cCount == 5
                line = strrep(line, ' ', '');
                line = split(line, ':');
                nNodes = line(2);
                nNodes = str2num(nNodes{1});
            end
        elseif tokens{1} == 'a'
            edge = str2double(tokens(2:end));
            edgeCount = edgeCount + 1;
            graph.edgeList(edgeCount, :) = edge;
        end
        line = fgetl(fid);
    end
    
    graph.nNodes = nNodes;
    fclose(fid);
end

function [adjacencyMatrix, b] = createAdjacencyMatrix(graph)
    edgeList = graph.edgeList;
    nNodes = graph.nNodes;
    adjacencyMatrix = zeros(nNodes);
    numEdges = size(edgeList, 1);

    b = [];

    for i = 1:numEdges
        vertex1 = edgeList(i, 1);
        vertex2 = edgeList(i, 2);
        if i <= nNodes
            b(end+1) = edgeList(i,3); %taking supply vector (c)
        else
            b(end+1) = edgeList(i,4); %taking vector of costs (b)
        end
        adjacencyMatrix(vertex1, vertex2) = adjacencyMatrix(vertex1, vertex2) + 1;
    end
end