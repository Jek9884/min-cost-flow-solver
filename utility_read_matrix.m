function [A, E, D, b,c] = utility_read_matrix(filename, seed)
    rng(seed);

    graph = readDimacsFile(filename);

    %[E, b] = createAdjacencyMatrix(graph);
    E = createIncidenceMatrix(graph);
    % Dimension of diagonal block
    m = size(E,2); 
    n = graph.nNodes-1;
    A = zeros(m+n, m+n);
    
    % Generate random elements
    diagonal = rand(m, 1);
    
    % Create diagonal matrix
    D = diag(diagonal);
    
    E = E(1:n, 1:end);
    
    A(1:m,1:m) = D;
    A(1:m, m+1:end) = E';
    A(m+1:end, 1:m) = E;
    
    [cost,flows] = get_b(graph);
    b = [flows;cost];
end

function [cost ,flows] = get_b(graph)
    edgeList = graph.edgeList;
    nNodes = graph.nNodes;

    flows = zeros(nNodes,1);
    for i=1:length(graph.nodesList)
        node_id = graph.nodesList(i,1);
        flows(node_id)=graph.nodesList(i,2);
    end

    cost = graph.edgeList(:,4);
end

function graph = readDimacsFile(filename)
    fid = fopen(filename, 'r');
    if fid == -1
        error('Unable to open file.');
    end

    graph = struct();
    graph.edgeList = [];
    graph.nodesList = [];
    
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
            graph.edgeList(end+1, :) = edge;
        elseif tokens{1} == 'n'
            node = str2double(tokens(2:end));
            graph.nodesList(end+1, :) = node;
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

function [indidenceMatrix] = createIncidenceMatrix(graph)
    edgeList = graph.edgeList;
    nNodes = graph.nNodes;
    numEdges = size(edgeList, 1);
    indidenceMatrix = zeros(nNodes, numEdges);

    for i = 1 : numEdges
        vertex1 = edgeList(i, 1);
        vertex2 = edgeList(i, 2);
        indidenceMatrix(vertex1, i) = 1;
        indidenceMatrix(vertex2, i) = -1;
    end
end