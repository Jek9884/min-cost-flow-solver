function [E, D, b] = utility_read_matrix(filename, seed, print_size)
    rng(seed); 

    graph = readDimacsFile(filename);

    %[E, b] = createAdjacencyMatrix(graph);
    E = createIncidenceMatrix(graph);

    %graph.nNodes = graph.nNodes-1;

    % Dimension of diagonal block
    m = size(E,2); 
    n = graph.nNodes;
    
    % Generate random elements
    D = rand(m, 1);
    %for i = 1:m
    %    D(i) = poissrnd(5);
    %end
 
    E = E(1:n-1, 1:end);
    
    [cost,flows] = get_b(graph);
    b = [cost;flows];

    if print_size
        fprintf("Number of nodes: %d\n Number of edges: %d\n", n-1, m)
    end
end

function [cost ,flows] = get_b(graph)
    
    nNodes = graph.nNodes;

    flows = zeros(nNodes-1,1); 
    for i=1:length(graph.nodesList)-1
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