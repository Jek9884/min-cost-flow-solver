seed = 42;
filename = "graphs/net8_8_1.dmx";
[A, E, D, b,c]  = utility_read_matrix(filename, seed);

disp(length(b))
disp(length(A))

spy([c;b]);