function averageTime = utility_calculate_avg_time(func, varargin)
    totalTime = 0;
    
    for i = 1:10
        tic;  % Start timer
        func(varargin{:});  % Call the input function with parameters
        elapsedTime = toc;  % Get elapsed time
        totalTime = totalTime + elapsedTime;
    end
    
    averageTime = totalTime / 10;
end