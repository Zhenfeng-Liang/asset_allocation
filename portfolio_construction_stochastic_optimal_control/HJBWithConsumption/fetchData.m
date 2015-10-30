function [data] = fetchData(stocks, f ,fromDate, endDate, period)

    c = yahoo;
    
    numStocks = length(stocks);
    
    tmp = fetch(c, stocks(1), f, fromDate, endDate, period);
    data = [tmp(:,2)];

    for i = 2:numStocks
        
        tmp = fetch(c, stocks(i), f, fromDate, endDate, period);
        data = [data, tmp(:,2)];
        
    end
        
    close(c);

end
