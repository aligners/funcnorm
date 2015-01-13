function p = findPercentile(s, pct)
    s2 = s.^2;
    cs2 = cumsum(s2)/sum(s2);
    locs = find(cs2 >= pct/100);
    p = locs(1);
end   
