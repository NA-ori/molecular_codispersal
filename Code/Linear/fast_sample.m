function y = fast_sample(v, w)
    f = find(w);
    w = w(f);
    w = w / sum(w);
    w = [0 cumsum(w)];
    % y = f(lookup(w, rand));
    y = f(find (w <= rand, 1, "last"));
    y = v(y);
end