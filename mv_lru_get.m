function [cache, value, hit] = mv_lru_get(cache, key)
%MV_LRU_GET Lookup and promote item in LRU cache.
    value = [];
    hit = false;
    if isKey(cache.map, key)
        value = cache.map(key);
        hit = true;
        cache.hits = cache.hits + 1;
        idx = find(strcmp(cache.order, key), 1, 'first');
        if ~isempty(idx)
            cache.order(idx) = [];
        end
        cache.order{end+1} = key;
    else
        cache.misses = cache.misses + 1;
    end
end
