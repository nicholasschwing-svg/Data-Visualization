function cache = mv_lru_put(cache, key, value)
%MV_LRU_PUT Insert or update an LRU item and evict if needed.
    if isKey(cache.map, key)
        cache.map(key) = value;
        idx = find(strcmp(cache.order, key), 1, 'first');
        if ~isempty(idx)
            cache.order(idx) = [];
        end
        cache.order{end+1} = key;
        return;
    end

    cache.map(key) = value;
    cache.order{end+1} = key;

    while numel(cache.order) > cache.maxItems
        oldKey = cache.order{1};
        cache.order(1) = [];
        if isKey(cache.map, oldKey)
            remove(cache.map, oldKey);
            cache.evictions = cache.evictions + 1;
        end
    end
end
