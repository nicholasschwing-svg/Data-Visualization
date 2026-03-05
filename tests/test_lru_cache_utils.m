function test_lru_cache_utils()
    c = mv_lru_create(2);
    c = mv_lru_put(c, 'a', 1);
    c = mv_lru_put(c, 'b', 2);

    [c, v, hit] = mv_lru_get(c, 'a');
    assert(hit);
    assert(v == 1);

    c = mv_lru_put(c, 'c', 3);
    [c, ~, hitB] = mv_lru_get(c, 'b');
    assert(~hitB, 'Expected b to be evicted as LRU');

    [c, vA, hitA] = mv_lru_get(c, 'a');
    [c, vC, hitC] = mv_lru_get(c, 'c');
    assert(hitA && vA == 1);
    assert(hitC && vC == 3);
end
