function [ result ] = csgDive( UE_idx, cs)
    if (isempty(UE_idx))
        get_cs_value(cs);
        result = 0;
        return
    end
    [piv, piv_idx] = min(UE_idx);
    for i = 1:size(UE_idx,2)
        if (UE_idx(i) == piv)
            continue;
        end
        pair = [piv UE_idx(i)];
        tmp = UE_idx;
        tmp([piv_idx i]) = [];
        
        cs_tmp = [cs pair];
%         pause;
        cs_set = csgDive( tmp , cs_tmp);
    end
    result = 0;
end

function [ pl ] = mapDist_PL( dist )
    d_ref = 400:100:1900;
    pl_ref = [116.118 119.149 121.626 123.72 125.534 127.134 128.566 129.861 131.043 132.13 133.137 134.074 134.951 135.775 136.551 137.286];
    id = find(d_ref==dist);
    pl = pl_ref(id);
end

