function [ value ] = get_cs_value( cs )
    global UE;
    global dat;
    global best;
    global best_cs;
    eff = 0;
    for j = 1:length(cs)
        a = UE([cs{j}],:);
        d1 = (a(1,1)^2 + a(1,2)^2)^0.5;
        d2 = (a(2,1)^2 + a(2,2)^2)^0.5;
        pair_position = floor([d1 d2]/100)*100;
%             pair_position;
        pl1 = mapDist_PL(pair_position(1));
        pl2 = mapDist_PL(pair_position(2));

        for k = 1:size(dat,1)
            if((dat(k,1) == pl1 && dat(k,2)==pl2) || (dat(k,2) == pl1 && dat(k,1)==pl2))

                pair_val = dat(k,3)*2*dat(k,5)+dat(k,4)*2*dat(k,6) ;               
                eff = eff + pair_val;
            end
        end
    end
    if (eff > best)
        best = eff;
        best_cs = cs;
    end
    value = eff;

end

function [ pl ] = mapDist_PL( dist )
    d_ref = 400:100:1900;
    pl_ref = [116.118 119.149 121.626 123.72 125.534 127.134 128.566 129.861 131.043 132.13 133.137 134.074 134.951 135.775 136.551 137.286];
    id = find(d_ref==dist);
    pl = pl_ref(id);
end
