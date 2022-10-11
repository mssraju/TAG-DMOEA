function wa = WeightAdapationIndicator(wa,Global,ZNad)
    if wa==0 && Global.gen>50
        if pdist2(ZNad(Global.gen,:),ZNad(Global.gen-50,:),'euclidean')<=1e-3 || Global.gen/Global.maxgen>0.5
                wa=1;
        end
    end    
end