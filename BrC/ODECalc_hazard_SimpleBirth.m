function v=ODECalc_hazard2(p,T,M,cell_num,Age,N0,bm,Menarche_Age)
    warning('off', 'all')
    if length(Age)<=1
        tt = 0:length(M)-1;
        INC_data = M;
    elseif length(Age)>1 && max(Age)>=T
        tt = min(Age):T;
        INC_data = M(1:length(tt));
    elseif length(Age)>1 && max(Age)<T
        tt = Age';
        T=max(Age);
        INC_data = M;
    end
    
    try 
        full_tt = 0:T;
        x_ic=[1 0 1 0 1 0];
        sol = ode15s(@(t,x,p) hazardfunc_multi_malignant_cells_SimpleBirth(t,x,p,cell_num,N0,bm,Menarche_Age),full_tt, x_ic, [], p);
        y = deval(sol,full_tt);
        INC = y(2,:);
        data_term = norm(abs(INC(ismember(full_tt, tt))' - INC_data));
        regularization_term = norm(INC(~ismember(full_tt, tt)));
        v = data_term;
    catch
        v=1e30;
    end

end
