function y = hazardfunc_multi_malignant_cells_SimpleBirth(t,x,p,cell_num,N0,lambda)
    muN = p(1); mu1 = p(2); mu2 = p(3); alpha1=p(4); alpha2=p(5); beta1=p(6); beta2=p(7);
    y = zeros(6,1);
    N = cell_num(ceil(t+1e-16));

    if t==0
        f=0;
        fprime=0;
    else
        f_log = (N-1)*log(1-exp(-lambda * t));
        fprime_log = log(N-1)+log(lambda)-lambda*t+(N-2)*log(1-exp(-lambda*t));
        f=exp(f_log);
        fprime = exp(fprime_log);
    end

    y(1) = muN*N0*x(1)*(x(3)-1);
    y(2) = -muN*N0*x(4); %Cumulative hazard. To get hazard rate you take the derivative
    y(3) = beta1 - (alpha1+beta1+mu1)*x(3)+mu1*x(3)*x(5)+alpha1*x(3)^2;
    y(4) = -(alpha1+beta1+mu1)*x(4)+mu1*x(4)*x(5)+mu1*x(3)*x(6)+2*alpha1*x(3)*x(4);
    y(5) = beta2-(alpha2+beta2+mu2*f)*x(5)+alpha2*x(5)^2;
    y(6) = -(alpha2+beta2+mu2*f)*x(6)+2*alpha2*x(6)*x(5)-mu2*fprime*x(5);
end
