%% solving differential equation



function dxdt = deriv(t,x,N,K,n_0,n_i,n_s,S,mu,ro,M_exc,M_int,theta,n,gamma,E_total)

last_term = 0;
    for i =1:N
        for j = 1:N 
            last_term = last_term + x(j)*( (1-theta)*M_exc(j,i) + theta * M_int(j,i));%*();
        end
        
        h = zeros(1,N*2);

        for k = 1:N*2
            for j = 1:3
            % h(i) = h(i) + S(j) * n_s * exp(-gamma*E_total(i,j));
             h(k) = h(k) + S(j) * n_s * exp(-gamma*E_total(k,j));
            end
        end

        
         dxdt(i) = x(i)*(1 - sum(x)/K) - (n_0 + n_i +h(i)+mu + ro)*x(i) + ro * last_term;
         dxdt(i+N) = x(i+N) * (1 - sum(x)/K) - (n_0 + h(i+N)) * x(i+N) + mu * x(i);

    end
    dxdt = dxdt';
end