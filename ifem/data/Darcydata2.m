function pde = Darcydata2
%% 
% p = x^3y+y^4 + sin(pi*x)cos(pi*y)

pde = struct('p',@exactp,'K',@K,'exactu',@exactu,'Dp',@gradp,...
             'f', @f,'g_D',@g_D,'g_N',@g_N);

    function s = K(pt)
       x = pt(:,1); y = pt(:,2);
       s(:,3) = 3*x.*y;
       s(:,1) = 1 + 4*(x.^2 + y.^2);
       s(:,2) = 1 + 11*(x.^2 + y.^2);
    end
    function s = exactp(pt)
       x = pt(:,1); y = pt(:,2);
       s = x.^3.*y + y.^4 + sin(pi*x).*cos(pi*y);
    end
    function s = gradp(pt)
       x = pt(:,1); y = pt(:,2);
       s(:,2) = x.^3 + 4*y.^3 - pi*sin(pi*x).*sin(pi*y);
       s(:,1) = 3*x.^2.*y + pi*cos(pi*x).*cos(pi*y);       
    end
    function s = exactu(pt)
       Dp = -gradp(pt);  % u = - K*grad(p)
       K = K(pt);
       s(:,2) = K(:,3).*Dp(:,1) + K(:,2).*Dp(:,2);
       s(:,1) = K(:,1).*Dp(:,1) + K(:,3).*Dp(:,2);
    end
    function s = f(pt)
       s = 1;
    end
    function s = g_D(pt)
       s = exactp(pt); 
    end
    function s = g_N(pt,n)
       s = dot(exactu(pt),n,2); 
    end
end