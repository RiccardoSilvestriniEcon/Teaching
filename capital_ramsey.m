%% VFI capital, Ramsey model

clear all
close all

%% Calibration

par.alpha=1/3;
par.betta=0.9;
par.delta=0.1;
par.gamma=2;

k_SS=(par.alpha*par.betta/(1-par.betta*(1-par.delta)))^(1/(1-par.alpha));

% par.Nk=101;
% par.Kini=0.25*k_SS;
% par.Kfin=1.75*k_SS;
% 
% k=zeros(1,par.Nk);
% 
% for i=1:par.Nk
%     k(i)=par.Kini+(i-1)*(par.Kfin-par.Kini)/par.Nk;
% end

k=1.484:0.01:2.484;
par.Nk=length(k);

V=zeros(1,par.Nk);

%% VFI

Vnew        = V;
k_ind       = zeros(1,par.Nk);
Vsupp       = zeros(par.Nk,par.Nk);

% Value function iteration
crit        = 100;
while crit  > 1e-6   
    for i = 1:par.Nk 
        Vsupp(i,:) = (((1-par.delta)*k(i)+k(i)^(par.alpha)-k).^(1-par.gamma)-1)/(1-par.gamma)...
            +par.betta*V;
    end
    %check feasibility
    for i = 1:par.Nk
        for j = 1:par.Nk
            if k(j)-(1-par.delta)*k(i)-k(i)^(par.alpha)>0
                Vsupp(i,j)=-10^9;
            end
        end
    end
    %
    for i = 1:par.Nk 
        [Vnew(1,i),k_ind(1,i)] = max(Vsupp(i,:));
    end
    % convergence criterion
    crit        = sqrt(mean((Vnew-V).^2));
    % update of firm value
    V           = Vnew;
    disp(crit)
end

save V V
