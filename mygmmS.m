function [SS, nn_all] = mygmmS(vars_tilde, t, ivf_awt, d20, d02, d11, d3, ...
                              delta_20, delta_02, delta_11, delta_3, ...
                              n20, n02, n11, n3)
    % Function mygmmS is for mainGMMIVF0410.m
    % Extract parameters
    theta_20 = vars_tilde(1);
    theta_02 = vars_tilde(2);
    theta_11 = vars_tilde(3);    
    gamma1   = vars_tilde(4);
    gamma2   = vars_tilde(5);
    gamma3   = vars_tilde(6);
    gamma4   = vars_tilde(7);
    
    % Compute theta_3
    theta_3 = 1 - theta_20 - theta_02 - theta_11;
    Score_t=[];

    for i = 1:length(t)
        % Choose correct gamma_t
        if t(i) == 1
            gamma_t = gamma1;
        elseif t(i) == 2
            gamma_t = gamma2;
        elseif t(i) == 3
            gamma_t = gamma3;
        elseif t(i) == 4
            gamma_t = gamma4;
        end

        % Compute Moment Conditions
        e20_1 = 1- ((1 - ivf_awt(i) * gamma_t) * d20(i) + theta_20 * ivf_awt(i) * gamma_t);        
        e02_1 = 1- ((1 - ivf_awt(i) * gamma_t) * d02(i) + theta_02 * ivf_awt(i) * gamma_t);
        e11_1 = 1- ((1 - ivf_awt(i) * gamma_t) * d11(i) + theta_11 * ivf_awt(i) * gamma_t);        
        e3_1  = 1-  ((1 - ivf_awt(i) * gamma_t) * d3(i)  + theta_3  * ivf_awt(i) * gamma_t);
        e20_0 = 0- ((1 - ivf_awt(i) * gamma_t) * d20(i) + theta_20 * ivf_awt(i) * gamma_t);
        e02_0 = 0- ((1 - ivf_awt(i) * gamma_t) * d02(i) + theta_02 * ivf_awt(i) * gamma_t);
        e11_0 = 0- ((1 - ivf_awt(i) * gamma_t) * d11(i) + theta_11 * ivf_awt(i) * gamma_t);
        e3_0  = 0-  ((1 - ivf_awt(i) * gamma_t) * d3(i)  + theta_3  * ivf_awt(i) * gamma_t);

        moment_20_1 = - e20_1 * ivf_awt(i) * gamma_t;      
        moment_02_1 = - e02_1 * ivf_awt(i) * gamma_t;
        moment_11_1 = - e11_1 * ivf_awt(i) * gamma_t;
        moment_3_1  =    e3_1   * ivf_awt(i) * gamma_t;
        moment_20_0 = - e20_0 * ivf_awt(i) * gamma_t;      
        moment_02_0 = - e02_0 * ivf_awt(i) * gamma_t;
        moment_11_0 = - e11_0 * ivf_awt(i) * gamma_t;
        moment_3_0  =    e3_0   * ivf_awt(i) * gamma_t;

        moment_r20_1 = (d20(i)-theta_20)*ivf_awt(i)*e20_1;
        moment_r02_1 = (d02(i)-theta_02)*ivf_awt(i)*e02_1;
        moment_r11_1 = (d11(i)-theta_11)*ivf_awt(i)*e11_1;
        moment_r3_1  =  (d3(i) -theta_3)  *ivf_awt(i) *e3_1;
        moment_r20_0 = (d20(i)-theta_20)*ivf_awt(i)*e20_0;
        moment_r02_0 = (d02(i)-theta_02)*ivf_awt(i)*e02_0;
        moment_r11_0 = (d11(i)-theta_11)*ivf_awt(i)*e11_0;
        moment_r3_0  =  (d3(i) -theta_3) *ivf_awt(i)*e3_0;
        
        
        if t(i) == 1
           ntot= n20(i) +n02(i) +n11(i) + n3(i);
           A1_1 = [moment_20_1, 0, 0, moment_r20_1, 0, 0, 0];
           A1_0 = [moment_20_0, 0, 0, moment_r20_0, 0, 0, 0];
           A2_1 = [0, moment_02_1, 0, moment_r02_1, 0, 0, 0];
           A2_0 = [0, moment_02_0, 0, moment_r02_0, 0, 0, 0];
           A3_1 = [0, 0, moment_11_1, moment_r11_1, 0, 0, 0];
           A3_0 = [0, 0, moment_11_0, moment_r11_0, 0, 0, 0];
           %A4_1 = [moment_3_1, moment_3_1, moment_3_1, moment_r3_1, 0, 0, 0];
           %A4_0 = [moment_3_0, moment_3_0, moment_3_0, moment_r3_0, 0, 0, 0];
           AA1_1 = repmat(A1_1, n20(i), 1);
           AA2_1 = repmat(A2_1, n02(i), 1);
           AA3_1 = repmat(A3_1, n11(i), 1);
           %AA4_1 = repmat(A4_1, n3(i), 1);
           AA1_0 = repmat(A1_0, ntot-n20(i), 1);
           AA2_0 = repmat(A2_0, ntot-n02(i), 1);
           AA3_0 = repmat(A3_0, ntot-n11(i), 1);
           %AA4_0 = repmat(A4_0, ntot-n3(i), 1);
           %AA=[AA1_1;AA2_1;AA3_1;AA4_1;AA1_0;AA2_0;AA3_0;AA4_0];
           AA=[AA1_1;AA2_1;AA3_1;AA1_0;AA2_0;AA3_0];
        elseif t(i) == 2
           ntot= n20(i) +n02(i) +n11(i)+n3(i);
           A1_1 = [moment_20_1, 0, 0, 0,  moment_r20_1, 0, 0];
           A1_0 = [moment_20_0, 0, 0, 0,  moment_r20_0, 0, 0];
           A2_1 = [0, moment_02_1, 0, 0,  moment_r02_1, 0, 0];
           A2_0 = [0, moment_02_0, 0, 0,  moment_r02_0, 0, 0];
           A3_1 = [0, 0, moment_11_1, 0,  moment_r11_1, 0, 0];
           A3_0 = [0, 0, moment_11_0, 0,  moment_r11_0, 0, 0];
           %A4_1 = [moment_3_1, moment_3_1, moment_3_1,  0, moment_r3_1, 0, 0];
           %A4_0 = [moment_3_0, moment_3_0, moment_3_0,  0, moment_r3_0, 0, 0];
           AA1_1 = repmat(A1_1, n20(i), 1);
           AA2_1 = repmat(A2_1, n02(i), 1);
           AA3_1 = repmat(A3_1, n11(i), 1);
           %AA4_1 = repmat(A4_1, n3(i), 1);
           AA1_0 = repmat(A1_0, ntot-n20(i), 1);
           AA2_0 = repmat(A2_0, ntot-n02(i), 1);
           AA3_0 = repmat(A3_0, ntot-n11(i), 1);
           %AA4_0 = repmat(A4_0, ntot-n3(i), 1);
           %AA=[AA1_1;AA2_1;AA3_1;AA4_1;AA1_0;AA2_0;AA3_0;AA4_0];
           AA=[AA1_1;AA2_1;AA3_1;AA1_0;AA2_0;AA3_0];
        elseif t(i) == 3
           ntot= n20(i) +n02(i) +n11(i)+n3(i);
           A1_1 = [moment_20_1, 0, 0, 0, 0, moment_r20_1, 0];
           A1_0 = [moment_20_0, 0, 0, 0, 0, moment_r20_0, 0];
           A2_1 = [0, moment_02_1, 0, 0, 0, moment_r02_1, 0];
           A2_0 = [0, moment_02_0, 0, 0, 0, moment_r02_0, 0];
           A3_1 = [0, 0, moment_11_1, 0, 0, moment_r11_1, 0];
           A3_0 = [0, 0, moment_11_0, 0, 0, moment_r11_0, 0];
           %A4_1 = [moment_3_1, moment_3_1, moment_3_1,  0,  0, moment_r3_1, 0];
           %A4_0 = [moment_3_0, moment_3_0, moment_3_0,  0,  0, moment_r3_0, 0];
           AA1_1 = repmat(A1_1, n20(i), 1);
           AA2_1 = repmat(A2_1, n02(i), 1);
           AA3_1 = repmat(A3_1, n11(i), 1);
           %AA4_1 = repmat(A4_1, n3(i), 1);
           AA1_0 = repmat(A1_0, ntot-n20(i), 1);
           AA2_0 = repmat(A2_0, ntot-n02(i), 1);
           AA3_0 = repmat(A3_0, ntot-n11(i), 1);
           %AA4_0 = repmat(A4_0, ntot-n3(i), 1);
           %AA=[AA1_1;AA2_1;AA3_1;AA4_1;AA1_0;AA2_0;AA3_0;AA4_0];
           AA=[AA1_1;AA2_1;AA3_1;AA1_0;AA2_0;AA3_0];
        elseif t(i) == 4
           ntot= n20(i) +n02(i) +n11(i)+n3(i);
           A1_1 = [moment_20_1, 0, 0, 0, 0, 0, moment_r20_1];
           A1_0 = [moment_20_0, 0, 0, 0, 0, 0, moment_r20_0];
           A2_1 = [0, moment_02_1, 0, 0, 0, 0, moment_r02_1];
           A2_0 = [0, moment_02_0, 0, 0, 0, 0, moment_r02_0];
           A3_1 = [0, 0, moment_11_1, 0, 0, 0, moment_r11_1];
           A3_0 = [0, 0, moment_11_0, 0, 0, 0, moment_r11_0];
           %A4_1 = [moment_3_1, moment_3_1, moment_3_1,  0,  0,  0, moment_r3_1];
           %A4_0 = [moment_3_0, moment_3_0, moment_3_0,  0,  0,  0, moment_r3_0];
           AA1_1 = repmat(A1_1, n20(i), 1);
           AA2_1 = repmat(A2_1, n02(i), 1);
           AA3_1 = repmat(A3_1, n11(i), 1);
           %AA4_1 = repmat(A4_1, n3(i), 1);
           AA1_0 = repmat(A1_0, ntot-n20(i), 1);
           AA2_0 = repmat(A2_0, ntot-n02(i), 1);
           AA3_0 = repmat(A3_0, ntot-n11(i), 1);
           %AA4_0 = repmat(A4_0, ntot-n3(i), 1);
           %AA=[AA1_1;AA2_1;AA3_1;AA4_1;AA1_0;AA2_0;AA3_0;AA4_0];
           AA=[AA1_1;AA2_1;AA3_1;AA1_0;AA2_0;AA3_0];
        end
        Score_t=[Score_t;AA];
    end
    BarS=mean(Score_t);
    nn_all=size(Score_t,1);
    Score_t=Score_t- ones(nn_all,1)*BarS;    
    SS=Score_t'*Score_t/nn_all;
end
