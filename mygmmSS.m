function SS7 = mygmmSS(vars_tilde, t, ivf_awt, d20, d02, d11, d3, ...
                              delta_20, delta_02, delta_11, delta_3, ...
                              n20, n02, n11, n3)
    % Function mygmmSS is for mainGMMIVF0410.m
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
    
    
    SS7 = zeros(7,7);  
    ntot_all=0
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

        %  Calculate SS

        g1_20 =   ivf_awt(i)^2 * gamma_t^2;      
        g2_02 =   ivf_awt(i)^2 * gamma_t^2;      
        g3_11 =   ivf_awt(i)^2 * gamma_t^2;      
        
        g1_t_1  = -((1-d20(i))* ivf_awt(i) - ... 
                        2*(d20(i)-theta_20) * ivf_awt(i)* gamma_t);
        g1_t_0  = -((0-d20(i))* ivf_awt(i) - ... 
                        2*(d20(i)-theta_20) * ivf_awt(i)* gamma_t);

        g2_t_1  = -((1-d02(i))* ivf_awt(i) - ... 
                        2*(d02(i)-theta_02) * ivf_awt(i)* gamma_t);
        g2_t_0  = -((0-d02(i))* ivf_awt(i) - ... 
                        2*(d02(i)-theta_02) * ivf_awt(i)* gamma_t);

        g3_t_1  = -((1-d11(i))* ivf_awt(i) - ... 
                        2*(d11(i)-theta_11) * ivf_awt(i)* gamma_t);
        g3_t_0  = -((0-d11(i))* ivf_awt(i) - ... 
                        2*(d11(i)-theta_11) * ivf_awt(i)* gamma_t);
        
        % Add to S(1,1), S(2,2), S(3,3)
        ntot= n20(i) +n02(i) +n11(i)+n3(i);

        SS7(1,1) = SS7(1,1) + g1_20 * ntot;        
        SS7(2,2) = SS7(2,2) + g2_02 * ntot;
        SS7(3,3) = SS7(3,3) + g3_11 * ntot;

        gt_20 = ((d20(i)-theta_20)*ivf_awt(i))^2;
        gt_02 = ((d02(i)-theta_02)*ivf_awt(i))^2;
        gt_11 = ((d11(i)-theta_11)*ivf_awt(i))^2;
        % gt_3  =  ((d3(i) -theta_3) *ivf_awt(i))^2;
        % gt_tt = gt_20 + gt_02 + gt_11 + gt_3;
        gt_tt = gt_20 + gt_02 + gt_11;

        gtp20_1 = -((1-d20(i))*ivf_awt(i)-...
                    2*(d20(i)-theta_20)*ivf_awt(i)^2*gamma_t);
        gtp20_0 = -((0-d20(i))*ivf_awt(i)-...
                    2*(d20(i)-theta_20)*ivf_awt(i)^2*gamma_t);
        gtp02_1 = -((1-d02(i))*ivf_awt(i)-...
                    2*(d02(i)-theta_02)*ivf_awt(i)^2*gamma_t);
        gtp02_0 = -((0-d02(i))*ivf_awt(i)-...
                    2*(d02(i)-theta_02)*ivf_awt(i)^2*gamma_t);        
        gtp11_1 = -((1-d11(i))*ivf_awt(i)-...
                    2*(d11(i)-theta_11)*ivf_awt(i)^2*gamma_t);
         gtp11_0 = -((0-d11(i))*ivf_awt(i)-...
                    2*(d11(i)-theta_11)*ivf_awt(i)^2*gamma_t);

         
        if t(i) == 1     
            SS7(1,4) = SS7(1,4) + g1_t_1* n20(i) + g1_t_0*(ntot-n20(i));
            SS7(2,4) = SS7(2,4) + g2_t_1* n02(i) + g2_t_0*(ntot-n02(i));
            SS7(3,4) = SS7(3,4) + g3_t_1* n11(i) + g3_t_0*(ntot-n11(i));
            SS7(4,4) = SS7(4,4) + gt_tt* ntot;
            SS7(4,1) = SS7(4,1) + gtp20_1*n20(i) + gtp20_0*(ntot-n20(i));
            SS7(4,2) = SS7(4,2) + gtp02_1*n02(i) + gtp02_0*(ntot-n02(i));
            SS7(4,3) = SS7(4,3) + gtp11_1*n11(i) + gtp11_0*(ntot-n11(i));
        elseif t(i) == 2
            SS7(1,5) = SS7(1,5) + g1_t_1* n20(i) + g1_t_0*(ntot-n20(i));
            SS7(2,5) = SS7(2,5) + g2_t_1* n02(i) + g2_t_0*(ntot-n02(i));
            SS7(3,5) = SS7(3,5) + g3_t_1* n11(i) + g3_t_0*(ntot-n11(i));
            SS7(5,5) = SS7(5,5) + gt_tt* ntot;
            SS7(5,1) = SS7(5,1) + gtp20_1*n20(i) + gtp20_0*(ntot-n20(i));
            SS7(5,2) = SS7(5,2) + gtp02_1*n02(i) + gtp02_0*(ntot-n02(i));
            SS7(5,3) = SS7(5,3) + gtp11_1*n11(i) + gtp11_0*(ntot-n11(i));
        elseif t(i) == 3
            SS7(1,6) = SS7(1,6) + g1_t_1* n20(i) + g1_t_0*(ntot-n20(i));
            SS7(2,6) = SS7(2,6) + g2_t_1* n02(i) + g2_t_0*(ntot-n02(i));
            SS7(3,6) = SS7(3,6) + g3_t_1* n11(i) + g3_t_0*(ntot-n11(i));
            SS7(6,6) = SS7(6,6) + gt_tt* ntot;
            SS7(6,1) = SS7(6,1) + gtp20_1*n20(i) + gtp20_0*(ntot-n20(i));
            SS7(6,2) = SS7(6,2) + gtp02_1*n02(i) + gtp02_0*(ntot-n02(i));
            SS7(6,3) = SS7(6,3) + gtp11_1*n11(i) + gtp11_0*(ntot-n11(i));
        elseif t(i) == 4
            SS7(1,7) = SS7(1,7) + g1_t_1* n20(i) + g1_t_0*(ntot-n20(i));
            SS7(2,7) = SS7(2,7) + g2_t_1* n02(i) + g2_t_0*(ntot-n02(i));
            SS7(3,7) = SS7(3,7) + g3_t_1* n11(i) + g3_t_0*(ntot-n11(i));
            SS7(7,7) = SS7(7,7) + gt_tt* ntot;
            SS7(7,1) = SS7(7,1) + gtp20_1*n20(i) + gtp20_0*(ntot-n20(i));
            SS7(7,2) = SS7(7,2) + gtp02_1*n02(i) + gtp02_0*(ntot-n02(i));
            SS7(7,3) = SS7(7,3) + gtp11_1*n11(i) + gtp11_0*(ntot-n11(i));
        end
        ntot_all=ntot_all + ntot * 3;
    end
    SS7=SS7/ntot_all;
end
