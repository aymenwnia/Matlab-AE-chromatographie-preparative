% Simulation ISMB classique : 2 fractions, 4 zones, 4 sous-étapes
clc
clear variables
close all

% Chargement conditions initiales
load_CI = true;                                         % true ou false (à modifier)
nom_fichier_CI_old = 'CI_ISMB_reglage2_2cycles';         % (à modifier)
nom_fichier_CI_new = 'CI_ISMB_reglage2_3cycles';         % (à modifier)

% Nom du fichier pour la sauvegarde des données
nom_fichier = 'donnees_a_supprimer';

% Chargement des points expérimentaux
load_exp = false;                                           % true ou false
nom_fichier_exp = 'Exp_ISMB_CIMV';

% Chargement des paramètres
    % Conditions opératoires et paramètres ISMB
    Ncol_tot = 4;                                           % -
    Ratio_AF_RES = 1;                                       % - (0.001 si CF ou 1 si AF)
    NETcol_AF = 180;                                        % - (à modifier d'après ensemble des résultats !!!)
    Porosite_AF = 0.39;                                     % BVcol (à modifier d'après résultats avec NaCl !!!)
    
    NETcol_CF = 130;                                        % - 
    Porosite_CF = 0.37;                                     % BVcol
    
    Volume_tuyau = 0.05;                                    % BVcol
    NETcol = round(NETcol_CF*(1-Ratio_AF_RES)+NETcol_AF*Ratio_AF_RES);
    Porosite = Porosite_CF*(1-Ratio_AF_RES)+Porosite_AF*Ratio_AF_RES;

    Ncycles_avant_opti = 1;                                 % - 
    Ncycles_opti = 0;                                       % - 
    Ncycles = Ncycles_avant_opti + Ncycles_opti;            % -
    Load_BV_opti = false;                                   % -
    Dernier_BV_opti = true;                                 % -
    
    % Réglage du procédé ISMB : ATTENTION (BV1 et BV3) > BV2 >= BV4 !!!
    BV1_PT = 1.1;                                          % - (à optimiser !!!)
    BV3_PT = 1.2;                                          % - (à optimiser !!!)
    BV2_PT = 0.7;                                          % - (à optimiser !!!)
    BV4_PT = 0.5;                                          % - (à optimiser !!!)
    
    % Paramètres pour l'optimisation
    % 1 = sucres, 2 = acétique, 3 = formique, 4 = sels 
    Rapides = [1,4];                                        % -
    Lents = [2,3];                                          % - 
    
    Perte_lent_Z3_droite_limite = 2;                        % %
    Perte_lent_Z3_droite_cible = 1;                         % %
    Perte_lent_Z4_droite_limite = 2;                        % %
    Perte_lent_Z4_droite_cible = 1;                         % %
    Perte_rapide_Z1_gauche_limite = 5;                      % %
    Perte_rapide_Z1_gauche_cible = 2;                       % %
    Perte_rapide_Z1_droite_limite = 5;                      % %
    Perte_rapide_Z1_droite_cible = 2;                       % %
    
    Taux_penetration_Z4_lent_cible = 25;                    % %
    Taux_penetration_Z1_gauche_rapide_cible = 25;           % %
    Taux_penetration_Z1_droite_rapide_cible = 25;           % %
    DeltaV_est_max = 0.02;                                  % BV
    DeltaBV_max = 2.5;                                      % %  
    
    if Load_BV_opti == true
        load(nom_fichier_CI_old)
        if Dernier_BV_opti == true
            BV1_PT = BV1_fig3(end);                             % -
            BV3_PT = BV3_fig3(end);                             % -
            BV2_PT = BV2_fig3(end);                             % -
            BV4_PT = BV4_fig3(end);                             % -
        
        else
            BV1_PT = BV1_fig3(end-1);                             % -
            BV3_PT = BV3_fig3(end-1);                             % -
            BV2_PT = BV2_fig3(end-1);                             % -
            BV4_PT = BV4_fig3(end-1);                             % -
        end        
    end  
    
    V_ech_est = BV3_PT-BV2_PT;                              % BV
    V_elu1_est = BV1_PT-BV2_PT;                             % BV
    V_elu2_est = BV2_PT-BV4_PT;                             % BV
    V_boucle_est = BV4_PT-Volume_tuyau;                     % BV

    V_etape_est = V_ech_est+V_elu1_est+V_elu2_est+V_boucle_est; % BV
    V_cycle_est = Ncol_tot*V_etape_est;                     % BV
    V_total_est = Ncycles*V_cycle_est;                      % BV
    
    NETtot = NETcol*Ncol_tot;                               % - 
    NET_SE1 = NETcol;                                       % -
    NET_SE2 = NETcol;                                       % -
    NET_SE3 = NETcol*3;                                     % -
    NET_SE4 = NETcol*Ncol_tot;                              % -
    
    % Paramètres résine CF-Af-Cf 
    qmax_s_CF = 20;                                         % eq/L_rés
    Ks_A1H_mat_CF = 0.024;                                  % L_sol/mol
    Ks_A2H_mat_CF = 0.036;                                  % L_sol/mol
    Ks_A3H_mat_CF = 0.040;                                  % L_sol/mol
    N_mat_A1H_CF = 1;                                       % -  
    N_mat_A2H_CF = 1;                                       % -  
    N_mat_A3H_CF = 1;                                       % -   

    % Paramètres résine AF
    qmax_ei_AF = 2.4;                                       % eq/L_rés
    Keq_HSO4_OH = 85;                                       % -
    Keq_SO4_OH = 150;                                       % L_sol/L_rés
    
    Keq_A1_OH = 0;                                          % -
    Keq_A2_OH = 5;                                          % -
    Keq_A3_OH = 0;                                          % -
    % valeur Keq pour info : glu = xyl = ara = 0, form = 3, acet = 5, lact= 7 
    
    Ks_A1H = 0;                                            % L_sol/mol
    Ks_A2H = 0.215;                                        % L_sol/mol (à modifier d'après résultats avec acide acétique !!!)
    Ks_A3H = 0;                                            % L_sol/mol
    % valeur Ks contre-anions(= 0 pour glu, xyl, ara, mais > 0 pour form, acét et lact)
    
    Ks_A1H_HSO4 = Ks_A1H;                                    % L_sol/mol
    Ks_A2H_HSO4 = Ks_A2H;                                    % L_sol/mol
    Ks_A3H_HSO4 = Ks_A3H;                                    % L_sol/mol
    Ks_A1H_SO4 = Ks_A1H;                                     % L_sol/mol
    Ks_A2H_SO4 = Ks_A2H;                                     % L_sol/mol 
    Ks_A3H_SO4 = Ks_A3H;                                     % L_sol/mol
    
    N_HSO4_A1H = 1;                                         % -
    N_HSO4_A2H = 1/4;                                       % - 
    N_HSO4_A3H = 1;                                         % -  
    N_SO4_A1H = 1;                                          % - 
    N_SO4_A2H = 1/4;                                        % -
    N_SO4_A3H = 1;                                          % -
    % valeur N pour info : glu = xyl = ara = 1, form = 1/8, acet = 1/4, lact = 1/4
    
    qmax_s_AF = 20;                                         % eq/L_rés
    Ks_A1H_mat_AF = 0.000;                                  % L_sol/mol
    Ks_A2H_mat_AF = 0.000;                                  % L_sol/mol
    Ks_A3H_mat_AF = 0.0185;                                 % L_sol/mol (à modifier d'après résultats avec glucose !!!)
    N_mat_A1H_AF = 1;                                       % -    
    N_mat_A2H_AF = 1;                                       % -    
    N_mat_A3H_AF = 1;                                       % -    
    % valeur Ks (> 0 pour glu, xyl, ara, mais = 0 pour form, acét et lact)
    
    % Paramètres moyens résines
    qmax_ei = qmax_ei_AF*Ratio_AF_RES;
    qmax_s = qmax_s_AF*Ratio_AF_RES+qmax_s_CF*(1-Ratio_AF_RES);
    Ks_A1H_mat = (qmax_s_AF*Ratio_AF_RES*Ks_A1H_mat_AF+qmax_s_CF*(1-Ratio_AF_RES)*Ks_A1H_mat_CF)/qmax_s;
    Ks_A2H_mat = (qmax_s_AF*Ratio_AF_RES*Ks_A2H_mat_AF+qmax_s_CF*(1-Ratio_AF_RES)*Ks_A2H_mat_CF)/qmax_s;
    Ks_A3H_mat = (qmax_s_AF*Ratio_AF_RES*Ks_A3H_mat_AF+qmax_s_CF*(1-Ratio_AF_RES)*Ks_A3H_mat_CF)/qmax_s;
    N_mat_A1H = N_mat_A1H_AF*Ratio_AF_RES+N_mat_A1H_CF*(1-Ratio_AF_RES);    
    N_mat_A2H = N_mat_A2H_AF*Ratio_AF_RES+N_mat_A2H_CF*(1-Ratio_AF_RES);
    N_mat_A3H = N_mat_A3H_AF*Ratio_AF_RES+N_mat_A3H_CF*(1-Ratio_AF_RES);

    % Paramètres solutés
    pKa_HSO4_0 = 1.99;                                      % -
    pKa_A1H_0 = 14;                                         % -   
    pKa_A2H_0 = 4.76;                                       % -   
    pKa_A3H_0 = 14;                                         % - 
    % valeur pKa pour info : glu = xyl = ara = 14, form = 3.75, acet = 4.76, lact = 3.86
    
    Ka_HSO4_0 = 10^(-pKa_HSO4_0);                           % mol/L_sol
    Ka_A1H_0 = 10^(-pKa_A1H_0);                             % mol/L_sol
    Ka_A2H_0 = 10^(-pKa_A2H_0);                             % mol/L_sol
    Ka_A3H_0 = 10^(-pKa_A3H_0);                             % mol/L_sol
    Keq_HSO4_SO4 = Keq_HSO4_OH^2/Keq_SO4_OH;                % L_sol/L_rés
    Keq_A1_HSO4 = Keq_A1_OH/Keq_HSO4_OH;                    % -    
    Keq_A1_SO4 = Keq_HSO4_SO4*Keq_A1_HSO4^2;                % L_sol/L_rés
    Keq_A2_HSO4 = Keq_A2_OH/Keq_HSO4_OH;                    % -    
    Keq_A2_SO4 = Keq_HSO4_SO4*Keq_A2_HSO4^2;                % L_sol/L_rés    
    Keq_A3_HSO4 = Keq_A3_OH/Keq_HSO4_OH;                    % -    
    Keq_A3_SO4 = Keq_HSO4_SO4*Keq_A3_HSO4^2;                % L_sol/L_rés 
    N_HSO4_min = min([N_HSO4_A1H,N_HSO4_A2H,N_HSO4_A3H]);   % -
    N_SO4_min = min([N_SO4_A1H,N_SO4_A2H,N_SO4_A3H]);       % -
    N_mat_min = min([N_mat_A1H,N_mat_A2H,N_mat_A3H]);       % -
    
    % Composition du produit
    pH_ech = 2.00;                                          % -          (à vérifier !!!)
    A1tot_ech = 1.711;                                      % mol/L_sol  (à modifier !!!)
    A2tot_ech = 1.667;                                      % mol/L_sol  (à modifier !!!)
    A3tot_ech = 0.556;                                      % mol/L_sol  (à modifier !!!)
    SO4tot_ech = 0.02;                                      % mol/L_sol  (à vérifier !!!)
    H_ech = 10^(-pH_ech);                                   % mol/L_sol 
    Ka_A1H_ech = Ka_A1H_0;                                  % mol/L_sol
    Ka_A2H_ech = Ka_A2H_0;                                  % mol/L_sol
    Ka_A3H_ech = Ka_A3H_0;                                  % mol/L_sol
    Ka_HSO4_ech = Ka_HSO4_0;                                % mol/L_sol
    Variation_negligeable = 0;
    while Variation_negligeable == 0   
        A1H_ech = A1tot_ech/(1+Ka_A1H_ech/H_ech);           % mol/L_sol     
        A1_ech = A1H_ech*Ka_A1H_ech/H_ech;                  % mol/L_sol 
        A2H_ech = A2tot_ech/(1+Ka_A2H_ech/H_ech);           % mol/L_sol     
        A2_ech = A2H_ech*Ka_A2H_ech/H_ech;                  % mol/L_sol        
        A3H_ech = A3tot_ech/(1+Ka_A3H_ech/H_ech);           % mol/L_sol     
        A3_ech = A3H_ech*Ka_A3H_ech/H_ech;                  % mol/L_sol              
        HSO4_ech = SO4tot_ech/(1+Ka_HSO4_ech/H_ech);        % mol/L_sol
        SO4_ech = HSO4_ech*Ka_HSO4_ech/H_ech;               % mol/L_sol
        F_ionique_ech = 0.5*(A1_ech+A2_ech+A3_ech+HSO4_ech+4*SO4_ech);
        Activite_AH_ech = 10^(0.51*F_ionique_ech^0.5/(1+F_ionique_ech^0.5));
        Activite_HSO4_ech = 10^(0.51*3*F_ionique_ech^0.5/(1+F_ionique_ech^0.5));
        Ka_A1H_ech_new = Ka_A1H_0*Activite_AH_ech;
        Ka_A2H_ech_new = Ka_A2H_0*Activite_AH_ech;
        Ka_A3H_ech_new = Ka_A3H_0*Activite_AH_ech;
        Ka_HSO4_ech_new = Ka_HSO4_0*Activite_HSO4_ech;
        Variation_negligeable = (Ka_A1H_ech_new-Ka_A1H_ech)/Ka_A1H_ech < 0.001 && (Ka_A2H_ech_new-Ka_A2H_ech)/Ka_A2H_ech < 0.001 && (Ka_A3H_ech_new-Ka_A3H_ech)/Ka_A3H_ech < 0.001 && (Ka_HSO4_ech_new-Ka_HSO4_ech)/Ka_HSO4_ech < 0.001;
        Ka_A1H_ech = Ka_A1H_ech_new;                         % mol/L_sol
        Ka_A2H_ech = Ka_A2H_ech_new;                         % mol/L_sol
        Ka_A3H_ech = Ka_A3H_ech_new;                         % mol/L_sol
        Ka_HSO4_ech = Ka_HSO4_ech_new;                       % mol/L_sol        
    end

    % Composition de l'éluant
    pH_elu = 2.00;                                          %           (à vérifier !!!)
    A1tot_elu = 0;                                          % mol/L_sol 
    A2tot_elu = 0;                                          % mol/L_sol    
    A3tot_elu = 0;                                          % mol/L_sol   
    SO4tot_elu = 0.02;                                      % mol/L_sol (à vérifier !!!)
    H_elu = 10^(-pH_elu);                                   % mol/L_sol
    Ka_A1H_elu = Ka_A1H_0;                                  % mol/L_sol
    Ka_A2H_elu = Ka_A2H_0;                                  % mol/L_sol
    Ka_A3H_elu = Ka_A3H_0;                                  % mol/L_sol
    Ka_HSO4_elu = Ka_HSO4_0;                                % mol/L_sol
    Variation_negligeable = 0;
    while Variation_negligeable == 0   
        A1H_elu = A1tot_elu/(1+Ka_A1H_elu/H_elu);           % mol/L_sol     
        A1_elu = A1H_elu*Ka_A1H_elu/H_elu;                  % mol/L_sol 
        A2H_elu = A2tot_elu/(1+Ka_A2H_elu/H_elu);           % mol/L_sol     
        A2_elu = A2H_elu*Ka_A2H_elu/H_elu;                  % mol/L_sol        
        A3H_elu = A3tot_elu/(1+Ka_A3H_elu/H_elu);           % mol/L_sol     
        A3_elu = A3H_elu*Ka_A3H_elu/H_elu;                  % mol/L_sol              
        HSO4_elu = SO4tot_elu/(1+Ka_HSO4_elu/H_elu);        % mol/L_sol
        SO4_elu = HSO4_elu*Ka_HSO4_elu/H_elu;               % mol/L_sol
        F_ionique_elu = 0.5*(A1_elu+A2_elu+A3_elu+HSO4_elu+4*SO4_elu);
        Activite_AH_elu = 10^(0.51*F_ionique_elu^0.5/(1+F_ionique_elu^0.5));
        Activite_HSO4_elu = 10^(0.51*3*F_ionique_elu^0.5/(1+F_ionique_elu^0.5));
        Ka_A1H_elu_new = Ka_A1H_0*Activite_AH_elu;
        Ka_A2H_elu_new = Ka_A2H_0*Activite_AH_elu;
        Ka_A3H_elu_new = Ka_A3H_0*Activite_AH_elu;
        Ka_HSO4_elu_new = Ka_HSO4_0*Activite_HSO4_elu;
        Variation_negligeable = (Ka_A1H_elu_new-Ka_A1H_elu)/Ka_A1H_elu < 0.001 && (Ka_A2H_elu_new-Ka_A2H_elu)/Ka_A2H_elu < 0.001 && (Ka_A3H_elu_new-Ka_A3H_elu)/Ka_A3H_elu < 0.001 && (Ka_HSO4_elu_new-Ka_HSO4_elu)/Ka_HSO4_elu < 0.001;
        Ka_A1H_elu = Ka_A1H_elu_new;                         % mol/L_sol
        Ka_A2H_elu = Ka_A2H_elu_new;                         % mol/L_sol
        Ka_A3H_elu = Ka_A3H_elu_new;                         % mol/L_sol
        Ka_HSO4_elu = Ka_HSO4_elu_new;                       % mol/L_sol        
    end 

    % Conditions initiales en phase mobile
    pH_ini = 2.00;                                          % -          (à vérifier !!!)
    A1tot_mob_ini = 0.00;                                   % mol/L_sol 
    A2tot_mob_ini = 0.00;                                   % mol/L_sol 
    A3tot_mob_ini = 0.00;                                   % mol/L_sol  
    SO4tot_mob_ini = 0.02;                                  % mol/L_sol  (à vérifier !!!)
    H_mob_ini = 10^(-pH_ini);                               % mol/L_sol
    Ka_A1H_ini = Ka_A1H_0;                                  % mol/L_sol
    Ka_A2H_ini = Ka_A2H_0;                                  % mol/L_sol
    Ka_A3H_ini = Ka_A3H_0;                                  % mol/L_sol
    Ka_HSO4_ini = Ka_HSO4_0;                                % mol/L_sol
    Variation_negligeable = 0;
    while Variation_negligeable == 0   
        A1H_mob_ini = A1tot_mob_ini/(1+Ka_A1H_ini/H_mob_ini);           % mol/L_sol     
        A1_mob_ini = A1H_mob_ini*Ka_A1H_ini/H_mob_ini;                  % mol/L_sol 
        A2H_mob_ini = A2tot_mob_ini/(1+Ka_A2H_ini/H_mob_ini);           % mol/L_sol     
        A2_mob_ini = A2H_mob_ini*Ka_A2H_ini/H_mob_ini;                  % mol/L_sol        
        A3H_mob_ini = A3tot_mob_ini/(1+Ka_A3H_ini/H_mob_ini);           % mol/L_sol     
        A3_mob_ini = A3H_mob_ini*Ka_A3H_ini/H_mob_ini;                  % mol/L_sol              
        HSO4_mob_ini = SO4tot_mob_ini/(1+Ka_HSO4_ini/H_mob_ini);        % mol/L_sol
        SO4_mob_ini = HSO4_mob_ini*Ka_HSO4_ini/H_mob_ini;               % mol/L_sol
        F_ionique_ini = 0.5*(A1_mob_ini+A2_mob_ini+A3_mob_ini+HSO4_mob_ini+4*SO4_mob_ini);
        Activite_AH_ini = 10^(0.51*F_ionique_ini^0.5/(1+F_ionique_ini^0.5));
        Activite_HSO4_ini = 10^(0.51*3*F_ionique_ini^0.5/(1+F_ionique_ini^0.5));
        Ka_A1H_ini_new = Ka_A1H_0*Activite_AH_ini;
        Ka_A2H_ini_new = Ka_A2H_0*Activite_AH_ini;
        Ka_A3H_ini_new = Ka_A3H_0*Activite_AH_ini;
        Ka_HSO4_ini_new = Ka_HSO4_0*Activite_HSO4_ini;
        Variation_negligeable = (Ka_A1H_ini_new-Ka_A1H_ini)/Ka_A1H_ini < 0.001 && (Ka_A2H_ini_new-Ka_A2H_ini)/Ka_A2H_ini < 0.001 && (Ka_A3H_ini_new-Ka_A3H_ini)/Ka_A3H_ini < 0.001 && (Ka_HSO4_ini_new-Ka_HSO4_ini)/Ka_HSO4_ini < 0.001;
        Ka_A1H_ini = Ka_A1H_ini_new;                         % mol/L_sol
        Ka_A2H_ini = Ka_A2H_ini_new;                         % mol/L_sol
        Ka_A3H_ini = Ka_A3H_ini_new;                         % mol/L_sol
        Ka_HSO4_ini = Ka_HSO4_ini_new;                       % mol/L_sol        
    end    

    % Conditions initiales en phase stationnaire
        % Calculs intermédiaires
        A_HSO4_stat = 2*Ka_HSO4_ini/Keq_HSO4_SO4;
        B_HSO4_stat = H_mob_ini*HSO4_mob_ini+Ka_A1H_ini*Keq_A1_HSO4*A1H_mob_ini+Ka_A2H_ini*Keq_A2_HSO4*A2H_mob_ini+Ka_A3H_ini*Keq_A3_HSO4*A3H_mob_ini;
        C_HSO4_stat = -H_mob_ini*HSO4_mob_ini*qmax_ei;
        D_HSO4_stat = B_HSO4_stat^2-4*A_HSO4_stat*C_HSO4_stat;
    HSO4_stat_ini = (sqrt(D_HSO4_stat)-B_HSO4_stat)/(2*A_HSO4_stat);                                % mol/L_res    
    SO4_stat_ini = HSO4_stat_ini^2*SO4_mob_ini/(HSO4_mob_ini^2*Keq_HSO4_SO4);                       % mol/L_res
    SO4tot_stat_ini = SO4_stat_ini+HSO4_stat_ini;                                                   % mol/L_res
    Site_HSO4_libre_ini = HSO4_stat_ini/(N_HSO4_min+N_HSO4_A1H*Ks_A1H_HSO4*A1H_mob_ini+N_HSO4_A2H*Ks_A2H_HSO4*A2H_mob_ini+N_HSO4_A3H*Ks_A3H_HSO4*A3H_mob_ini);  % mol/L_res
    Site_SO4_libre_ini = SO4_stat_ini/(N_SO4_min+N_SO4_A1H*Ks_A1H_SO4*A1H_mob_ini+N_SO4_A2H*Ks_A2H_SO4*A2H_mob_ini+N_SO4_A3H*Ks_A3H_SO4*A3H_mob_ini);           % mol/L_res
    Site_mat_libre_ini = qmax_s/(N_mat_min+N_mat_A1H*Ks_A1H_mat*A1H_mob_ini+N_mat_A2H*Ks_A2H_mat*A2H_mob_ini+N_mat_A3H*Ks_A3H_mat*A3H_mob_ini);                 % mol/L_res
    A1H_HSO4_stat_ini = Ks_A1H_HSO4*A1H_mob_ini*Site_HSO4_libre_ini;                                % mol/L_res
    A1H_SO4_stat_ini = Ks_A1H_SO4*A1H_mob_ini*Site_SO4_libre_ini;                                   % mol/L_res
    A1H_mat_stat_ini = Ks_A1H_mat*A1H_mob_ini*Site_mat_libre_ini;                                   % mol/L_res
    A1_stat_ini = Keq_A1_HSO4*HSO4_stat_ini*A1_mob_ini/HSO4_mob_ini;                                % mol/L_res
    A1tot_stat_ini =A1H_HSO4_stat_ini+A1H_SO4_stat_ini+A1H_mat_stat_ini+A1_stat_ini;                % mol/L_res
    A2H_HSO4_stat_ini = Ks_A2H_HSO4*A2H_mob_ini*Site_HSO4_libre_ini;                                % mol/L_res
    A2H_SO4_stat_ini = Ks_A2H_SO4*A2H_mob_ini*Site_SO4_libre_ini;                                   % mol/L_res
    A2H_mat_stat_ini = Ks_A2H_mat*A2H_mob_ini*Site_mat_libre_ini;                                   % mol/L_res
    A2_stat_ini = Keq_A2_HSO4*HSO4_stat_ini*A2_mob_ini/HSO4_mob_ini;                                % mol/L_res
    A2tot_stat_ini =A2H_HSO4_stat_ini+A2H_SO4_stat_ini+A2H_mat_stat_ini+A2_stat_ini;                % mol/L_res   
    A3H_HSO4_stat_ini = Ks_A3H_HSO4*A3H_mob_ini*Site_HSO4_libre_ini;                                % mol/L_res
    A3H_SO4_stat_ini = Ks_A3H_SO4*A3H_mob_ini*Site_SO4_libre_ini;                                   % mol/L_res
    A3H_mat_stat_ini = Ks_A3H_mat*A3H_mob_ini*Site_mat_libre_ini;                                   % mol/L_res
    A3_stat_ini = Keq_A3_HSO4*HSO4_stat_ini*A3_mob_ini/HSO4_mob_ini;                                % mol/L_res
    A3tot_stat_ini =A3H_HSO4_stat_ini+A3H_SO4_stat_ini+A3H_mat_stat_ini+A3_stat_ini;                % mol/L_res
    
    % Paramètres numériques
    Pas_V_calcul = 0.25*Porosite/NETcol;                    % BV
    Pas_V_figure = 0.1;                                     % BV (à modifier si besoin)
    Pas_V_memoire = 0.10;                                   % BV
    Pas_n_memoire = 4;                                      % -
    N_iter_max_Newton = 25;                                 % -
    N_iter_max_F = 10;                                      % -
    N_iter_max_pas_V = 3;                                   % -
    Facteur_diminution_pas_V = 2;                           % -    
    Facteur_freinage = 2;                                   % -
    Facteur_acceleration = 1.5;                             % -
    Variation_rel_max = 3;                                  % -
    Precision_Conc_abs = 1e-7;                              % mol/L
    Precision_Conc_rel = 0.01/100;                          % %
    Precision_Bilan_abs = 1e-7/NETcol;                      % mol/L
    Precision_Bilan_rel = 0.01/100;                         % %
    Methode_initialisation = 3;                             % -
    % 1 = composition à l'instant t-1
    % 2 = composition à l'instant t si aucune réaction
    % 3 = composition en acides organiques à l'instant t en considérant [H+], [HSO4-] et qHSO4- à l'instant t-1
    
% Calculs intermédiaires
C1 = (1-Porosite)/Porosite;
C2 = Ka_HSO4_0/Keq_HSO4_SO4;
C3_1 = Ka_A1H_0*Keq_A1_HSO4;
C3_2 = Ka_A2H_0*Keq_A2_HSO4;
C3_3 = Ka_A3H_0*Keq_A3_HSO4;
C6_1 = N_SO4_A1H*Ks_A1H_SO4;
C6_2 = N_SO4_A2H*Ks_A2H_SO4;
C6_3 = N_SO4_A3H*Ks_A3H_SO4;
C7_1 = N_HSO4_A1H*Ks_A1H_HSO4;
C7_2 = N_HSO4_A2H*Ks_A2H_HSO4;
C7_3 = N_HSO4_A3H*Ks_A3H_HSO4;
C38_1 = N_mat_A1H*Ks_A1H_mat;
C38_2 = N_mat_A2H*Ks_A2H_mat;
C38_3 = N_mat_A3H*Ks_A3H_mat;
C8_1 = Ks_A1H_mat*qmax_s;
C8_2 = Ks_A2H_mat*qmax_s;
C8_3 = Ks_A3H_mat*qmax_s;
C9 = 2*C2/qmax_ei;
C10 = 1/qmax_ei;
C11_1 = C3_1/qmax_ei;
C11_2 = C3_2/qmax_ei;
C11_3 = C3_3/qmax_ei;
C12 = 2*C9;
C15 = C1*C2;
C16_1 = C15*Ks_A1H_SO4;
C16_2 = C15*Ks_A2H_SO4;
C16_3 = C15*Ks_A3H_SO4;
C17_1 = C1*Ks_A1H_HSO4;
C17_2 = C1*Ks_A2H_HSO4;
C17_3 = C1*Ks_A3H_HSO4;
C18_1 = C1*C8_1;
C18_2 = C1*C8_2;
C18_3 = C1*C8_3;
C19_1 = C1*C3_1;
C19_2 = C1*C3_2;
C19_3 = C1*C3_3;

Entree_Z1 = zeros(Ncol_tot,1);
Entree_Z2 = zeros(Ncol_tot,1);
Entree_Z3 = zeros(Ncol_tot,1);
Entree_Z4 = zeros(Ncol_tot,1);
Sortie_Z1 = zeros(Ncol_tot,1);
Sortie_Z2 = zeros(Ncol_tot,1);
Sortie_Z3 = zeros(Ncol_tot,1);
Sortie_Z4 = zeros(Ncol_tot,1);
Milieu_Z1 = zeros(Ncol_tot,1);
Milieu_Z2 = zeros(Ncol_tot,1);
Milieu_Z3 = zeros(Ncol_tot,1);
Milieu_Z4 = zeros(Ncol_tot,1);

for i=1:Ncol_tot
    Entree_Z1(i,1) = 1+(i-1)*NETcol;
    Entree_Z2(i,1) = (1+(1+(i-1))*NETcol)-NETtot*floor((1+(1+(i-1))*NETcol)/(NETtot+1));
    Entree_Z3(i,1) = (1+(2+(i-1))*NETcol)-NETtot*floor((1+(2+(i-1))*NETcol)/(NETtot+1));
    Entree_Z4(i,1) = (1+(3+(i-1))*NETcol)-NETtot*floor((1+(3+(i-1))*NETcol)/(NETtot+1));
    Sortie_Z1(i,1) = (1+(i-1))*NETcol-NETtot*floor((1+(i-1))*NETcol/(NETtot+1));
    Sortie_Z2(i,1) = (2+(i-1))*NETcol-NETtot*floor((2+(i-1))*NETcol/(NETtot+1));
    Sortie_Z3(i,1) = (3+(i-1))*NETcol-NETtot*floor((3+(i-1))*NETcol/(NETtot+1));
    Sortie_Z4(i,1) = (Ncol_tot+(i-1))*NETcol-NETtot*floor((Ncol_tot+(i-1))*NETcol/(NETtot+1));
    Milieu_Z1(i,1) = Entree_Z1(i,1)+NETcol/2;
    Milieu_Z2(i,1) = Entree_Z2(i,1)+NETcol/2;
    Milieu_Z3(i,1) = Entree_Z3(i,1)+NETcol/2;
    Milieu_Z4(i,1) = Entree_Z4(i,1)+NETcol/2;  
end

Correction_F = ones(Ncycles,1);
F_V_ech = ones(Ncycles,1);
F_V_elu1 = ones(Ncycles,1);
F_V_elu2 = ones(Ncycles,1);
F_V_boucle = ones(Ncycles,1);

% Initialisation des variables mémorisées
N_volume = floor(V_total_est/Pas_V_memoire);
N_etages = floor(NETtot/Pas_n_memoire);

Cycle = zeros(N_volume,1);
Etape = zeros(N_volume,1);
Sous_Etape = zeros(N_volume,1);

V_inj_total = zeros(N_volume,1);
V_inj_cycle = zeros(N_volume,1);
V_inj_etape = zeros(N_volume,1);
V_inj_sous_etape = zeros(N_volume,1);
V_ech_total = zeros(N_volume,1);
V_ech_cycle = zeros(N_volume,1);
V_ech_etape = zeros(N_volume,1);
V_ech_sous_etape = zeros(N_volume,1);
V_elu_total = zeros(N_volume,1);
V_elu_cycle = zeros(N_volume,1);
V_elu_etape = zeros(N_volume,1);
V_elu_sous_etape = zeros(N_volume,1);
V_ext_total = zeros(N_volume,1);
V_ext_cycle = zeros(N_volume,1);
V_ext_etape = zeros(N_volume,1);
V_ext_sous_etape = zeros(N_volume,1);
V_raff_total = zeros(N_volume,1);
V_raff_cycle = zeros(N_volume,1);
V_raff_etape = zeros(N_volume,1);
V_raff_sous_etape = zeros(N_volume,1);

H_mob = zeros(N_volume,N_etages);
HSO4_mob = zeros(N_volume,N_etages);
SO4_mob = zeros(N_volume,N_etages);
SO4tot_mob = zeros(N_volume,N_etages);
A1H_mob = zeros(N_volume,N_etages);
A2H_mob = zeros(N_volume,N_etages);
A3H_mob = zeros(N_volume,N_etages);
A1_mob = zeros(N_volume,N_etages);
A2_mob = zeros(N_volume,N_etages);
A3_mob = zeros(N_volume,N_etages);
A1tot_mob = zeros(N_volume,N_etages);
A2tot_mob = zeros(N_volume,N_etages);
A3tot_mob = zeros(N_volume,N_etages);

HSO4_stat = zeros(N_volume,N_etages);
SO4_stat = zeros(N_volume,N_etages);
SO4tot_stat = zeros(N_volume,N_etages);
A1H_HSO4_stat = zeros(N_volume,N_etages);
A2H_HSO4_stat = zeros(N_volume,N_etages);
A3H_HSO4_stat = zeros(N_volume,N_etages);
A1H_SO4_stat = zeros(N_volume,N_etages);
A2H_SO4_stat = zeros(N_volume,N_etages);
A3H_SO4_stat = zeros(N_volume,N_etages);
A1H_mat_stat = zeros(N_volume,N_etages);
A2H_mat_stat = zeros(N_volume,N_etages);
A3H_mat_stat = zeros(N_volume,N_etages);
A1_stat = zeros(N_volume,N_etages);
A2_stat = zeros(N_volume,N_etages);
A3_stat = zeros(N_volume,N_etages);
A1tot_stat = zeros(N_volume,N_etages);
A2tot_stat = zeros(N_volume,N_etages);
A3tot_stat = zeros(N_volume,N_etages);

F_ionique = zeros(N_volume,N_etages);
Activite_AH = zeros(N_volume,N_etages);
Activite_HSO4 = zeros(N_volume,N_etages);

Vecteur_new1 = 1:Pas_n_memoire:NETtot;
Vecteur_new2 = 1:NETtot;
if load_CI == true
    load(nom_fichier_CI_old)
    NETtot_load = length(H_mob_load);
    Vecteur_old = linspace (1,NETtot,NETtot_load);
    
    H_mob(1,:) = pchip(Vecteur_old,H_mob_load,Vecteur_new1);
    HSO4_mob(1,:) = pchip(Vecteur_old,HSO4_mob_load,Vecteur_new1);
    SO4_mob(1,:) = pchip(Vecteur_old,SO4_mob_load,Vecteur_new1);
    SO4tot_mob(1,:) = pchip(Vecteur_old,SO4tot_mob_load,Vecteur_new1);
    A1H_mob(1,:) = pchip(Vecteur_old,A1H_mob_load,Vecteur_new1);
    A2H_mob(1,:) = pchip(Vecteur_old,A2H_mob_load,Vecteur_new1);
    A3H_mob(1,:) = pchip(Vecteur_old,A3H_mob_load,Vecteur_new1);
    A1_mob(1,:) = pchip(Vecteur_old,A1_mob_load,Vecteur_new1);
    A2_mob(1,:) = pchip(Vecteur_old,A2_mob_load,Vecteur_new1);
    A3_mob(1,:) = pchip(Vecteur_old,A3_mob_load,Vecteur_new1);    
    A1tot_mob(1,:) = pchip(Vecteur_old,A1tot_mob_load,Vecteur_new1);
    A2tot_mob(1,:) = pchip(Vecteur_old,A2tot_mob_load,Vecteur_new1);
    A3tot_mob(1,:) = pchip(Vecteur_old,A3tot_mob_load,Vecteur_new1);
    
    H_mob_new = pchip(Vecteur_old,H_mob_load,Vecteur_new2);
    HSO4_mob_new = pchip(Vecteur_old,HSO4_mob_load,Vecteur_new2);
    SO4_mob_new = pchip(Vecteur_old,SO4_mob_load,Vecteur_new2);
    SO4tot_mob_new = pchip(Vecteur_old,SO4tot_mob_load,Vecteur_new2);
    A1H_mob_new = pchip(Vecteur_old,A1H_mob_load,Vecteur_new2);
    A2H_mob_new = pchip(Vecteur_old,A2H_mob_load,Vecteur_new2);
    A3H_mob_new = pchip(Vecteur_old,A3H_mob_load,Vecteur_new2);
    A1_mob_new = pchip(Vecteur_old,A1_mob_load,Vecteur_new2);
    A2_mob_new = pchip(Vecteur_old,A2_mob_load,Vecteur_new2);
    A3_mob_new = pchip(Vecteur_old,A3_mob_load,Vecteur_new2);
    A1tot_mob_new = pchip(Vecteur_old,A1tot_mob_load,Vecteur_new2);
    A2tot_mob_new = pchip(Vecteur_old,A2tot_mob_load,Vecteur_new2);
    A3tot_mob_new = pchip(Vecteur_old,A3tot_mob_load,Vecteur_new2);
    
    H_mob_old = pchip(Vecteur_old,H_mob_load,Vecteur_new2);
    HSO4_mob_old = pchip(Vecteur_old,HSO4_mob_load,Vecteur_new2);
    SO4_mob_old = pchip(Vecteur_old,SO4_mob_load,Vecteur_new2);
    SO4tot_mob_old = pchip(Vecteur_old,SO4tot_mob_load,Vecteur_new2);
    A1H_mob_old = pchip(Vecteur_old,A1H_mob_load,Vecteur_new2);
    A2H_mob_old = pchip(Vecteur_old,A2H_mob_load,Vecteur_new2);
    A3H_mob_old = pchip(Vecteur_old,A3H_mob_load,Vecteur_new2);
    A1_mob_old = pchip(Vecteur_old,A1_mob_load,Vecteur_new2);
    A2_mob_old = pchip(Vecteur_old,A2_mob_load,Vecteur_new2);
    A3_mob_old = pchip(Vecteur_old,A3_mob_load,Vecteur_new2);
    A1tot_mob_old = pchip(Vecteur_old,A1tot_mob_load,Vecteur_new2);
    A2tot_mob_old = pchip(Vecteur_old,A2tot_mob_load,Vecteur_new2);
    A3tot_mob_old = pchip(Vecteur_old,A3tot_mob_load,Vecteur_new2);
    
    HSO4_stat(1,:) = pchip(Vecteur_old,HSO4_stat_load,Vecteur_new1);
    SO4_stat(1,:) = pchip(Vecteur_old,SO4_stat_load,Vecteur_new1);
    SO4tot_stat(1,:) = pchip(Vecteur_old,SO4tot_stat_load,Vecteur_new1);
    A1H_HSO4_stat(1,:) = pchip(Vecteur_old,A1H_HSO4_stat_load,Vecteur_new1);
    A2H_HSO4_stat(1,:) = pchip(Vecteur_old,A2H_HSO4_stat_load,Vecteur_new1);
    A3H_HSO4_stat(1,:) = pchip(Vecteur_old,A3H_HSO4_stat_load,Vecteur_new1);
    A1H_SO4_stat(1,:) = pchip(Vecteur_old,A1H_SO4_stat_load,Vecteur_new1);
    A2H_SO4_stat(1,:) = pchip(Vecteur_old,A2H_SO4_stat_load,Vecteur_new1);
    A3H_SO4_stat(1,:) = pchip(Vecteur_old,A3H_SO4_stat_load,Vecteur_new1); 
    A1H_mat_stat(1,:) = pchip(Vecteur_old,A1H_mat_stat_load,Vecteur_new1);
    A2H_mat_stat(1,:) = pchip(Vecteur_old,A2H_mat_stat_load,Vecteur_new1);
    A3H_mat_stat(1,:) = pchip(Vecteur_old,A3H_mat_stat_load,Vecteur_new1);
    A1_stat(1,:) = pchip(Vecteur_old,A1_stat_load,Vecteur_new1);
    A2_stat(1,:) = pchip(Vecteur_old,A2_stat_load,Vecteur_new1);
    A3_stat(1,:) = pchip(Vecteur_old,A3_stat_load,Vecteur_new1); 
    A1tot_stat(1,:) = pchip(Vecteur_old,A1tot_stat_load,Vecteur_new1);
    A2tot_stat(1,:) = pchip(Vecteur_old,A2tot_stat_load,Vecteur_new1);
    A3tot_stat(1,:) = pchip(Vecteur_old,A3tot_stat_load,Vecteur_new1);

    HSO4_stat_new = pchip(Vecteur_old,HSO4_stat_load,Vecteur_new2);
    SO4_stat_new = pchip(Vecteur_old,SO4_stat_load,Vecteur_new2);
    SO4tot_stat_new = pchip(Vecteur_old,SO4tot_stat_load,Vecteur_new2);
    A1H_HSO4_stat_new = pchip(Vecteur_old,A1H_HSO4_stat_load,Vecteur_new2);
    A2H_HSO4_stat_new = pchip(Vecteur_old,A2H_HSO4_stat_load,Vecteur_new2);
    A3H_HSO4_stat_new = pchip(Vecteur_old,A3H_HSO4_stat_load,Vecteur_new2);
    A1H_SO4_stat_new = pchip(Vecteur_old,A1H_SO4_stat_load,Vecteur_new2);
    A2H_SO4_stat_new = pchip(Vecteur_old,A2H_SO4_stat_load,Vecteur_new2);
    A3H_SO4_stat_new = pchip(Vecteur_old,A3H_SO4_stat_load,Vecteur_new2);
    A1H_mat_stat_new = pchip(Vecteur_old,A1H_mat_stat_load,Vecteur_new2);
    A2H_mat_stat_new = pchip(Vecteur_old,A2H_mat_stat_load,Vecteur_new2);
    A3H_mat_stat_new = pchip(Vecteur_old,A3H_mat_stat_load,Vecteur_new2);
    A1_stat_new = pchip(Vecteur_old,A1_stat_load,Vecteur_new2);
    A2_stat_new = pchip(Vecteur_old,A2_stat_load,Vecteur_new2);
    A3_stat_new = pchip(Vecteur_old,A3_stat_load,Vecteur_new2);
    A1tot_stat_new = pchip(Vecteur_old,A1tot_stat_load,Vecteur_new2);
    A2tot_stat_new = pchip(Vecteur_old,A2tot_stat_load,Vecteur_new2);
    A3tot_stat_new = pchip(Vecteur_old,A3tot_stat_load,Vecteur_new2);

    HSO4_stat_old = pchip(Vecteur_old,HSO4_stat_load,Vecteur_new2);
    SO4_stat_old = pchip(Vecteur_old,SO4_stat_load,Vecteur_new2);
    SO4tot_stat_old = pchip(Vecteur_old,SO4tot_stat_load,Vecteur_new2);
    A1H_HSO4_stat_old = pchip(Vecteur_old,A1H_HSO4_stat_load,Vecteur_new2);
    A2H_HSO4_stat_old = pchip(Vecteur_old,A2H_HSO4_stat_load,Vecteur_new2);
    A3H_HSO4_stat_old = pchip(Vecteur_old,A3H_HSO4_stat_load,Vecteur_new2);
    A1H_SO4_stat_old = pchip(Vecteur_old,A1H_SO4_stat_load,Vecteur_new2);
    A2H_SO4_stat_old = pchip(Vecteur_old,A2H_SO4_stat_load,Vecteur_new2);
    A3H_SO4_stat_old = pchip(Vecteur_old,A3H_SO4_stat_load,Vecteur_new2);
    A1H_mat_stat_old = pchip(Vecteur_old,A1H_mat_stat_load,Vecteur_new2);
    A2H_mat_stat_old = pchip(Vecteur_old,A2H_mat_stat_load,Vecteur_new2);
    A3H_mat_stat_old = pchip(Vecteur_old,A3H_mat_stat_load,Vecteur_new2);
    A1_stat_old = pchip(Vecteur_old,A1_stat_load,Vecteur_new2);
    A2_stat_old = pchip(Vecteur_old,A2_stat_load,Vecteur_new2);
    A3_stat_old = pchip(Vecteur_old,A3_stat_load,Vecteur_new2);
    A1tot_stat_old = pchip(Vecteur_old,A1tot_stat_load,Vecteur_new2);
    A2tot_stat_old = pchip(Vecteur_old,A2tot_stat_load,Vecteur_new2);
    A3tot_stat_old = pchip(Vecteur_old,A3tot_stat_load,Vecteur_new2);

    F_ionique(1,:) = pchip(Vecteur_old,F_ionique_load,Vecteur_new1);
    Activite_AH(1,:) = pchip(Vecteur_old,Activite_AH_load,Vecteur_new1);
    Activite_HSO4(1,:) = pchip(Vecteur_old,Activite_HSO4_load,Vecteur_new1);
    
    F_ionique_new = pchip(Vecteur_old,F_ionique_load,Vecteur_new2);
    Activite_AH_new = pchip(Vecteur_old,Activite_AH_load,Vecteur_new2);
    Activite_HSO4_new = pchip(Vecteur_old,Activite_HSO4_load,Vecteur_new2);
    
else
    Vecteur_old = linspace (1,NETtot,N_etages);
    Vecteur_new2 = 1:NETtot;
    
    H_mob(1,:) = H_mob_ini*ones(1,N_etages);
    HSO4_mob(1,:) = HSO4_mob_ini*ones(1,N_etages);
    SO4_mob(1,:) = SO4_mob_ini*ones(1,N_etages);
    SO4tot_mob(1,:) = SO4tot_mob_ini*ones(1,N_etages);
    A1H_mob(1,:) = A1H_mob_ini*ones(1,N_etages);
    A2H_mob(1,:) = A2H_mob_ini*ones(1,N_etages);
    A3H_mob(1,:) = A3H_mob_ini*ones(1,N_etages);
    A1_mob(1,:) = A1_mob_ini*ones(1,N_etages);
    A2_mob(1,:) = A2_mob_ini*ones(1,N_etages);
    A3_mob(1,:) = A3_mob_ini*ones(1,N_etages);    
    A1tot_mob(1,:) = A1tot_mob_ini*ones(1,N_etages);
    A2tot_mob(1,:) = A2tot_mob_ini*ones(1,N_etages);
    A3tot_mob(1,:) = A3tot_mob_ini*ones(1,N_etages);
    
    H_mob_new = H_mob_ini*ones(1,NETtot);
    HSO4_mob_new = HSO4_mob_ini*ones(1,NETtot);
    SO4_mob_new = SO4_mob_ini*ones(1,NETtot);
    SO4tot_mob_new = SO4tot_mob_ini*ones(1,NETtot);
    A1H_mob_new = A1H_mob_ini*ones(1,NETtot);
    A2H_mob_new = A2H_mob_ini*ones(1,NETtot);
    A3H_mob_new = A3H_mob_ini*ones(1,NETtot);
    A1_mob_new = A1_mob_ini*ones(1,NETtot);
    A2_mob_new = A2_mob_ini*ones(1,NETtot);
    A3_mob_new = A3_mob_ini*ones(1,NETtot);
    A1tot_mob_new = A1tot_mob_ini*ones(1,NETtot);
    A2tot_mob_new = A2tot_mob_ini*ones(1,NETtot);
    A3tot_mob_new = A3tot_mob_ini*ones(1,NETtot);
    
    H_mob_old = H_mob_ini*ones(1,NETtot);
    HSO4_mob_old = HSO4_mob_ini*ones(1,NETtot);
    SO4_mob_old = SO4_mob_ini*ones(1,NETtot);
    SO4tot_mob_old = SO4tot_mob_ini*ones(1,NETtot);
    A1H_mob_old = A1H_mob_ini*ones(1,NETtot);
    A2H_mob_old = A2H_mob_ini*ones(1,NETtot);
    A3H_mob_old = A3H_mob_ini*ones(1,NETtot);
    A1_mob_old = A1_mob_ini*ones(1,NETtot);
    A2_mob_old = A2_mob_ini*ones(1,NETtot);
    A3_mob_old = A3_mob_ini*ones(1,NETtot);
    A1tot_mob_old = A1tot_mob_ini*ones(1,NETtot);
    A2tot_mob_old = A2tot_mob_ini*ones(1,NETtot);
    A3tot_mob_old = A3tot_mob_ini*ones(1,NETtot);
    
    HSO4_stat(1,:) = HSO4_stat_ini*ones(1,N_etages);
    SO4_stat(1,:) = SO4_stat_ini*ones(1,N_etages);
    SO4tot_stat(1,:) = SO4tot_stat_ini*ones(1,N_etages);
    A1H_HSO4_stat(1,:) = A1H_HSO4_stat_ini*ones(1,N_etages);
    A2H_HSO4_stat(1,:) = A2H_HSO4_stat_ini*ones(1,N_etages);
    A3H_HSO4_stat(1,:) = A3H_HSO4_stat_ini*ones(1,N_etages);
    A1H_SO4_stat(1,:) = A1H_SO4_stat_ini*ones(1,N_etages);
    A2H_SO4_stat(1,:) = A2H_SO4_stat_ini*ones(1,N_etages);
    A3H_SO4_stat(1,:) = A3H_SO4_stat_ini*ones(1,N_etages); 
    A1H_mat_stat(1,:) = A1H_mat_stat_ini*ones(1,N_etages);
    A2H_mat_stat(1,:) = A2H_mat_stat_ini*ones(1,N_etages);
    A3H_mat_stat(1,:) = A3H_mat_stat_ini*ones(1,N_etages);
    A1_stat(1,:) = A1_stat_ini*ones(1,N_etages);
    A2_stat(1,:) = A2_stat_ini*ones(1,N_etages);
    A3_stat(1,:) = A3_stat_ini*ones(1,N_etages); 
    A1tot_stat(1,:) = A1tot_stat_ini*ones(1,N_etages);
    A2tot_stat(1,:) = A2tot_stat_ini*ones(1,N_etages);
    A3tot_stat(1,:) = A3tot_stat_ini*ones(1,N_etages);

    HSO4_stat_new = HSO4_stat_ini*ones(1,NETtot);
    SO4_stat_new = SO4_stat_ini*ones(1,NETtot);
    SO4tot_stat_new = SO4tot_stat_ini*ones(1,NETtot);
    A1H_HSO4_stat_new = A1H_HSO4_stat_ini*ones(1,NETtot);
    A2H_HSO4_stat_new = A2H_HSO4_stat_ini*ones(1,NETtot);
    A3H_HSO4_stat_new = A3H_HSO4_stat_ini*ones(1,NETtot);
    A1H_SO4_stat_new = A1H_SO4_stat_ini*ones(1,NETtot);
    A2H_SO4_stat_new = A2H_SO4_stat_ini*ones(1,NETtot);
    A3H_SO4_stat_new = A3H_SO4_stat_ini*ones(1,NETtot);
    A1H_mat_stat_new = A1H_mat_stat_ini*ones(1,NETtot);
    A2H_mat_stat_new = A2H_mat_stat_ini*ones(1,NETtot);
    A3H_mat_stat_new = A3H_mat_stat_ini*ones(1,NETtot);
    A1_stat_new = A1_stat_ini*ones(1,NETtot);
    A2_stat_new = A2_stat_ini*ones(1,NETtot);
    A3_stat_new = A3_stat_ini*ones(1,NETtot);
    A1tot_stat_new = A1tot_stat_ini*ones(1,NETtot);
    A2tot_stat_new = A2tot_stat_ini*ones(1,NETtot);
    A3tot_stat_new = A3tot_stat_ini*ones(1,NETtot);

    HSO4_stat_old = HSO4_stat_ini*ones(1,NETtot);
    SO4_stat_old = SO4_stat_ini*ones(1,NETtot);
    SO4tot_stat_old = SO4tot_stat_ini*ones(1,NETtot);
    A1H_HSO4_stat_old = A1H_HSO4_stat_ini*ones(1,NETtot);
    A2H_HSO4_stat_old = A2H_HSO4_stat_ini*ones(1,NETtot);
    A3H_HSO4_stat_old = A3H_HSO4_stat_ini*ones(1,NETtot);
    A1H_SO4_stat_old = A1H_SO4_stat_ini*ones(1,NETtot);
    A2H_SO4_stat_old = A2H_SO4_stat_ini*ones(1,NETtot);
    A3H_SO4_stat_old = A3H_SO4_stat_ini*ones(1,NETtot);
    A1H_mat_stat_old = A1H_mat_stat_ini*ones(1,NETtot);
    A2H_mat_stat_old = A2H_mat_stat_ini*ones(1,NETtot);
    A3H_mat_stat_old = A3H_mat_stat_ini*ones(1,NETtot);
    A1_stat_old = A1_stat_ini*ones(1,NETtot);
    A2_stat_old = A2_stat_ini*ones(1,NETtot);
    A3_stat_old = A3_stat_ini*ones(1,NETtot);
    A1tot_stat_old = A1tot_stat_ini*ones(1,NETtot);
    A2tot_stat_old = A2tot_stat_ini*ones(1,NETtot);
    A3tot_stat_old = A3tot_stat_ini*ones(1,NETtot);

    F_ionique(1,:) = F_ionique_ini*ones(1,N_etages);
    F_ionique_new = F_ionique_ini*ones(1,NETtot);
    Activite_AH(1,:) = Activite_AH_ini*ones(1,N_etages);
    Activite_AH_new = Activite_AH_ini*ones(1,NETtot);
    Activite_HSO4(1,:) = Activite_HSO4_ini*ones(1,N_etages);
    Activite_HSO4_new = Activite_HSO4_ini*ones(1,NETtot);
end

H_ext_total = zeros(N_volume,1);
H_ext_cycle = zeros(N_volume,1);
H_ext_etape = zeros(N_volume,1);
H_ext_sous_etape = zeros(N_volume,1);
HSO4_ext_total = zeros(N_volume,1);
HSO4_ext_cycle = zeros(N_volume,1);
HSO4_ext_etape = zeros(N_volume,1);
HSO4_ext_sous_etape = zeros(N_volume,1);
SO4_ext_total = zeros(N_volume,1);
SO4_ext_cycle = zeros(N_volume,1);
SO4_ext_etape = zeros(N_volume,1);
SO4_ext_sous_etape = zeros(N_volume,1);
A1H_ext_total = zeros(N_volume,1);
A1H_ext_cycle = zeros(N_volume,1);
A1H_ext_etape = zeros(N_volume,1);
A1H_ext_sous_etape = zeros(N_volume,1);
A1_ext_total = zeros(N_volume,1);
A1_ext_cycle = zeros(N_volume,1);
A1_ext_etape = zeros(N_volume,1);
A1_ext_sous_etape = zeros(N_volume,1);
A2H_ext_total = zeros(N_volume,1);
A2H_ext_cycle = zeros(N_volume,1);
A2H_ext_etape = zeros(N_volume,1);
A2H_ext_sous_etape = zeros(N_volume,1);
A2_ext_total = zeros(N_volume,1);
A2_ext_cycle = zeros(N_volume,1);
A2_ext_etape = zeros(N_volume,1);
A2_ext_sous_etape = zeros(N_volume,1);
A3H_ext_total = zeros(N_volume,1);
A3H_ext_cycle = zeros(N_volume,1);
A3H_ext_etape = zeros(N_volume,1);
A3H_ext_sous_etape = zeros(N_volume,1);
A3_ext_total = zeros(N_volume,1);
A3_ext_cycle = zeros(N_volume,1);
A3_ext_etape = zeros(N_volume,1);
A3_ext_sous_etape = zeros(N_volume,1);

H_raff_total = zeros(N_volume,1);
H_raff_cycle = zeros(N_volume,1);
H_raff_etape = zeros(N_volume,1);
H_raff_sous_etape = zeros(N_volume,1);
HSO4_raff_total = zeros(N_volume,1);
HSO4_raff_cycle = zeros(N_volume,1);
HSO4_raff_etape = zeros(N_volume,1);
HSO4_raff_sous_etape = zeros(N_volume,1);
SO4_raff_total = zeros(N_volume,1);
SO4_raff_cycle = zeros(N_volume,1);
SO4_raff_etape = zeros(N_volume,1);
SO4_raff_sous_etape = zeros(N_volume,1);
A1H_raff_total = zeros(N_volume,1);
A1H_raff_cycle = zeros(N_volume,1);
A1H_raff_etape = zeros(N_volume,1);
A1H_raff_sous_etape = zeros(N_volume,1);
A1_raff_total = zeros(N_volume,1);
A1_raff_cycle = zeros(N_volume,1);
A1_raff_etape = zeros(N_volume,1);
A1_raff_sous_etape = zeros(N_volume,1);
A2H_raff_total = zeros(N_volume,1);
A2H_raff_cycle = zeros(N_volume,1);
A2H_raff_etape = zeros(N_volume,1);
A2H_raff_sous_etape = zeros(N_volume,1);
A2_raff_total = zeros(N_volume,1);
A2_raff_cycle = zeros(N_volume,1);
A2_raff_etape = zeros(N_volume,1);
A2_raff_sous_etape = zeros(N_volume,1);
A3H_raff_total = zeros(N_volume,1);
A3H_raff_cycle = zeros(N_volume,1);
A3H_raff_etape = zeros(N_volume,1);
A3H_raff_sous_etape = zeros(N_volume,1);
A3_raff_total = zeros(N_volume,1);
A3_raff_cycle = zeros(N_volume,1);
A3_raff_etape = zeros(N_volume,1);
A3_raff_sous_etape = zeros(N_volume,1);

N_iter_loc_1 = zeros(N_volume,N_etages);
N_iter_loc_1(1,:) = ones(1,N_etages);
N_iter_loc_1_tot = zeros(N_volume,1);
N_iter_loc_1_tot(1) = N_etages;
N_iter_loc_2 = zeros(N_volume,N_etages);
N_iter_loc_2(1,:) = ones(1,N_etages);
N_iter_loc_2_tot = zeros(N_volume,1);
N_iter_loc_2_tot(1) = N_etages;
N_iter_glob = ones(N_volume,1);
Ecart_abs_max = ones(N_volume,N_etages);
Ecart_abs_max(1,:) = zeros(1,N_etages);
Ecart_rel_max = ones(N_volume,N_etages);
Ecart_rel_max(1,:) = zeros(1,N_etages);
Residu_abs_max = ones(N_volume,N_etages);
Residu_abs_max(1,:) = zeros(1,N_etages);
Residu_rel_max = ones(N_volume,N_etages);
Residu_rel_max(1,:) = zeros(1,N_etages);
Convergence = zeros(N_volume,N_etages);
Convergence(1,:) = ones(1,N_etages);
Modele = zeros(N_volume,N_etages);

% Création du graphique
figure(1)

V_etape_fig = num2str(V_inj_etape(1),3);
V_sous_etape_fig = num2str(V_inj_sous_etape(1),3);
Cycle_fig = num2str(Cycle(1),3);
Etape_fig = num2str(Etape(1),3);
Sous_Etape_fig = num2str(Sous_Etape(1),3);

    % H+, HSO4- et SO4 2- en phase mobile
    subplot(4,1,1)
    plot(Vecteur_new1,H_mob(1,:),'LineStyle','none','Marker','o','MarkerSize',2,'Color',[0.8 0.8 0.8]);
    hold on
    plot(Vecteur_new1,HSO4_mob(1,:),'LineStyle','none','Marker','o','MarkerSize',2,'Color',[0.8 0.8 0]);
    plot(Vecteur_new1,SO4_mob(1,:),'LineStyle','none','Marker','o','MarkerSize',2,'Color',[0.8 0 0.8]);
    
    Titre = ['Cycle = ',Cycle_fig,', Etape = ',Etape_fig,', Sous-Etape = ',Sous_Etape_fig,', Vétape = ',V_etape_fig,' BV, Vsous-étape = ',V_sous_etape_fig,' BV'];
    title(Titre);
    xlabel('Etage théorique');
    ylabel('Concentration (mol/L)');
    xlim([0 NETtot*1.09]);
    Ymax_mob_1 = max([max(H_mob(1,:));max(HSO4_mob(1,:));max(SO4_mob(1,:))]);
    ylim([0 Ymax_mob_1]);
    
    plot(Entree_Z1(1)*ones(3,1),linspace(0,Ymax_mob_1,3),'r--','LineWidth',2);
    plot(Entree_Z2(1)*ones(3,1),linspace(0,Ymax_mob_1,3),'g--','LineWidth',2);
    plot(Entree_Z3(1)*ones(3,1),linspace(0,Ymax_mob_1,3),'c--','LineWidth',2);
    plot(Entree_Z4(1)*ones(3,1),linspace(0,Ymax_mob_1,3),'b--','LineWidth',2);
    plot(Sortie_Z1(1)*ones(3,1),linspace(0,Ymax_mob_1,3),'r-.','MarkerSize',2);
    plot(Sortie_Z2(1)*ones(3,1),linspace(0,Ymax_mob_1,3),'g-.','MarkerSize',2);
    plot(Sortie_Z3(1)*ones(3,1),linspace(0,Ymax_mob_1,3),'c-.','MarkerSize',2);
    plot(Sortie_Z4(1)*ones(3,1),linspace(0,Ymax_mob_1,3),'b-.','MarkerSize',2);
    
    hold off
    hleg1 = legend('H+','HSO4-','SO4 2-');
    set(hleg1,'FontSize',9);
    text(Milieu_Z1(1),Ymax_mob_1/2,'Z1','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','red')
    text(Milieu_Z2(1),Ymax_mob_1/2,'Z2','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','green')
    text(Milieu_Z3(1),Ymax_mob_1/2,'Z3','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','cyan')
    text(Milieu_Z4(1),Ymax_mob_1/2,'Z4','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','blue')
    
    % A1tot, A2tot, A3tot en phase mobile
    subplot(4,1,2)
    plot(Vecteur_new1,A1tot_mob(1,:),'LineStyle','none','Marker','o','MarkerSize',2,'Color',[1 0 0]);
    hold on
    plot(Vecteur_new1,A2tot_mob(1,:),'LineStyle','none','Marker','o','MarkerSize',2,'Color',[0 1 0]);
    plot(Vecteur_new1,A3tot_mob(1,:),'LineStyle','none','Marker','o','MarkerSize',2,'Color',[0 0 1]);
    
    xlabel('Etage théorique');
    ylabel('Concentration (mol/L)');
    xlim([0 NETtot*1.09]);
    Ymax_mob_2 = max([max(A1tot_mob(1,:));max(A2tot_mob(1,:));max(A3tot_mob(1,:))]);
    if Ymax_mob_2>0
        ylim([0 Ymax_mob_2]);
    else
        Ymax_mob_2 = Ymax_mob_1;
        ylim([0 Ymax_mob_1]);
    end
    
    plot(Entree_Z1(1)*ones(3,1),linspace(0,Ymax_mob_2,3),'r--','LineWidth',2);
    plot(Entree_Z2(1)*ones(3,1),linspace(0,Ymax_mob_2,3),'g--','LineWidth',2);
    plot(Entree_Z3(1)*ones(3,1),linspace(0,Ymax_mob_2,3),'c--','LineWidth',2);
    plot(Entree_Z4(1)*ones(3,1),linspace(0,Ymax_mob_2,3),'b--','LineWidth',2);
    plot(Sortie_Z1(1)*ones(3,1),linspace(0,Ymax_mob_2,3),'r-.','MarkerSize',2);
    plot(Sortie_Z2(1)*ones(3,1),linspace(0,Ymax_mob_2,3),'g-.','MarkerSize',2);
    plot(Sortie_Z3(1)*ones(3,1),linspace(0,Ymax_mob_2,3),'c-.','MarkerSize',2);
    plot(Sortie_Z4(1)*ones(3,1),linspace(0,Ymax_mob_2,3),'b-.','MarkerSize',2);
    
    hold off
    hleg2 = legend('A1 tot','A2 tot','A3 tot');
    set(hleg2,'FontSize',9);
    text(Milieu_Z1(1),Ymax_mob_2/2,'Z1','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','red')
    text(Milieu_Z2(1),Ymax_mob_2/2,'Z2','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','green')
    text(Milieu_Z3(1),Ymax_mob_2/2,'Z3','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','cyan')
    text(Milieu_Z4(1),Ymax_mob_2/2,'Z4','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','blue')
    
    % HSO4- et SO4 2- en phase stationnaire
    subplot(4,1,3)
    plot(Vecteur_new1,HSO4_stat(1,:),'LineStyle','none','Marker','o','MarkerSize',2,'Color',[0.8 0.8 0]);
    hold on    
    plot(Vecteur_new1,SO4_stat(1,:),'LineStyle','none','Marker','o','MarkerSize',2,'Color',[0.8 0 0.8]);     

    xlabel('Etage théorique');
    ylabel('Concentration (mol/Lrésine)');
    xlim([0 NETtot*1.09]);
    Ymax_stat_1 = max([max(HSO4_stat(1,:));max(SO4_stat(1,:))]);
    ylim([0 Ymax_stat_1]);
    
    plot(Entree_Z1(1)*ones(3,1),linspace(0,Ymax_stat_1,3),'r--','LineWidth',2);
    plot(Entree_Z2(1)*ones(3,1),linspace(0,Ymax_stat_1,3),'g--','LineWidth',2);
    plot(Entree_Z3(1)*ones(3,1),linspace(0,Ymax_stat_1,3),'c--','LineWidth',2);
    plot(Entree_Z4(1)*ones(3,1),linspace(0,Ymax_stat_1,3),'b--','LineWidth',2);
    plot(Sortie_Z1(1)*ones(3,1),linspace(0,Ymax_stat_1,3),'r-.','MarkerSize',2);
    plot(Sortie_Z2(1)*ones(3,1),linspace(0,Ymax_stat_1,3),'g-.','MarkerSize',2);
    plot(Sortie_Z3(1)*ones(3,1),linspace(0,Ymax_stat_1,3),'c-.','MarkerSize',2);
    plot(Sortie_Z4(1)*ones(3,1),linspace(0,Ymax_stat_1,3),'b-.','MarkerSize',2);
    
    hold off
    hleg3 = legend('HSO4-','SO4 2-');
    set(hleg3,'FontSize',9);
    text(Milieu_Z1(1),Ymax_stat_1/2,'Z1','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','red')
    text(Milieu_Z2(1),Ymax_stat_1/2,'Z2','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','green')
    text(Milieu_Z3(1),Ymax_stat_1/2,'Z3','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','cyan')
    text(Milieu_Z4(1),Ymax_stat_1/2,'Z4','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','blue')
    
    % A1adsorbé, A1échangé, A2adsorbé, A2échangé, A3adsorbé et A3échangé en phase stationnaire
    subplot(4,1,4)
    A1H_stat=A1tot_stat(1,:)-A1_stat(1,:);
    A2H_stat=A2tot_stat(1,:)-A2_stat(1,:);
    A3H_stat=A3tot_stat(1,:)-A3_stat(1,:);
    plot(Vecteur_new1,A1H_stat,'LineStyle','none','Marker','o','MarkerSize',2,'Color',[1 0 0]);
    hold on
    plot(Vecteur_new1,A1_stat(1,:),'LineStyle','none','Marker','o','MarkerSize',2,'Color',[0.6 0 0]);
    plot(Vecteur_new1,A2H_stat,'LineStyle','none','Marker','o','MarkerSize',2,'Color',[0 1 0]);
    plot(Vecteur_new1,A2_stat(1,:),'LineStyle','none','Marker','o','MarkerSize',2,'Color',[0 0.6 0]);
    plot(Vecteur_new1,A3H_stat,'LineStyle','none','Marker','o','MarkerSize',2,'Color',[0 0 1]);
    plot(Vecteur_new1,A3_stat(1,:),'LineStyle','none','Marker','o','MarkerSize',2,'Color',[0 0 0.6]);

    xlabel('Etage théorique');
    ylabel('Concentration (mol/Lrésine)');
    xlim([0 NETtot*1.09]);
    Ymax_stat_2 = max([max(A1H_stat(1,:));max(A1_stat(1,:));max(A2H_stat(1,:));max(A2_stat(1,:));max(A3H_stat(1,:));max(A3_stat(1,:))]);
    if Ymax_stat_2>0
        ylim([0 Ymax_stat_2]);
    else
        Ymax_stat_2 = Ymax_stat_1;
        ylim([0 Ymax_stat_1]);
    end
    
    plot(Entree_Z1(1)*ones(3,1),linspace(0,Ymax_stat_2,3),'r--','LineWidth',2);
    plot(Entree_Z2(1)*ones(3,1),linspace(0,Ymax_stat_2,3),'g--','LineWidth',2);
    plot(Entree_Z3(1)*ones(3,1),linspace(0,Ymax_stat_2,3),'c--','LineWidth',2);
    plot(Entree_Z4(1)*ones(3,1),linspace(0,Ymax_stat_2,3),'b--','LineWidth',2);
    plot(Sortie_Z1(1)*ones(3,1),linspace(0,Ymax_stat_2,3),'r-.','MarkerSize',2);
    plot(Sortie_Z2(1)*ones(3,1),linspace(0,Ymax_stat_2,3),'g-.','MarkerSize',2);
    plot(Sortie_Z3(1)*ones(3,1),linspace(0,Ymax_stat_2,3),'c-.','MarkerSize',2);
    plot(Sortie_Z4(1)*ones(3,1),linspace(0,Ymax_stat_2,3),'b-.','MarkerSize',2);
    
    hold off
    hleg4 = legend('A1H','A1-','A2H','A2-','A3H','A3-');
    set(hleg4,'FontSize',9);
    text(Milieu_Z1(1),Ymax_stat_2/2,'Z1','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','red')
    text(Milieu_Z2(1),Ymax_stat_2/2,'Z2','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','green')
    text(Milieu_Z3(1),Ymax_stat_2/2,'Z3','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','cyan')
    text(Milieu_Z4(1),Ymax_stat_2/2,'Z4','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','blue')
    
    pause(10)
    
% Simulation en fonction du volume injecté
% Initialisation SSMB
i = 1;
V_inj_total_new = 0;
V_ech_total_new = 0;
V_elu_total_new = 0;
V_raff_total_new = 0;
V_ext_total_new = 0;
H_raff_total_new = 0;
HSO4_raff_total_new = 0;
SO4_raff_total_new = 0;
A1H_raff_total_new = 0;
A1_raff_total_new = 0;
A2H_raff_total_new = 0;
A2_raff_total_new = 0;
A3H_raff_total_new = 0;
A3_raff_total_new = 0;
H_ext_total_new = 0;
HSO4_ext_total_new = 0;
SO4_ext_total_new = 0;
A1H_ext_total_new = 0;
A1_ext_total_new = 0;
A2H_ext_total_new = 0;
A2_ext_total_new = 0;
A3H_ext_total_new = 0;
A3_ext_total_new = 0;
Cycle_new = 1;

while Cycle_new <= Ncycles
    % Initialisation cycle
    V_inj_cycle_new = 0;
    V_ech_cycle_new = 0;
    V_elu_cycle_new = 0;
    V_raff_cycle_new = 0;
    V_ext_cycle_new = 0;
    H_raff_cycle_new = 0;
    HSO4_raff_cycle_new = 0;
    SO4_raff_cycle_new = 0;
    A1H_raff_cycle_new = 0;
    A1_raff_cycle_new = 0;
    A2H_raff_cycle_new = 0;
    A2_raff_cycle_new = 0;
    A3H_raff_cycle_new = 0;
    A3_raff_cycle_new = 0;
    H_ext_cycle_new = 0;
    HSO4_ext_cycle_new = 0;
    SO4_ext_cycle_new = 0;
    A1H_ext_cycle_new = 0;
    A1_ext_cycle_new = 0;
    A2H_ext_cycle_new = 0;
    A2_ext_cycle_new = 0;
    A3H_ext_cycle_new = 0;
    A3_ext_cycle_new = 0;  
    Etape_new = 1;
    
    while Etape_new <= Ncol_tot
        % Initialisation étape
        V_inj_etape_new = 0;
        V_ech_etape_new = 0;
        V_elu_etape_new = 0;
        V_raff_etape_new = 0;
        V_ext_etape_new = 0;
        H_raff_etape_new = 0;
        HSO4_raff_etape_new = 0;
        SO4_raff_etape_new = 0;
        A1H_raff_etape_new = 0;
        A1_raff_etape_new = 0;
        A2H_raff_etape_new = 0;
        A2_raff_etape_new = 0;
        A3H_raff_etape_new = 0;
        A3_raff_etape_new = 0;
        H_ext_etape_new = 0;
        HSO4_ext_etape_new = 0;
        SO4_ext_etape_new = 0;
        A1H_ext_etape_new = 0;
        A1_ext_etape_new = 0;
        A2H_ext_etape_new = 0;
        A2_ext_etape_new = 0;
        A3H_ext_etape_new = 0;
        A3_ext_etape_new = 0;
        Sous_Etape_new = 1;
        
        while Sous_Etape_new <= 4
            % Initialisation sous-étape
            V_inj_sous_etape_old = 0;
            V_inj_sous_etape_new = 0;
            V_ech_sous_etape_new = 0;
            V_elu_sous_etape_new = 0;
            V_raff_sous_etape_new = 0;
            V_ext_sous_etape_new = 0;
            H_raff_sous_etape_new = 0;
            HSO4_raff_sous_etape_new = 0;
            SO4_raff_sous_etape_new = 0;
            A1H_raff_sous_etape_new = 0;
            A1_raff_sous_etape_new = 0;
            A2H_raff_sous_etape_new = 0;
            A2_raff_sous_etape_new = 0;
            A3H_raff_sous_etape_new = 0;
            A3_raff_sous_etape_new = 0;
            H_ext_sous_etape_new = 0;
            HSO4_ext_sous_etape_new = 0;
            SO4_ext_sous_etape_new = 0;
            A1H_ext_sous_etape_new = 0;
            A1_ext_sous_etape_new = 0;
            A2H_ext_sous_etape_new = 0;
            A2_ext_sous_etape_new = 0;
            A3H_ext_sous_etape_new = 0;
            A3_ext_sous_etape_new = 0;
            
            % Définition de la zone de travail et de l'alimentation
            if Sous_Etape_new == 1
                n_debut = Entree_Z3(Etape_new,1);
                n_fin = Sortie_Z3(Etape_new,1);
                NET_SE = NET_SE1;
                V_alim = F_V_ech(Cycle_new)*V_ech_est;
                H_alim = H_ech;
                HSO4_alim = HSO4_ech;
                SO4_alim = SO4_ech;                       
                A1H_alim = A1H_ech;
                A1_alim = A1_ech;
                A2H_alim = A2H_ech;
                A2_alim = A2_ech;
                A3H_alim = A3H_ech;
                A3_alim = A3_ech;                        
                
            elseif Sous_Etape_new == 2
                n_debut = Entree_Z1(Etape_new,1);
                n_fin = Sortie_Z1(Etape_new,1);
                NET_SE = NET_SE2;
                V_alim = F_V_elu1(Cycle_new)*V_elu1_est;
                H_alim = H_elu;
                HSO4_alim = HSO4_elu;
                SO4_alim = SO4_elu;                       
                A1H_alim = A1H_elu;
                A1_alim = A1_elu;
                A2H_alim = A2H_elu;
                A2_alim = A2_elu;
                A3H_alim = A3H_elu;
                A3_alim = A3_elu;
                
            elseif Sous_Etape_new == 3
                n_debut = Entree_Z1(Etape_new,1);
                n_fin = Sortie_Z3(Etape_new,1);
                NET_SE = NET_SE3;
                V_alim = F_V_elu2(Cycle_new)*V_elu2_est;            
                
            else
                n_debut = Entree_Z1(Etape_new,1);
                n_fin = Sortie_Z4(Etape_new,1);
                NET_SE = NET_SE4;
                V_alim = F_V_boucle(Cycle_new)*V_boucle_est; 
            end
            
            while V_inj_sous_etape_new < V_alim
                % Pas par défaut sur le volume injecté
                Pas_V = Pas_V_calcul;
                V_inj_sous_etape_new = V_inj_sous_etape_new + Pas_V; 

                % Calculs dans les étages théoriques des colonnes de la zone de travail
                n = n_debut;
                N_iter_loc_1_tot_new = 0;
                N_iter_loc_2_tot_new = 0;                    
                for k=1:NET_SE
                    N_iter_glob_new = 1;
                    N_iter_loc_1_new = 0;
                    N_iter_loc_2_new = 0;
                    Convergence_new = 0;                  
                    Ecart_abs_max_new = 2 * Precision_Conc_abs;
                    Ecart_rel_max_new = 2 * Precision_Conc_rel;
                    Residu_abs_max_new = 2 * Precision_Bilan_abs;
                    Residu_rel_max_new = 2 * Precision_Bilan_rel;
                    while ((Ecart_abs_max_new >= Precision_Conc_abs && Ecart_rel_max_new >= Precision_Conc_rel) || (Residu_abs_max_new >= Precision_Bilan_abs && Residu_rel_max_new >= Precision_Bilan_rel)) && (N_iter_glob_new < N_iter_max_pas_V || N_iter_loc_1_new < N_iter_max_Newton)
                        N_iter_loc_1_new = N_iter_loc_1_new + 1;
                        N_iter_loc_1_tot_new = N_iter_loc_1_tot_new + 1;

                        if N_iter_loc_1_new > N_iter_max_Newton
                            % Diminution du pas sur le volume injecté pour converger
                            V_inj_sous_etape_new = V_inj_sous_etape_new - Pas_V;
                            Pas_V = Pas_V / Facteur_diminution_pas_V;
                            V_inj_sous_etape_new = V_inj_sous_etape_new + Pas_V;

                            % Recalcul à partir du 1er étage théorique (n_debut)
                            N_iter_glob_new = N_iter_glob_new + 1;
                            N_iter_loc_1_new = 1;
                            N_iter_loc_2_new = 0;
                            Convergence_new = 0;
                            n = n_debut;
                        end

                        if N_iter_loc_1_new == 1
                            if n == n_debut
                                % Ajustement de V_inj_sous_etape_new et Pas_V par rapport à V_alim
                                if V_inj_sous_etape_new > V_alim
                                    V_inj_sous_etape_new = V_inj_sous_etape_new - Pas_V;
                                    Pas_V = V_alim - V_inj_sous_etape_new;
                                    V_inj_sous_etape_new = V_alim;
                                end

                                % Calculs intermédiaires
                                C4 = Pas_V*NETcol/Porosite;
                                C5 = 1 + C4;
                                C13 = C5*Ka_HSO4_0;
                                C14_1 = C5*Ka_A1H_0;
                                C14_2 = C5*Ka_A2H_0;
                                C14_3 = C5*Ka_A3H_0;
                    
                                % Définition composition en entrée                                
                                if Sous_Etape_new <= 3
                                    H_entree = H_alim;
                                    HSO4_entree = HSO4_alim;
                                    SO4_entree = SO4_alim;  
                                    A1H_entree = A1H_alim;
                                    A2H_entree = A2H_alim;
                                    A3H_entree = A3H_alim;
                                    A1_entree = A1_alim;
                                    A2_entree = A2_alim;
                                    A3_entree = A3_alim;
                                    
                                else
                                    H_entree = H_mob_old(n_fin);
                                    HSO4_entree = HSO4_mob_old(n_fin);
                                    SO4_entree = SO4_mob_old(n_fin);                       
                                    A1H_entree = A1H_mob_old(n_fin);
                                    A1_entree = A1_mob_old(n_fin);
                                    A2H_entree = A2H_mob_old(n_fin);
                                    A2_entree = A2_mob_old(n_fin);
                                    A3H_entree = A3H_mob_old(n_fin);
                                    A3_entree = A3_mob_old(n_fin);                                     
                                end
                                
                            else
                                if n == 1
                                    H_entree = H_mob_new(NETtot);
                                    HSO4_entree = HSO4_mob_new(NETtot);
                                    SO4_entree = SO4_mob_new(NETtot);
                                    A1H_entree = A1H_mob_new(NETtot);
                                    A2H_entree = A2H_mob_new(NETtot);
                                    A3H_entree = A3H_mob_new(NETtot);
                                    A1_entree = A1_mob_new(NETtot);
                                    A2_entree = A2_mob_new(NETtot);
                                    A3_entree = A3_mob_new(NETtot);
                                else
                                    H_entree = H_mob_new(n-1);
                                    HSO4_entree = HSO4_mob_new(n-1);
                                    SO4_entree = SO4_mob_new(n-1);
                                    A1H_entree = A1H_mob_new(n-1);
                                    A2H_entree = A2H_mob_new(n-1);
                                    A3H_entree = A3H_mob_new(n-1);
                                    A1_entree = A1_mob_new(n-1);
                                    A2_entree = A2_mob_new(n-1);
                                    A3_entree = A3_mob_new(n-1);                                    
                                end
                            end
                            
                            % Calculs intermédiaires
                            Htot_entree = H_entree + HSO4_entree + A1H_entree + A2H_entree + A3H_entree;
                            SO4tot_entree = HSO4_entree + SO4_entree; 
                            A1tot_entree = A1H_entree + A1_entree;
                            A2tot_entree = A2H_entree + A2_entree;
                            A3tot_entree = A3H_entree + A3_entree;

                            C_H = C4 * Htot_entree + (H_mob_old(n) + HSO4_mob_old(n) + A1H_mob_old(n) + A2H_mob_old(n) + A3H_mob_old(n)) + C1 * (HSO4_stat_old(n) + A1H_HSO4_stat_old(n) + A1H_SO4_stat_old(n) + A1H_mat_stat_old(n) + A2H_HSO4_stat_old(n) + A2H_SO4_stat_old(n) + A2H_mat_stat_old(n) + A3H_HSO4_stat_old(n) + A3H_SO4_stat_old(n) + A3H_mat_stat_old(n));
                            C_SO4 = C4 * SO4tot_entree + SO4tot_mob_old(n) + C1 * SO4tot_stat_old(n);   
                            C_A1H = C4 * A1tot_entree + A1tot_mob_old(n) + C1 * A1tot_stat_old(n);
                            C_A2H = C4 * A2tot_entree + A2tot_mob_old(n) + C1 * A2tot_stat_old(n);
                            C_A3H = C4 * A3tot_entree + A3tot_mob_old(n) + C1 * A3tot_stat_old(n);

                            F_ionique_new(n) = 0.5*(A1_mob_old(n) + A2_mob_old(n) + A3_mob_old(n) + HSO4_mob_old(n) + 4*SO4_mob_old(n));               
                            Activite_AH_new(n) = 10^(0.51*F_ionique_new(n)^0.5/(1+F_ionique_new(n)^0.5));
                            Activite_HSO4_new(n) = 10^(0.51*3*F_ionique_new(n)^0.5/(1+F_ionique_new(n)^0.5));               

                            % Simplification du modèle si concentration négligeable en acides organiques
                            Modele_new = 0;                 
                            if C_A1H < Precision_Conc_abs
                                if C_A2H < Precision_Conc_abs
                                    if C_A3H < Precision_Conc_abs       
                                        Modele_new = 7;    
                                    else
                                        Modele_new = 4;                
                                    end
                                else
                                    if C_A3H < Precision_Conc_abs       
                                        Modele_new = 5;                
                                    else
                                        Modele_new = 1;                
                                    end                             
                                end
                            elseif C_A2H < Precision_Conc_abs
                                if C_A3H < Precision_Conc_abs
                                    Modele_new = 6; 
                                else
                                    Modele_new = 2;                
                                end
                            elseif C_A3H < Precision_Conc_abs
                                Modele_new = 3;                
                            end

                            % Calculs intermédiaires
                            C9_corr = C9*Activite_HSO4_new(n);
                            C12_corr = C12*Activite_HSO4_new(n);
                            C13_corr = C13*Activite_HSO4_new(n);             
                            C15_corr = C15*Activite_HSO4_new(n);                
                            C16_1_corr = C16_1*Activite_HSO4_new(n);
                            C16_2_corr = C16_2*Activite_HSO4_new(n);
                            C16_3_corr = C16_3*Activite_HSO4_new(n);
                            C20 = C5/C_SO4;
                            C21 = C13_corr/C_SO4;
                            C22 = C1/C_SO4;
                            C23 = C15_corr/C_SO4;
                            C24 = 2*C23;
                            C25 = C5/C_H;
                            C26 = C1/C_H;

                            if  Modele_new == 0 || Modele_new == 2 || Modele_new == 3 || Modele_new == 6;
                                C11_1_corr = C11_1*Activite_AH_new(n);
                                C14_1_corr = C14_1*Activite_AH_new(n);                  
                                C19_1_corr = C19_1*Activite_AH_new(n);                

                                C27_1 = C16_1_corr/C_H;
                                C28_1 = C17_1/C_H;
                                C29_1 = C18_1/C_H;
                                C30_1 = 2*C27_1;
                                C31_1 = C5/C_A1H;
                                C32_1 = C14_1_corr/C_A1H;
                                C33_1 = C19_1_corr/C_A1H;
                                C34_1 = C16_1_corr/C_A1H;
                                C35_1 = C17_1/C_A1H;
                                C36_1 = C18_1/C_A1H;
                                C37_1 = 2*C34_1;

                                C27_1_1 = C27_1*C6_1;
                                C28_1_1 = C28_1*C7_1;
                                C29_1_1 = C29_1*C38_1;
                                C34_1_1 = C34_1*C6_1;
                                C35_1_1 = C35_1*C7_1;
                                C36_1_1 = C36_1*C38_1;

                                if  Modele_new == 0 || Modele_new == 3;
                                    C27_1_2 = C27_1*C6_2;
                                    C28_1_2 = C28_1*C7_2;
                                    C29_1_2 = C29_1*C38_2;
                                    C34_1_2 = C34_1*C6_2;
                                    C35_1_2 = C35_1*C7_2;
                                    C36_1_2 = C36_1*C38_2;
                                end

                                if  Modele_new == 0 || Modele_new == 2;
                                    C27_1_3 = C27_1*C6_3;
                                    C28_1_3 = C28_1*C7_3;
                                    C29_1_3 = C29_1*C38_3;
                                    C34_1_3 = C34_1*C6_3;
                                    C35_1_3 = C35_1*C7_3;
                                    C36_1_3 = C36_1*C38_3;
                                end
                            end

                            if  Modele_new == 0 || Modele_new == 1 || Modele_new == 3 || Modele_new == 5;
                                C11_2_corr = C11_2*Activite_AH_new(n);
                                C14_2_corr = C14_2*Activite_AH_new(n);                  
                                C19_2_corr = C19_2*Activite_AH_new(n); 

                                C27_2 = C16_2_corr/C_H;
                                C28_2 = C17_2/C_H;
                                C29_2 = C18_2/C_H;
                                C30_2 = 2*C27_2;
                                C31_2 = C5/C_A2H;
                                C32_2 = C14_2_corr/C_A2H;
                                C33_2 = C19_2_corr/C_A2H;
                                C34_2 = C16_2_corr/C_A2H;
                                C35_2 = C17_2/C_A2H;
                                C36_2 = C18_2/C_A2H;
                                C37_2 = 2*C34_2;

                                if  Modele_new == 0 || Modele_new == 3;
                                    C27_2_1 = C27_2*C6_1;
                                    C28_2_1 = C28_2*C7_1;
                                    C29_2_1 = C29_2*C38_1;
                                    C34_2_1 = C34_2*C6_1;
                                    C35_2_1 = C35_2*C7_1;
                                    C36_2_1 = C36_2*C38_1;
                                end

                                C27_2_2 = C27_2*C6_2;
                                C28_2_2 = C28_2*C7_2;
                                C29_2_2 = C29_2*C38_2;
                                C34_2_2 = C34_2*C6_2;
                                C35_2_2 = C35_2*C7_2;
                                C36_2_2 = C36_2*C38_2;

                                if  Modele_new == 0 || Modele_new == 1;
                                    C27_2_3 = C27_2*C6_3;
                                    C28_2_3 = C28_2*C7_3;
                                    C29_2_3 = C29_2*C38_3;
                                    C34_2_3 = C34_2*C6_3;
                                    C35_2_3 = C35_2*C7_3;
                                    C36_2_3 = C36_2*C38_3;
                                end
                            end

                            if  Modele_new == 0 || Modele_new == 1 || Modele_new == 2 || Modele_new == 4;
                                C11_3_corr = C11_3*Activite_AH_new(n);
                                C14_3_corr = C14_3*Activite_AH_new(n);                  
                                C19_3_corr = C19_3*Activite_AH_new(n);                    

                                C27_3 = C16_3_corr/C_H;
                                C28_3 = C17_3/C_H;
                                C29_3 = C18_3/C_H;
                                C30_3 = 2*C27_3;
                                C31_3 = C5/C_A3H;
                                C32_3 = C14_3_corr/C_A3H;
                                C33_3 = C19_3_corr/C_A3H;
                                C34_3 = C16_3_corr/C_A3H;
                                C35_3 = C17_3/C_A3H;
                                C36_3 = C18_3/C_A3H;
                                C37_3 = 2*C34_3;

                                if  Modele_new == 0 || Modele_new == 2;   
                                    C27_3_1 = C27_3*C6_1;
                                    C28_3_1 = C28_3*C7_1;
                                    C29_3_1 = C29_3*C38_1;
                                    C34_3_1 = C34_3*C6_1;
                                    C35_3_1 = C35_3*C7_1;
                                    C36_3_1 = C36_3*C38_1;
                                end

                                if  Modele_new == 0 || Modele_new == 1;  
                                    C27_3_2 = C27_3*C6_2;
                                    C28_3_2 = C28_3*C7_2;
                                    C29_3_2 = C29_3*C38_2;
                                    C34_3_2 = C34_3*C6_2;
                                    C35_3_2 = C35_3*C7_2;
                                    C36_3_2 = C36_3*C38_2;
                                end

                                C27_3_3 = C27_3*C6_3;
                                C28_3_3 = C28_3*C7_3;
                                C29_3_3 = C29_3*C38_3;
                                C34_3_3 = C34_3*C6_3;
                                C35_3_3 = C35_3*C7_3;
                                C36_3_3 = C36_3*C38_3;
                            end

                            % Initialisation de la méthode numérique
                            HSO4_stat_bef_ini = HSO4_stat_old(n);
                            if Methode_initialisation == 1
                                % composition à l'instant t-1
                                HSO4_mob_bef_ini = HSO4_mob_old(n);
                                H_mob_bef_ini = H_mob_old(n);
                                A1H_mob_bef_ini = A1H_mob_old(n);
                                A2H_mob_bef_ini = A2H_mob_old(n);
                                A3H_mob_bef_ini = A3H_mob_old(n);

                            elseif Methode_initialisation == 2
                                % composition à l'instant t si aucune réaction
                                HSO4_mob_bef_ini = (HSO4_mob_old(n) + C4*HSO4_entree)/C5;
                                H_mob_bef_ini = (H_mob_old(n) + C4*H_entree)/C5;
                                A1H_mob_bef_ini = (A1H_mob_old(n) + C4*A1H_entree)/C5;
                                A2H_mob_bef_ini = (A2H_mob_old(n) + C4*A2H_entree)/C5;
                                A3H_mob_bef_ini = (A3H_mob_old(n) + C4*A3H_entree)/C5;                    

                            elseif Methode_initialisation == 3
                                % composition en acides organiques à l'instant t en considérant [H+], [HSO4-] et qHSO4- à l'instant t-1
                                HSO4_mob_bef_ini = HSO4_mob_old(n);
                                H_mob_bef_ini = H_mob_old(n);
                                Facteur_Site_HSO4_libre_bef_ini = (N_HSO4_min + C7_1*A1H_mob_old(n) + C7_2*A2H_mob_old(n) + C7_3*A3H_mob_old(n));
                                Facteur_Site_SO4_libre_bef_ini = (N_SO4_min + C6_1*A1H_mob_old(n) + C6_2*A2H_mob_old(n) + C6_3*A3H_mob_old(n));
                                Facteur_Site_mat_libre_bef_ini = (N_mat_min + C38_1*A1H_mob_old(n) + C38_2*A2H_mob_old(n) + C38_3*A3H_mob_old(n));                    
                                A1H_mob_bef_ini = C_A1H/(C5*(1+Ka_A1H_0*Activite_AH_new(n)/H_mob_old(n))+C19_1*Activite_AH_new(n)*HSO4_stat_old(n)/(H_mob_old(n)*HSO4_mob_old(n))+C18_1/Facteur_Site_mat_libre_bef_ini+C17_1*HSO4_stat_old(n)/Facteur_Site_HSO4_libre_bef_ini+C1*Ks_A1H_SO4*SO4_stat_old(n)/Facteur_Site_SO4_libre_bef_ini);
                                A2H_mob_bef_ini = C_A2H/(C5*(1+Ka_A2H_0*Activite_AH_new(n)/H_mob_old(n))+C19_2*Activite_AH_new(n)*HSO4_stat_old(n)/(H_mob_old(n)*HSO4_mob_old(n))+C18_2/Facteur_Site_mat_libre_bef_ini+C17_2*HSO4_stat_old(n)/Facteur_Site_HSO4_libre_bef_ini+C1*Ks_A2H_SO4*SO4_stat_old(n)/Facteur_Site_SO4_libre_bef_ini);
                                A3H_mob_bef_ini = C_A3H/(C5*(1+Ka_A3H_0*Activite_AH_new(n)/H_mob_old(n))+C19_3*Activite_AH_new(n)*HSO4_stat_old(n)/(H_mob_old(n)*HSO4_mob_old(n))+C18_3/Facteur_Site_mat_libre_bef_ini+C17_3*HSO4_stat_old(n)/Facteur_Site_HSO4_libre_bef_ini+C1*Ks_A3H_SO4*SO4_stat_old(n)/Facteur_Site_SO4_libre_bef_ini);
                            end

                            if Modele_new > 0
                                if Modele_new == 1 || Modele_new == 4 || Modele_new == 5 || Modele_new == 7
                                    A1H_mob_bef_ini = 0;
                                end
                                if Modele_new == 2 || Modele_new == 4 || Modele_new == 6 || Modele_new == 7
                                    A2H_mob_bef_ini = 0;
                                end
                                if Modele_new == 3 || Modele_new == 5 || Modele_new == 6 || Modele_new == 7
                                    A3H_mob_bef_ini = 0;
                                end
                            end

                            HSO4_stat_bef = HSO4_stat_bef_ini; 
                            HSO4_mob_bef = HSO4_mob_bef_ini;
                            H_mob_bef = H_mob_bef_ini;
                            A1H_mob_bef = A1H_mob_bef_ini;
                            A2H_mob_bef = A2H_mob_bef_ini;
                            A3H_mob_bef = A3H_mob_bef_ini;   

                            % Calcul du vecteur "bilans de matière"
                            V1 = HSO4_mob_bef*H_mob_bef;
                            V2 = HSO4_stat_bef^2;
                            V4 = V2/V1;
                            V6 = HSO4_mob_bef/H_mob_bef;

                            if Modele_new < 7
                                V7 = N_SO4_min + C6_1*A1H_mob_bef + C6_2*A2H_mob_bef + C6_3*A3H_mob_bef;
                                V8 = N_HSO4_min + C7_1*A1H_mob_bef + C7_2*A2H_mob_bef + C7_3*A3H_mob_bef;
                                V9 = N_mat_min + C38_1*A1H_mob_bef + C38_2*A2H_mob_bef + C38_3*A3H_mob_bef;

                                if  Modele_new == 0 || Modele_new == 2 || Modele_new == 3 || Modele_new == 6;
                                    V3_1 = HSO4_stat_bef*A1H_mob_bef;
                                    V5_1 = V3_1/V1;
                                    V10_1 = A1H_mob_bef/V7;
                                    V11_1 = A1H_mob_bef/V9;
                                    V12_1 = V4*V10_1;
                                    V13_1 = V3_1/V8;
                                    V14_1 = A1H_mob_bef/H_mob_bef;
                                end

                                if  Modele_new == 0 || Modele_new == 1 || Modele_new == 3 || Modele_new == 5;
                                    V3_2 = HSO4_stat_bef*A2H_mob_bef;
                                    V5_2 = V3_2/V1;
                                    V10_2 = A2H_mob_bef/V7;
                                    V11_2 = A2H_mob_bef/V9;
                                    V12_2 = V4*V10_2;
                                    V13_2 = V3_2/V8;
                                    V14_2 = A2H_mob_bef/H_mob_bef;
                                end

                                if  Modele_new == 0 || Modele_new == 1 || Modele_new == 2 || Modele_new == 4;
                                    V3_3 = HSO4_stat_bef*A3H_mob_bef;
                                    V5_3 = V3_3/V1;
                                    V10_3 = A3H_mob_bef/V7;
                                    V11_3 = A3H_mob_bef/V9;
                                    V12_3 = V4*V10_3;
                                    V13_3 = V3_3/V8;
                                    V14_3 = A3H_mob_bef/H_mob_bef;
                                end
                            end

                            if Modele_new == 0                    
                                Bilan_Site_bef = C9_corr*V4 + C10*HSO4_stat_bef + C11_1_corr*V5_1 + C11_2_corr*V5_2 + C11_3_corr*V5_3 - 1;
                                Bilan_SO4_bef = C20* HSO4_mob_bef + C21*V6 + C22*HSO4_stat_bef + C23*V4 - 1;
                                Bilan_H_bef = C25*(H_mob_bef+HSO4_mob_bef+A1H_mob_bef+A2H_mob_bef+A3H_mob_bef) + C26*HSO4_stat_bef + C27_1*V12_1 + C27_2*V12_2 + C27_3*V12_3 + C28_1*V13_1 + C28_2*V13_2 + C28_3*V13_3 + C29_1*V11_1 + C29_2*V11_2 + C29_3*V11_3 - 1;             
                                Bilan_A1_bef = C31_1*A1H_mob_bef + C32_1*V14_1 + C33_1*V5_1 + C34_1*V12_1 + C35_1*V13_1 + C36_1*V11_1 - 1;               
                                Bilan_A2_bef = C31_2*A2H_mob_bef + C32_2*V14_2 + C33_2*V5_2 + C34_2*V12_2 + C35_2*V13_2 + C36_2*V11_2 - 1;
                                Bilan_A3_bef = C31_3*A3H_mob_bef + C32_3*V14_3 + C33_3*V5_3 + C34_3*V12_3 + C35_3*V13_3 + C36_3*V11_3 - 1;

                            elseif Modele_new == 1                    
                                Bilan_Site_bef = C9_corr*V4 + C10*HSO4_stat_bef + C11_2_corr*V5_2 + C11_3_corr*V5_3 - 1;
                                Bilan_SO4_bef = C20* HSO4_mob_bef + C21*V6 + C22*HSO4_stat_bef + C23*V4 - 1;
                                Bilan_H_bef = C25*(H_mob_bef+HSO4_mob_bef+A2H_mob_bef+A3H_mob_bef) + C26*HSO4_stat_bef + C27_2*V12_2 + C27_3*V12_3 + C28_2*V13_2 + C28_3*V13_3 + C29_2*V11_2 + C29_3*V11_3 - 1;             
                                Bilan_A1_bef = 0;               
                                Bilan_A2_bef = C31_2*A2H_mob_bef + C32_2*V14_2 + C33_2*V5_2 + C34_2*V12_2 + C35_2*V13_2 + C36_2*V11_2 - 1;
                                Bilan_A3_bef = C31_3*A3H_mob_bef + C32_3*V14_3 + C33_3*V5_3 + C34_3*V12_3 + C35_3*V13_3 + C36_3*V11_3 - 1;

                            elseif Modele_new == 2                    
                                Bilan_Site_bef = C9_corr*V4 + C10*HSO4_stat_bef + C11_1_corr*V5_1 + C11_3_corr*V5_3 - 1;
                                Bilan_SO4_bef = C20* HSO4_mob_bef + C21*V6 + C22*HSO4_stat_bef + C23*V4 - 1;
                                Bilan_H_bef = C25*(H_mob_bef+HSO4_mob_bef+A1H_mob_bef+A3H_mob_bef) + C26*HSO4_stat_bef + C27_1*V12_1 + C27_3*V12_3 + C28_1*V13_1 + C28_3*V13_3 + C29_1*V11_1 + C29_3*V11_3 - 1;             
                                Bilan_A1_bef = C31_1*A1H_mob_bef + C32_1*V14_1 + C33_1*V5_1 + C34_1*V12_1 + C35_1*V13_1 + C36_1*V11_1 - 1;               
                                Bilan_A2_bef = 0;
                                Bilan_A3_bef = C31_3*A3H_mob_bef + C32_3*V14_3 + C33_3*V5_3 + C34_3*V12_3 + C35_3*V13_3 + C36_3*V11_3 - 1;

                            elseif Modele_new == 3                    
                                Bilan_Site_bef = C9_corr*V4 + C10*HSO4_stat_bef + C11_1_corr*V5_1 + C11_2_corr*V5_2 - 1;
                                Bilan_SO4_bef = C20* HSO4_mob_bef + C21*V6 + C22*HSO4_stat_bef + C23*V4 - 1;
                                Bilan_H_bef = C25*(H_mob_bef+HSO4_mob_bef+A1H_mob_bef+A2H_mob_bef) + C26*HSO4_stat_bef + C27_1*V12_1 + C27_2*V12_2 + C28_1*V13_1 + C28_2*V13_2 + C29_1*V11_1 + C29_2*V11_2 - 1;             
                                Bilan_A1_bef = C31_1*A1H_mob_bef + C32_1*V14_1 + C33_1*V5_1 + C34_1*V12_1 + C35_1*V13_1 + C36_1*V11_1 - 1;               
                                Bilan_A2_bef = C31_2*A2H_mob_bef + C32_2*V14_2 + C33_2*V5_2 + C34_2*V12_2 + C35_2*V13_2 + C36_2*V11_2 - 1;
                                Bilan_A3_bef = 0;

                            elseif Modele_new == 4                    
                                Bilan_Site_bef = C9_corr*V4 + C10*HSO4_stat_bef + C11_3_corr*V5_3 - 1;
                                Bilan_SO4_bef = C20* HSO4_mob_bef + C21*V6 + C22*HSO4_stat_bef + C23*V4 - 1;
                                Bilan_H_bef = C25*(H_mob_bef+HSO4_mob_bef+A3H_mob_bef) + C26*HSO4_stat_bef + C27_3*V12_3 + C28_3*V13_3 + C29_3*V11_3 - 1;             
                                Bilan_A1_bef = 0;               
                                Bilan_A2_bef = 0;
                                Bilan_A3_bef = C31_3*A3H_mob_bef + C32_3*V14_3 + C33_3*V5_3 + C34_3*V12_3 + C35_3*V13_3 + C36_3*V11_3 - 1;

                            elseif Modele_new == 5                    
                                Bilan_Site_bef = C9_corr*V4 + C10*HSO4_stat_bef + C11_2_corr*V5_2 - 1;
                                Bilan_SO4_bef = C20* HSO4_mob_bef + C21*V6 + C22*HSO4_stat_bef + C23*V4 - 1;
                                Bilan_H_bef = C25*(H_mob_bef+HSO4_mob_bef+A2H_mob_bef) + C26*HSO4_stat_bef + C27_2*V12_2 + C28_2*V13_2 + C29_2*V11_2 - 1;             
                                Bilan_A1_bef = 0;               
                                Bilan_A2_bef = C31_2*A2H_mob_bef + C32_2*V14_2 + C33_2*V5_2 + C34_2*V12_2 + C35_2*V13_2 + C36_2*V11_2 - 1;
                                Bilan_A3_bef = 0;

                            elseif Modele_new == 6                    
                                Bilan_Site_bef = C9_corr*V4 + C10*HSO4_stat_bef + C11_1_corr*V5_1 - 1;
                                Bilan_SO4_bef = C20* HSO4_mob_bef + C21*V6 + C22*HSO4_stat_bef + C23*V4 - 1;
                                Bilan_H_bef = C25*(H_mob_bef+HSO4_mob_bef+A1H_mob_bef) + C26*HSO4_stat_bef + C27_1*V12_1 + C28_1*V13_1 + C29_1*V11_1 - 1;             
                                Bilan_A1_bef = C31_1*A1H_mob_bef + C32_1*V14_1 + C33_1*V5_1 + C34_1*V12_1 + C35_1*V13_1 + C36_1*V11_1 - 1;               
                                Bilan_A2_bef = 0;
                                Bilan_A3_bef = 0;

                            else
                                Bilan_Site_bef = C9_corr*V4 + C10*HSO4_stat_bef - 1;
                                Bilan_SO4_bef = C20* HSO4_mob_bef + C21*V6 + C22*HSO4_stat_bef + C23*V4 - 1;
                                Bilan_H_bef = C25*(H_mob_bef+HSO4_mob_bef) + C26*HSO4_stat_bef - 1;
                                Bilan_A1_bef = 0;               
                                Bilan_A2_bef = 0;
                                Bilan_A3_bef = 0;
                            end

                            Bilans_bef = [Bilan_Site_bef; Bilan_SO4_bef; Bilan_H_bef; Bilan_A1_bef; Bilan_A2_bef; Bilan_A3_bef];
                            N_Bilans_bef = norm(Bilans_bef);
                        end

                        % Calcul de la matrice jacobienne
                        V15 = HSO4_stat_bef/V1;
                        V16 = V4/H_mob_bef;           
                        V17 = V4/HSO4_mob_bef;
                        V31 = V6/H_mob_bef;

                        if Modele_new < 7
                            V21 = V7^2;
                            V22 = V8^2;
                            V23 = V9^2;
                            V32 = V4/V7;
                            V33 = HSO4_stat_bef/V8;

                            if  Modele_new == 0 || Modele_new == 2 || Modele_new == 3 || Modele_new == 6;
                                V18_1 = V5_1/H_mob_bef;
                                V19_1 = V5_1/HSO4_mob_bef;
                                V20_1 = A1H_mob_bef/V1;
                                V21_1 = A1H_mob_bef/V21;
                                V22_1 = A1H_mob_bef/V22;
                                V23_1 = A1H_mob_bef/V23;
                                V24_1 = V4*V21_1;
                                V25_1 = HSO4_stat_bef*V22_1;
                                V26_1 = V12_1/H_mob_bef;
                                V27_1 = V12_1/HSO4_mob_bef;
                                V28_1 = V5_1/V7;
                                V29_1 = A1H_mob_bef/V8;
                                V30_1 = V14_1/H_mob_bef;
                            end

                            if  Modele_new == 0 || Modele_new == 1 || Modele_new == 3 || Modele_new == 5;
                                V18_2 = V5_2/H_mob_bef;
                                V19_2 = V5_2/HSO4_mob_bef;
                                V20_2 = A2H_mob_bef/V1;
                                V21_2 = A2H_mob_bef/V21;
                                V22_2 = A2H_mob_bef/V22;
                                V23_2 = A2H_mob_bef/V23;
                                V24_2 = V4*V21_2;
                                V25_2 = HSO4_stat_bef*V22_2;
                                V26_2 = V12_2/H_mob_bef;
                                V27_2 = V12_2/HSO4_mob_bef;
                                V28_2 = V5_2/V7;
                                V29_2 = A2H_mob_bef/V8;
                                V30_2 = V14_2/H_mob_bef;
                            end

                            if  Modele_new == 0 || Modele_new == 1 || Modele_new == 2 || Modele_new == 4;
                                V18_3 = V5_3/H_mob_bef;
                                V19_3 = V5_3/HSO4_mob_bef;
                                V20_3 = A3H_mob_bef/V1;
                                V21_3 = A3H_mob_bef/V21;
                                V22_3 = A3H_mob_bef/V22;
                                V23_3 = A3H_mob_bef/V23;
                                V24_3 = V4*V21_3;
                                V25_3 = HSO4_stat_bef*V22_3;
                                V26_3 = V12_3/H_mob_bef;
                                V27_3 = V12_3/HSO4_mob_bef;
                                V28_3 = V5_3/V7;                    
                                V29_3 = A3H_mob_bef/V8;
                                V30_3 = V14_3/H_mob_bef;
                            end
                        end

                        if Modele_new == 0
                            Jac_B_Site_HSO4_stat = C12_corr*V15 + C11_1_corr*V20_1 + C11_2_corr*V20_2 + C11_3_corr*V20_3 + C10;
                            Jac_B_Site_HSO4_mob = -(C9_corr*V17 + C11_1_corr*V19_1 + C11_2_corr*V19_2 + C11_3_corr*V19_3);
                            Jac_B_Site_H_mob = -(C9_corr*V16 + C11_1_corr*V18_1 + C11_2_corr*V18_2 + C11_3_corr*V18_3);
                            Jac_B_Site_A1H_mob = C11_1_corr*V15;
                            Jac_B_Site_A2H_mob = C11_2_corr*V15;
                            Jac_B_Site_A3H_mob = C11_3_corr*V15;

                            Jac_B_SO4_HSO4_stat =  C24*V15 + C22;
                            Jac_B_SO4_HSO4_mob =  C20 + C21/H_mob_bef - C23*V17;
                            Jac_B_SO4_H_mob = -(C21*V31 + C23*V16);                
                            Jac_B_SO4_A1H_mob = 0;
                            Jac_B_SO4_A2H_mob = 0;
                            Jac_B_SO4_A3H_mob = 0;

                            Jac_B_H_HSO4_stat = C26 + C30_1*V28_1 + C30_2*V28_2 + C30_3*V28_3 + C28_1*V29_1 + C28_2*V29_2 + C28_3*V29_3;
                            Jac_B_H_HSO4_mob = C25 - (C27_1*V27_1 + C27_2*V27_2 + C27_3*V27_3);
                            Jac_B_H_H_mob = C25 - (C27_1*V26_1 + C27_2*V26_2 + C27_3*V26_3);
                            Jac_B_H_A1H_mob = C25 + C27_1*V32 - C27_1_1*V24_1 - C27_2_1*V24_2 - C27_3_1*V24_3 + C28_1*V33 - C28_1_1*V25_1 - C28_2_1*V25_2 - C28_3_1*V25_3 + C29_1/V9 - C29_1_1*V23_1 - C29_2_1*V23_2 - C29_3_1*V23_3;
                            Jac_B_H_A2H_mob = C25 + C27_2*V32 - C27_1_2*V24_1 - C27_2_2*V24_2 - C27_3_2*V24_3 + C28_2*V33 - C28_1_2*V25_1 - C28_2_2*V25_2 - C28_3_2*V25_3 + C29_2/V9 - C29_1_2*V23_1 - C29_2_2*V23_2 - C29_3_2*V23_3;
                            Jac_B_H_A3H_mob = C25 + C27_3*V32 - C27_1_3*V24_1 - C27_2_3*V24_2 - C27_3_3*V24_3 + C28_3*V33 - C28_1_3*V25_1 - C28_2_3*V25_2 - C28_3_3*V25_3 + C29_3/V9 - C29_1_3*V23_1 - C29_2_3*V23_2 - C29_3_3*V23_3;

                            Jac_B_A1_HSO4_stat = C33_1*V20_1 + C37_1*V28_1 + C35_1*V29_1;
                            Jac_B_A1_HSO4_mob = -(C33_1*V19_1 + C34_1*V27_1);
                            Jac_B_A1_H_mob = -(C32_1*V30_1 + C33_1*V18_1 + C34_1*V26_1);                
                            Jac_B_A1_A1H_mob = C31_1 + C32_1/H_mob_bef + C33_1*V15 + C34_1*V32 - C34_1_1*V24_1 + C35_1*V33 - C35_1_1*V25_1 + C36_1/V9 - C36_1_1*V23_1;
                            Jac_B_A1_A2H_mob = -(C34_1_2*V24_1 + C35_1_2*V25_1 + C36_1_2*V23_1);
                            Jac_B_A1_A3H_mob = -(C34_1_3*V24_1 + C35_1_3*V25_1 + C36_1_3*V23_1);

                            Jac_B_A2_HSO4_stat = C33_2*V20_2 + C37_2*V28_2 + C35_2*V29_2;
                            Jac_B_A2_HSO4_mob = -(C33_2*V19_2 + C34_2*V27_2);
                            Jac_B_A2_H_mob = -(C32_2*V30_2 + C33_2*V18_2 + C34_2*V26_2); 
                            Jac_B_A2_A1H_mob = -(C34_2_1*V24_2 + C35_2_1*V25_2 + C36_2_1*V23_2);
                            Jac_B_A2_A2H_mob = C31_2 + C32_2/H_mob_bef + C33_2*V15 + C34_2*V32 - C34_2_2*V24_2 + C35_2*V33 - C35_2_2*V25_2 + C36_2/V9 - C36_2_2*V23_2;
                            Jac_B_A2_A3H_mob = -(C34_2_3*V24_2 + C35_2_3*V25_2 + C36_2_3*V23_2);

                            Jac_B_A3_HSO4_stat = C33_3*V20_3 + C37_3*V28_3 + C35_3*V29_3;
                            Jac_B_A3_HSO4_mob = -(C33_3*V19_3 + C34_3*V27_3);
                            Jac_B_A3_H_mob = -(C32_3*V30_3 + C33_3*V18_3 + C34_3*V26_3);               
                            Jac_B_A3_A1H_mob = -(C34_3_1*V24_3 + C35_3_1*V25_3 + C36_3_1*V23_3);
                            Jac_B_A3_A2H_mob = -(C34_3_2*V24_3 + C35_3_2*V25_3 + C36_3_2*V23_3);
                            Jac_B_A3_A3H_mob = C31_3 + C32_3/H_mob_bef + C33_3*V15 + C34_3*V32 - C34_3_3*V24_3 + C35_3*V33 - C35_3_3*V25_3 + C36_3/V9 - C36_3_3*V23_3;

                            Jac_6Bilans = [Jac_B_Site_HSO4_stat Jac_B_Site_HSO4_mob Jac_B_Site_H_mob Jac_B_Site_A1H_mob Jac_B_Site_A2H_mob Jac_B_Site_A3H_mob;
                                Jac_B_SO4_HSO4_stat Jac_B_SO4_HSO4_mob Jac_B_SO4_H_mob Jac_B_SO4_A1H_mob Jac_B_SO4_A2H_mob Jac_B_SO4_A3H_mob;
                                Jac_B_H_HSO4_stat Jac_B_H_HSO4_mob Jac_B_H_H_mob Jac_B_H_A1H_mob Jac_B_H_A2H_mob Jac_B_H_A3H_mob;
                                Jac_B_A1_HSO4_stat Jac_B_A1_HSO4_mob Jac_B_A1_H_mob Jac_B_A1_A1H_mob Jac_B_A1_A2H_mob Jac_B_A1_A3H_mob;         
                                Jac_B_A2_HSO4_stat Jac_B_A2_HSO4_mob Jac_B_A2_H_mob Jac_B_A2_A1H_mob Jac_B_A2_A2H_mob Jac_B_A2_A3H_mob;
                                Jac_B_A3_HSO4_stat Jac_B_A3_HSO4_mob Jac_B_A3_H_mob Jac_B_A3_A1H_mob Jac_B_A3_A2H_mob Jac_B_A3_A3H_mob];

                            % Calcul des variations à appliquer aux solutions
                            Variations6 = -(Jac_6Bilans\Bilans_bef);
                            HSO4_stat_variation = Variations6(1);
                            HSO4_mob_variation = Variations6(2);
                            H_mob_variation = Variations6(3);
                            A1H_mob_variation = Variations6(4);
                            A2H_mob_variation = Variations6(5);
                            A3H_mob_variation = Variations6(6);

                        elseif Modele_new == 1                    
                            Jac_B_Site_HSO4_stat = C12_corr*V15 + C11_2_corr*V20_2 + C11_3_corr*V20_3 + C10;
                            Jac_B_Site_HSO4_mob = -(C9_corr*V17 + C11_2_corr*V19_2 + C11_3_corr*V19_3);
                            Jac_B_Site_H_mob = -(C9_corr*V16 + C11_2_corr*V18_2 + C11_3_corr*V18_3);
                            Jac_B_Site_A2H_mob = C11_2_corr*V15;
                            Jac_B_Site_A3H_mob = C11_3_corr*V15;

                            Jac_B_SO4_HSO4_stat =  C24*V15 + C22;
                            Jac_B_SO4_HSO4_mob =  C20 + C21/H_mob_bef - C23*V17;
                            Jac_B_SO4_H_mob = -(C21*V31 + C23*V16);
                            Jac_B_SO4_A2H_mob = 0;
                            Jac_B_SO4_A3H_mob = 0;

                            Jac_B_H_HSO4_stat = C26 + C30_2*V28_2 + C30_3*V28_3 + C28_2*V29_2 + C28_3*V29_3;
                            Jac_B_H_HSO4_mob = C25 - (C27_2*V27_2 + C27_3*V27_3);
                            Jac_B_H_H_mob = C25 - (C27_2*V26_2 + C27_3*V26_3);
                            Jac_B_H_A2H_mob = C25 + C27_2*V32 - C27_2_2*V24_2 - C27_3_2*V24_3 + C28_2*V33 - C28_2_2*V25_2 - C28_3_2*V25_3 + C29_2/V9 - C29_2_2*V23_2 - C29_3_2*V23_3;
                            Jac_B_H_A3H_mob = C25 + C27_3*V32 - C27_2_3*V24_2 - C27_3_3*V24_3 + C28_3*V33 - C28_2_3*V25_2 - C28_3_3*V25_3 + C29_3/V9 - C29_2_3*V23_2 - C29_3_3*V23_3;

                            Jac_B_A2_HSO4_stat = C33_2*V20_2 + C37_2*V28_2 + C35_2*V29_2;
                            Jac_B_A2_HSO4_mob = -(C33_2*V19_2 + C34_2*V27_2);
                            Jac_B_A2_H_mob = -(C32_2*V30_2 + C33_2*V18_2 + C34_2*V26_2); 
                            Jac_B_A2_A2H_mob = C31_2 + C32_2/H_mob_bef + C33_2*V15 + C34_2*V32 - C34_2_2*V24_2 + C35_2*V33 - C35_2_2*V25_2 + C36_2/V9 - C36_2_2*V23_2;
                            Jac_B_A2_A3H_mob = -(C34_2_3*V24_2 + C35_2_3*V25_2 + C36_2_3*V23_2);

                            Jac_B_A3_HSO4_stat = C33_3*V20_3 + C37_3*V28_3 + C35_3*V29_3;
                            Jac_B_A3_HSO4_mob = -(C33_3*V19_3 + C34_3*V27_3);
                            Jac_B_A3_H_mob = -(C32_3*V30_3 + C33_3*V18_3 + C34_3*V26_3);               
                            Jac_B_A3_A2H_mob = -(C34_3_2*V24_3 + C35_3_2*V25_3 + C36_3_2*V23_3);
                            Jac_B_A3_A3H_mob = C31_3 + C32_3/H_mob_bef + C33_3*V15 + C34_3*V32 - C34_3_3*V24_3 + C35_3*V33 - C35_3_3*V25_3 + C36_3/V9 - C36_3_3*V23_3;

                            Jac_5Bilans = [Jac_B_Site_HSO4_stat Jac_B_Site_HSO4_mob Jac_B_Site_H_mob Jac_B_Site_A2H_mob Jac_B_Site_A3H_mob;
                                Jac_B_SO4_HSO4_stat Jac_B_SO4_HSO4_mob Jac_B_SO4_H_mob Jac_B_SO4_A2H_mob Jac_B_SO4_A3H_mob;
                                Jac_B_H_HSO4_stat Jac_B_H_HSO4_mob Jac_B_H_H_mob Jac_B_H_A2H_mob Jac_B_H_A3H_mob;
                                Jac_B_A2_HSO4_stat Jac_B_A2_HSO4_mob Jac_B_A2_H_mob Jac_B_A2_A2H_mob Jac_B_A2_A3H_mob;
                                Jac_B_A3_HSO4_stat Jac_B_A3_HSO4_mob Jac_B_A3_H_mob Jac_B_A3_A2H_mob Jac_B_A3_A3H_mob];

                            % Calcul des variations à appliquer aux solutions
                            Variations5 = -(Jac_5Bilans\Bilans_bef([1:3,5,6]));
                            HSO4_stat_variation = Variations5(1);
                            HSO4_mob_variation = Variations5(2);
                            H_mob_variation = Variations5(3);
                            A1H_mob_variation = 0;
                            A2H_mob_variation = Variations5(4);
                            A3H_mob_variation = Variations5(5);

                        elseif Modele_new == 2                    
                            Jac_B_Site_HSO4_stat = C12_corr*V15 + C11_1_corr*V20_1 + C11_3_corr*V20_3 + C10;
                            Jac_B_Site_HSO4_mob = -(C9_corr*V17 + C11_1_corr*V19_1 + C11_3_corr*V19_3);
                            Jac_B_Site_H_mob = -(C9_corr*V16 + C11_1_corr*V18_1 + C11_3_corr*V18_3);
                            Jac_B_Site_A1H_mob = C11_1_corr*V15;
                            Jac_B_Site_A3H_mob = C11_3_corr*V15;

                            Jac_B_SO4_HSO4_stat =  C24*V15 + C22;
                            Jac_B_SO4_HSO4_mob =  C20 + C21/H_mob_bef - C23*V17;
                            Jac_B_SO4_H_mob = -(C21*V31 + C23*V16);                
                            Jac_B_SO4_A1H_mob = 0;
                            Jac_B_SO4_A3H_mob = 0;

                            Jac_B_H_HSO4_stat = C26 + C30_1*V28_1 + C30_3*V28_3 + C28_1*V29_1 + C28_3*V29_3;
                            Jac_B_H_HSO4_mob = C25 - (C27_1*V27_1 + C27_3*V27_3);
                            Jac_B_H_H_mob = C25 - (C27_1*V26_1 + C27_3*V26_3);
                            Jac_B_H_A1H_mob = C25 + C27_1*V32 - C27_1_1*V24_1 - C27_3_1*V24_3 + C28_1*V33 - C28_1_1*V25_1 - C28_3_1*V25_3 + C29_1/V9 - C29_1_1*V23_1 - C29_3_1*V23_3;
                            Jac_B_H_A3H_mob = C25 + C27_3*V32 - C27_1_3*V24_1 - C27_3_3*V24_3 + C28_3*V33 - C28_1_3*V25_1 - C28_3_3*V25_3 + C29_3/V9 - C29_1_3*V23_1 - C29_3_3*V23_3;

                            Jac_B_A1_HSO4_stat = C33_1*V20_1 + C37_1*V28_1 + C35_1*V29_1;
                            Jac_B_A1_HSO4_mob = -(C33_1*V19_1 + C34_1*V27_1);
                            Jac_B_A1_H_mob = -(C32_1*V30_1 + C33_1*V18_1 + C34_1*V26_1);                
                            Jac_B_A1_A1H_mob = C31_1 + C32_1/H_mob_bef + C33_1*V15 + C34_1*V32 - C34_1_1*V24_1 + C35_1*V33 - C35_1_1*V25_1 + C36_1/V9 - C36_1_1*V23_1;
                            Jac_B_A1_A3H_mob = -(C34_1_3*V24_1 + C35_1_3*V25_1 + C36_1_3*V23_1);

                            Jac_B_A3_HSO4_stat = C33_3*V20_3 + C37_3*V28_3 + C35_3*V29_3;
                            Jac_B_A3_HSO4_mob = -(C33_3*V19_3 + C34_3*V27_3);
                            Jac_B_A3_H_mob = -(C32_3*V30_3 + C33_3*V18_3 + C34_3*V26_3);               
                            Jac_B_A3_A1H_mob = -(C34_3_1*V24_3 + C35_3_1*V25_3 + C36_3_1*V23_3);
                            Jac_B_A3_A3H_mob = C31_3 + C32_3/H_mob_bef + C33_3*V15 + C34_3*V32 - C34_3_3*V24_3 + C35_3*V33 - C35_3_3*V25_3 + C36_3/V9 - C36_3_3*V23_3;

                            Jac_5Bilans = [Jac_B_Site_HSO4_stat Jac_B_Site_HSO4_mob Jac_B_Site_H_mob Jac_B_Site_A1H_mob Jac_B_Site_A3H_mob;
                                Jac_B_SO4_HSO4_stat Jac_B_SO4_HSO4_mob Jac_B_SO4_H_mob Jac_B_SO4_A1H_mob Jac_B_SO4_A3H_mob;
                                Jac_B_H_HSO4_stat Jac_B_H_HSO4_mob Jac_B_H_H_mob Jac_B_H_A1H_mob Jac_B_H_A3H_mob;
                                Jac_B_A1_HSO4_stat Jac_B_A1_HSO4_mob Jac_B_A1_H_mob Jac_B_A1_A1H_mob Jac_B_A1_A3H_mob;
                                Jac_B_A3_HSO4_stat Jac_B_A3_HSO4_mob Jac_B_A3_H_mob Jac_B_A3_A1H_mob Jac_B_A3_A3H_mob];

                            % Calcul des variations à appliquer aux solutions
                            Variations5 = -(Jac_5Bilans\Bilans_bef([1:4,6]));
                            HSO4_stat_variation = Variations5(1);
                            HSO4_mob_variation = Variations5(2);
                            H_mob_variation = Variations5(3);
                            A1H_mob_variation = Variations5(4);
                            A2H_mob_variation = 0;
                            A3H_mob_variation = Variations5(5);

                        elseif Modele_new == 3                    
                            Jac_B_Site_HSO4_stat = C12_corr*V15 + C11_1_corr*V20_1 + C11_2_corr*V20_2 + C10;
                            Jac_B_Site_HSO4_mob = -(C9_corr*V17 + C11_1_corr*V19_1 + C11_2_corr*V19_2);
                            Jac_B_Site_H_mob = -(C9_corr*V16 + C11_1_corr*V18_1 + C11_2_corr*V18_2);
                            Jac_B_Site_A1H_mob = C11_1_corr*V15;
                            Jac_B_Site_A2H_mob = C11_2_corr*V15;

                            Jac_B_SO4_HSO4_stat =  C24*V15 + C22;
                            Jac_B_SO4_HSO4_mob =  C20 + C21/H_mob_bef - C23*V17;
                            Jac_B_SO4_H_mob = -(C21*V31 + C23*V16);                
                            Jac_B_SO4_A1H_mob = 0;
                            Jac_B_SO4_A2H_mob = 0;

                            Jac_B_H_HSO4_stat = C26 + C30_1*V28_1 + C30_2*V28_2 + C28_1*V29_1 + C28_2*V29_2;
                            Jac_B_H_HSO4_mob = C25 - (C27_1*V27_1 + C27_2*V27_2);
                            Jac_B_H_H_mob = C25 - (C27_1*V26_1 + C27_2*V26_2);
                            Jac_B_H_A1H_mob = C25 + C27_1*V32 - C27_1_1*V24_1 - C27_2_1*V24_2 + C28_1*V33 - C28_1_1*V25_1 - C28_2_1*V25_2 + C29_1/V9 - C29_1_1*V23_1 - C29_2_1*V23_2;
                            Jac_B_H_A2H_mob = C25 + C27_2*V32 - C27_1_2*V24_1 - C27_2_2*V24_2 + C28_2*V33 - C28_1_2*V25_1 - C28_2_2*V25_2 + C29_2/V9 - C29_1_2*V23_1 - C29_2_2*V23_2;

                            Jac_B_A1_HSO4_stat = C33_1*V20_1 + C37_1*V28_1 + C35_1*V29_1;
                            Jac_B_A1_HSO4_mob = -(C33_1*V19_1 + C34_1*V27_1);
                            Jac_B_A1_H_mob = -(C32_1*V30_1 + C33_1*V18_1 + C34_1*V26_1);                
                            Jac_B_A1_A1H_mob = C31_1 + C32_1/H_mob_bef + C33_1*V15 + C34_1*V32 - C34_1_1*V24_1 + C35_1*V33 - C35_1_1*V25_1 + C36_1/V9 - C36_1_1*V23_1;
                            Jac_B_A1_A2H_mob = -(C34_1_2*V24_1 + C35_1_2*V25_1 + C36_1_2*V23_1);

                            Jac_B_A2_HSO4_stat = C33_2*V20_2 + C37_2*V28_2 + C35_2*V29_2;
                            Jac_B_A2_HSO4_mob = -(C33_2*V19_2 + C34_2*V27_2);
                            Jac_B_A2_H_mob = -(C32_2*V30_2 + C33_2*V18_2 + C34_2*V26_2); 
                            Jac_B_A2_A1H_mob = -(C34_2_1*V24_2 + C35_2_1*V25_2 + C36_2_1*V23_2);
                            Jac_B_A2_A2H_mob = C31_2 + C32_2/H_mob_bef + C33_2*V15 + C34_2*V32 - C34_2_2*V24_2 + C35_2*V33 - C35_2_2*V25_2 + C36_2/V9 - C36_2_2*V23_2;

                            Jac_5Bilans = [Jac_B_Site_HSO4_stat Jac_B_Site_HSO4_mob Jac_B_Site_H_mob Jac_B_Site_A1H_mob Jac_B_Site_A2H_mob;
                                Jac_B_SO4_HSO4_stat Jac_B_SO4_HSO4_mob Jac_B_SO4_H_mob Jac_B_SO4_A1H_mob Jac_B_SO4_A2H_mob;
                                Jac_B_H_HSO4_stat Jac_B_H_HSO4_mob Jac_B_H_H_mob Jac_B_H_A1H_mob Jac_B_H_A2H_mob;
                                Jac_B_A1_HSO4_stat Jac_B_A1_HSO4_mob Jac_B_A1_H_mob Jac_B_A1_A1H_mob Jac_B_A1_A2H_mob;         
                                Jac_B_A2_HSO4_stat Jac_B_A2_HSO4_mob Jac_B_A2_H_mob Jac_B_A2_A1H_mob Jac_B_A2_A2H_mob];

                            % Calcul des variations à appliquer aux solutions
                            Variations5 = -(Jac_5Bilans\Bilans_bef(1:5));
                            HSO4_stat_variation = Variations5(1);
                            HSO4_mob_variation = Variations5(2);
                            H_mob_variation = Variations5(3);
                            A1H_mob_variation = Variations5(4);
                            A2H_mob_variation = Variations5(5);
                            A3H_mob_variation = 0;

                        elseif Modele_new == 4                    
                            Jac_B_Site_HSO4_stat = C12_corr*V15 + C11_3_corr*V20_3 + C10;
                            Jac_B_Site_HSO4_mob = -(C9_corr*V17 + C11_3_corr*V19_3);
                            Jac_B_Site_H_mob = -(C9_corr*V16 + C11_3_corr*V18_3);
                            Jac_B_Site_A3H_mob = C11_3_corr*V15;

                            Jac_B_SO4_HSO4_stat =  C24*V15 + C22;
                            Jac_B_SO4_HSO4_mob =  C20 + C21/H_mob_bef - C23*V17;
                            Jac_B_SO4_H_mob = -(C21*V31 + C23*V16);                
                            Jac_B_SO4_A3H_mob = 0;

                            Jac_B_H_HSO4_stat = C26 + C30_3*V28_3 + C28_3*V29_3;
                            Jac_B_H_HSO4_mob = C25 - C27_3*V27_3;
                            Jac_B_H_H_mob = C25 - C27_3*V26_3;
                            Jac_B_H_A3H_mob = C25 + C27_3*V32 - C27_3_3*V24_3 + C28_3*V33 - C28_3_3*V25_3 + C29_3/V9 - C29_3_3*V23_3;

                            Jac_B_A3_HSO4_stat = C33_3*V20_3 + C37_3*V28_3 + C35_3*V29_3;
                            Jac_B_A3_HSO4_mob = -(C33_3*V19_3 + C34_3*V27_3);
                            Jac_B_A3_H_mob = -(C32_3*V30_3 + C33_3*V18_3 + C34_3*V26_3);               
                            Jac_B_A3_A3H_mob = C31_3 + C32_3/H_mob_bef + C33_3*V15 + C34_3*V32 - C34_3_3*V24_3 + C35_3*V33 - C35_3_3*V25_3 + C36_3/V9 - C36_3_3*V23_3;

                            Jac_4Bilans = [Jac_B_Site_HSO4_stat Jac_B_Site_HSO4_mob Jac_B_Site_H_mob Jac_B_Site_A3H_mob;
                                Jac_B_SO4_HSO4_stat Jac_B_SO4_HSO4_mob Jac_B_SO4_H_mob Jac_B_SO4_A3H_mob;
                                Jac_B_H_HSO4_stat Jac_B_H_HSO4_mob Jac_B_H_H_mob Jac_B_H_A3H_mob;
                                Jac_B_A3_HSO4_stat Jac_B_A3_HSO4_mob Jac_B_A3_H_mob Jac_B_A3_A3H_mob];

                            % Calcul des variations à appliquer aux solutions
                            Variations4 = -(Jac_4Bilans\Bilans_bef([1:3,6]));
                            HSO4_stat_variation = Variations4(1);
                            HSO4_mob_variation = Variations4(2);
                            H_mob_variation = Variations4(3);
                            A1H_mob_variation = 0;
                            A2H_mob_variation = 0;
                            A3H_mob_variation = Variations4(4);

                        elseif Modele_new == 5                    
                            Jac_B_Site_HSO4_stat = C12_corr*V15 + C11_2_corr*V20_2 + C10;
                            Jac_B_Site_HSO4_mob = -(C9_corr*V17 + C11_2_corr*V19_2);
                            Jac_B_Site_H_mob = -(C9_corr*V16 + C11_2_corr*V18_2);
                            Jac_B_Site_A2H_mob = C11_2_corr*V15;

                            Jac_B_SO4_HSO4_stat =  C24*V15 + C22;
                            Jac_B_SO4_HSO4_mob =  C20 + C21/H_mob_bef - C23*V17;
                            Jac_B_SO4_H_mob = -(C21*V31 + C23*V16);
                            Jac_B_SO4_A2H_mob = 0;

                            Jac_B_H_HSO4_stat = C26 + C30_2*V28_2 + C28_2*V29_2;
                            Jac_B_H_HSO4_mob = C25 - C27_2*V27_2;
                            Jac_B_H_H_mob = C25 - C27_2*V26_2;
                            Jac_B_H_A2H_mob = C25 + C27_2*V32 - C27_2_2*V24_2 + C28_2*V33 - C28_2_2*V25_2 + C29_2/V9 - C29_2_2*V23_2;

                            Jac_B_A2_HSO4_stat = C33_2*V20_2 + C37_2*V28_2 + C35_2*V29_2;
                            Jac_B_A2_HSO4_mob = -(C33_2*V19_2 + C34_2*V27_2);
                            Jac_B_A2_H_mob = -(C32_2*V30_2 + C33_2*V18_2 + C34_2*V26_2);
                            Jac_B_A2_A2H_mob = C31_2 + C32_2/H_mob_bef + C33_2*V15 + C34_2*V32 - C34_2_2*V24_2 + C35_2*V33 - C35_2_2*V25_2 + C36_2/V9 - C36_2_2*V23_2;

                            Jac_4Bilans = [Jac_B_Site_HSO4_stat Jac_B_Site_HSO4_mob Jac_B_Site_H_mob Jac_B_Site_A2H_mob;
                                Jac_B_SO4_HSO4_stat Jac_B_SO4_HSO4_mob Jac_B_SO4_H_mob Jac_B_SO4_A2H_mob;
                                Jac_B_H_HSO4_stat Jac_B_H_HSO4_mob Jac_B_H_H_mob Jac_B_H_A2H_mob;        
                                Jac_B_A2_HSO4_stat Jac_B_A2_HSO4_mob Jac_B_A2_H_mob Jac_B_A2_A2H_mob];

                            % Calcul des variations à appliquer aux solutions
                            Variations4 = -(Jac_4Bilans\Bilans_bef([1:3,5]));
                            HSO4_stat_variation = Variations4(1);
                            HSO4_mob_variation = Variations4(2);
                            H_mob_variation = Variations4(3);
                            A1H_mob_variation = 0;
                            A2H_mob_variation = Variations4(4);
                            A3H_mob_variation = 0;

                        elseif Modele_new == 6                    
                            Jac_B_Site_HSO4_stat = C12_corr*V15 + C11_1_corr*V20_1 + C10;
                            Jac_B_Site_HSO4_mob = -(C9_corr*V17 + C11_1_corr*V19_1);
                            Jac_B_Site_H_mob = -(C9_corr*V16 + C11_1_corr*V18_1);
                            Jac_B_Site_A1H_mob = C11_1_corr*V15;

                            Jac_B_SO4_HSO4_stat =  C24*V15 + C22;
                            Jac_B_SO4_HSO4_mob =  C20 + C21/H_mob_bef - C23*V17;
                            Jac_B_SO4_H_mob = -(C21*V31 + C23*V16);                
                            Jac_B_SO4_A1H_mob = 0;

                            Jac_B_H_HSO4_stat = C26 + C30_1*V28_1 + C28_1*V29_1;
                            Jac_B_H_HSO4_mob = C25 - C27_1*V27_1;
                            Jac_B_H_H_mob = C25 - C27_1*V26_1;
                            Jac_B_H_A1H_mob = C25 + C27_1*V32 - C27_1_1*V24_1 + C28_1*V33 - C28_1_1*V25_1 + C29_1/V9 - C29_1_1*V23_1;

                            Jac_B_A1_HSO4_stat = C33_1*V20_1 + C37_1*V28_1 + C35_1*V29_1;
                            Jac_B_A1_HSO4_mob = -(C33_1*V19_1 + C34_1*V27_1);
                            Jac_B_A1_H_mob = -(C32_1*V30_1 + C33_1*V18_1 + C34_1*V26_1);                
                            Jac_B_A1_A1H_mob = C31_1 + C32_1/H_mob_bef + C33_1*V15 + C34_1*V32 - C34_1_1*V24_1 + C35_1*V33 - C35_1_1*V25_1 + C36_1/V9 - C36_1_1*V23_1;

                            Jac_4Bilans = [Jac_B_Site_HSO4_stat Jac_B_Site_HSO4_mob Jac_B_Site_H_mob Jac_B_Site_A1H_mob;
                                Jac_B_SO4_HSO4_stat Jac_B_SO4_HSO4_mob Jac_B_SO4_H_mob Jac_B_SO4_A1H_mob;
                                Jac_B_H_HSO4_stat Jac_B_H_HSO4_mob Jac_B_H_H_mob Jac_B_H_A1H_mob;
                                Jac_B_A1_HSO4_stat Jac_B_A1_HSO4_mob Jac_B_A1_H_mob Jac_B_A1_A1H_mob];

                            % Calcul des variations à appliquer aux solutions
                            Variations4 = -(Jac_4Bilans\Bilans_bef(1:4));
                            HSO4_stat_variation = Variations4(1);
                            HSO4_mob_variation = Variations4(2);
                            H_mob_variation = Variations4(3);
                            A1H_mob_variation = Variations4(4);
                            A2H_mob_variation = 0;
                            A3H_mob_variation = 0;

                        else
                            Jac_B_Site_HSO4_stat = C12_corr*V15 + C10;
                            Jac_B_Site_HSO4_mob = -C9_corr*V17;
                            Jac_B_Site_H_mob = -C9_corr*V16;

                            Jac_B_SO4_HSO4_stat =  C24*V15 + C22;
                            Jac_B_SO4_HSO4_mob =  C20 + C21/H_mob_bef - C23*V17;
                            Jac_B_SO4_H_mob = -(C21*V31 + C23*V16);                

                            Jac_B_H_HSO4_stat = C26;
                            Jac_B_H_HSO4_mob = C25;
                            Jac_B_H_H_mob = C25;            

                            Jac_3Bilans = [Jac_B_Site_HSO4_stat Jac_B_Site_HSO4_mob Jac_B_Site_H_mob;
                                Jac_B_SO4_HSO4_stat Jac_B_SO4_HSO4_mob Jac_B_SO4_H_mob;
                                Jac_B_H_HSO4_stat Jac_B_H_HSO4_mob Jac_B_H_H_mob];

                            % Calcul des variations à appliquer aux solutions
                            Variations3 = -(Jac_3Bilans\Bilans_bef(1:3));
                            HSO4_stat_variation = Variations3(1);
                            HSO4_mob_variation = Variations3(2);
                            H_mob_variation = Variations3(3);
                            A1H_mob_variation = 0;
                            A2H_mob_variation = 0;
                            A3H_mob_variation = 0;
                        end

                        % Relaxation/Accélération de la méthode si nouvelle solution pire/meilleure
                            % Calcul du facteur de relaxation min et d'accélération max
                            Relaxation_min = 1;
                            Acceleration_max = Variation_rel_max*max([HSO4_stat_bef,HSO4_mob_bef,H_mob_bef,A1H_mob_bef,A2H_mob_bef,A3H_mob_bef])/min([max([0,abs(HSO4_stat_variation)]),max([0,abs(HSO4_mob_variation)]),max([0,abs(H_mob_variation)]),max([0,abs(A1H_mob_variation)]),max([0,abs(A2H_mob_variation)]),max([0,abs(A3H_mob_variation)])]);

                            HSO4_stat_variation_ini = HSO4_stat_variation;
                            if HSO4_stat_variation < 0
                                if HSO4_stat_bef > 0
                                    Relaxation_min = max(Relaxation_min, abs(HSO4_stat_variation) / HSO4_stat_bef);
                                    Acceleration_max = min(Acceleration_max, HSO4_stat_bef/abs(HSO4_stat_variation));
                                else
                                    HSO4_stat_variation_ini = 0;
                                end
                            elseif HSO4_stat_variation > 0
                                Acceleration_max = min(Acceleration_max, Variation_rel_max*HSO4_stat_bef/abs(HSO4_stat_variation));                    
                            end

                            HSO4_mob_variation_ini = HSO4_mob_variation;
                            if HSO4_mob_variation < 0
                                if HSO4_mob_bef > 0
                                    Relaxation_min = max(Relaxation_min, abs(HSO4_mob_variation)/HSO4_mob_bef);
                                    Acceleration_max = min(Acceleration_max, HSO4_mob_bef/abs(HSO4_mob_variation));
                                else
                                    HSO4_mob_variation_ini = 0;
                                end
                            elseif HSO4_mob_variation > 0
                                Acceleration_max = min(Acceleration_max, Variation_rel_max*HSO4_mob_bef/abs(HSO4_mob_variation));
                            end

                            H_mob_variation_ini = H_mob_variation;
                            if H_mob_variation < 0
                                if H_mob_bef > 0
                                    Relaxation_min = max(Relaxation_min, abs(H_mob_variation)/H_mob_bef);
                                    Acceleration_max = min(Acceleration_max, H_mob_bef/abs(H_mob_variation));
                                else
                                    H_mob_variation_ini = 0;
                                end
                            elseif H_mob_variation > 0
                                Acceleration_max = min(Acceleration_max, Variation_rel_max*H_mob_bef/abs(H_mob_variation));                   
                            end

                            A1H_mob_variation_ini = A1H_mob_variation;
                            if A1H_mob_variation < 0
                                if A1H_mob_bef > 0
                                    Relaxation_min = max(Relaxation_min, abs(A1H_mob_variation)/A1H_mob_bef);
                                    Acceleration_max = min(Acceleration_max, A1H_mob_bef/abs(A1H_mob_variation));
                                else
                                    A1H_mob_variation_ini = 0;
                                end
                            elseif A1H_mob_variation > 0
                                Acceleration_max = min(Acceleration_max, Variation_rel_max*A1H_mob_bef/abs(A1H_mob_variation));
                            end

                            A2H_mob_variation_ini = A2H_mob_variation;
                            if A2H_mob_variation < 0
                                if A2H_mob_bef > 0
                                    Relaxation_min = max(Relaxation_min, abs(A2H_mob_variation)/A2H_mob_bef);
                                    Acceleration_max = min(Acceleration_max, A2H_mob_bef/abs(A2H_mob_variation));
                                else
                                    A2H_mob_variation_ini = 0;
                                end
                            elseif A2H_mob_variation > 0
                                Acceleration_max = min(Acceleration_max, Variation_rel_max*A2H_mob_bef/abs(A2H_mob_variation));
                            end                

                            A3H_mob_variation_ini = A3H_mob_variation;
                            if A3H_mob_variation < 0
                                if A3H_mob_bef > 0
                                    Relaxation_min = max(Relaxation_min, abs(A3H_mob_variation)/A3H_mob_bef);
                                    Acceleration_max = min(Acceleration_max, A3H_mob_bef/abs(A3H_mob_variation));
                                else
                                    A3H_mob_variation_ini = 0;
                                end
                            elseif A3H_mob_variation > 0
                                Acceleration_max = min(Acceleration_max, Variation_rel_max*A3H_mob_bef/abs(A3H_mob_variation));
                            end

                            N_iter_loc_2_new = N_iter_loc_2_new + 1;
                            N_iter_loc_2_tot_new  = N_iter_loc_2_tot_new + 1;
                            Relaxation_ini = ceil(Relaxation_min * Facteur_freinage) / Facteur_freinage;
                            Acceleration_max = floor(Acceleration_max);

                            % Nouvelle solution avec le facteur de relaxation minimum               
                            HSO4_stat_next = HSO4_stat_bef + HSO4_stat_variation_ini / Relaxation_ini;
                            HSO4_mob_next = HSO4_mob_bef + HSO4_mob_variation_ini / Relaxation_ini;
                            H_mob_next = H_mob_bef + H_mob_variation_ini / Relaxation_ini;
                            A1H_mob_next = A1H_mob_bef + A1H_mob_variation_ini / Relaxation_ini;
                            A2H_mob_next = A2H_mob_bef + A2H_mob_variation_ini / Relaxation_ini;
                            A3H_mob_next = A3H_mob_bef + A3H_mob_variation_ini / Relaxation_ini;

                            % Calcul du nouveau vecteur "bilans de matière" avec facteur de relaxation minimum
                            V1 = HSO4_mob_next*H_mob_next;
                            V2 = HSO4_stat_next^2;
                            V4 = V2/V1;
                            V6 = HSO4_mob_next/H_mob_next;

                            if Modele_new < 7
                                V7 = N_SO4_min + C6_1*A1H_mob_next + C6_2*A2H_mob_next + C6_3*A3H_mob_next;
                                V8 = N_HSO4_min + C7_1*A1H_mob_next + C7_2*A2H_mob_next + C7_3*A3H_mob_next;
                                V9 = N_mat_min + C38_1*A1H_mob_next + C38_2*A2H_mob_next + C38_3*A3H_mob_next;

                                if  Modele_new == 0 || Modele_new == 2 || Modele_new == 3 || Modele_new == 6;
                                    V3_1 = HSO4_stat_next*A1H_mob_next;
                                    V5_1 = V3_1/V1;
                                    V10_1 = A1H_mob_next/V7;
                                    V11_1 = A1H_mob_next/V9;
                                    V12_1 = V4*V10_1;
                                    V13_1 = V3_1/V8;
                                    V14_1 = A1H_mob_next/H_mob_next;
                                end

                                if  Modele_new == 0 || Modele_new == 1 || Modele_new == 3 || Modele_new == 5;
                                    V3_2 = HSO4_stat_next*A2H_mob_next;
                                    V5_2 = V3_2/V1;
                                    V10_2 = A2H_mob_next/V7;
                                    V11_2 = A2H_mob_next/V9;
                                    V12_2 = V4*V10_2;
                                    V13_2 = V3_2/V8;
                                    V14_2 = A2H_mob_next/H_mob_next;
                                end

                                if  Modele_new == 0 || Modele_new == 1 || Modele_new == 2 || Modele_new == 4;
                                    V3_3 = HSO4_stat_next*A3H_mob_next;
                                    V5_3 = V3_3/V1;
                                    V10_3 = A3H_mob_next/V7;
                                    V11_3 = A3H_mob_next/V9;
                                    V12_3 = V4*V10_3;
                                    V13_3 = V3_3/V8;
                                    V14_3 = A3H_mob_next/H_mob_next;
                                end
                            end

                            if Modele_new == 0                    
                                Bilan_Site_next = C9_corr*V4 + C10*HSO4_stat_next + C11_1_corr*V5_1 + C11_2_corr*V5_2 + C11_3_corr*V5_3 - 1;
                                Bilan_SO4_next = C20* HSO4_mob_next + C21*V6 + C22*HSO4_stat_next + C23*V4 - 1;
                                Bilan_H_next = C25*(H_mob_next+HSO4_mob_next+A1H_mob_next+A2H_mob_next+A3H_mob_next) + C26*HSO4_stat_next + C27_1*V12_1 + C27_2*V12_2 + C27_3*V12_3 + C28_1*V13_1 + C28_2*V13_2 + C28_3*V13_3 + C29_1*V11_1 + C29_2*V11_2 + C29_3*V11_3 - 1;             
                                Bilan_A1_next = C31_1*A1H_mob_next + C32_1*V14_1 + C33_1*V5_1 + C34_1*V12_1 + C35_1*V13_1 + C36_1*V11_1 - 1;               
                                Bilan_A2_next = C31_2*A2H_mob_next + C32_2*V14_2 + C33_2*V5_2 + C34_2*V12_2 + C35_2*V13_2 + C36_2*V11_2 - 1;
                                Bilan_A3_next = C31_3*A3H_mob_next + C32_3*V14_3 + C33_3*V5_3 + C34_3*V12_3 + C35_3*V13_3 + C36_3*V11_3 - 1;

                            elseif Modele_new == 1                    
                                Bilan_Site_next = C9_corr*V4 + C10*HSO4_stat_next + C11_2_corr*V5_2 + C11_3_corr*V5_3 - 1;
                                Bilan_SO4_next = C20* HSO4_mob_next + C21*V6 + C22*HSO4_stat_next + C23*V4 - 1;
                                Bilan_H_next = C25*(H_mob_next+HSO4_mob_next+A2H_mob_next+A3H_mob_next) + C26*HSO4_stat_next + C27_2*V12_2 + C27_3*V12_3 + C28_2*V13_2 + C28_3*V13_3 + C29_2*V11_2 + C29_3*V11_3 - 1;             
                                Bilan_A1_next = 0;               
                                Bilan_A2_next = C31_2*A2H_mob_next + C32_2*V14_2 + C33_2*V5_2 + C34_2*V12_2 + C35_2*V13_2 + C36_2*V11_2 - 1;
                                Bilan_A3_next = C31_3*A3H_mob_next + C32_3*V14_3 + C33_3*V5_3 + C34_3*V12_3 + C35_3*V13_3 + C36_3*V11_3 - 1;

                            elseif Modele_new == 2                    
                                Bilan_Site_next = C9_corr*V4 + C10*HSO4_stat_next + C11_1_corr*V5_1 + C11_3_corr*V5_3 - 1;
                                Bilan_SO4_next = C20* HSO4_mob_next + C21*V6 + C22*HSO4_stat_next + C23*V4 - 1;
                                Bilan_H_next = C25*(H_mob_next+HSO4_mob_next+A1H_mob_next+A3H_mob_next) + C26*HSO4_stat_next + C27_1*V12_1 + C27_3*V12_3 + C28_1*V13_1 + C28_3*V13_3 + C29_1*V11_1 + C29_3*V11_3 - 1;             
                                Bilan_A1_next = C31_1*A1H_mob_next + C32_1*V14_1 + C33_1*V5_1 + C34_1*V12_1 + C35_1*V13_1 + C36_1*V11_1 - 1;               
                                Bilan_A2_next = 0;
                                Bilan_A3_next = C31_3*A3H_mob_next + C32_3*V14_3 + C33_3*V5_3 + C34_3*V12_3 + C35_3*V13_3 + C36_3*V11_3 - 1;

                            elseif Modele_new == 3                    
                                Bilan_Site_next = C9_corr*V4 + C10*HSO4_stat_next + C11_1_corr*V5_1 + C11_2_corr*V5_2 - 1;
                                Bilan_SO4_next = C20* HSO4_mob_next + C21*V6 + C22*HSO4_stat_next + C23*V4 - 1;
                                Bilan_H_next = C25*(H_mob_next+HSO4_mob_next+A1H_mob_next+A2H_mob_next) + C26*HSO4_stat_next + C27_1*V12_1 + C27_2*V12_2 + C28_1*V13_1 + C28_2*V13_2 + C29_1*V11_1 + C29_2*V11_2 - 1;             
                                Bilan_A1_next = C31_1*A1H_mob_next + C32_1*V14_1 + C33_1*V5_1 + C34_1*V12_1 + C35_1*V13_1 + C36_1*V11_1 - 1;               
                                Bilan_A2_next = C31_2*A2H_mob_next + C32_2*V14_2 + C33_2*V5_2 + C34_2*V12_2 + C35_2*V13_2 + C36_2*V11_2 - 1;
                                Bilan_A3_next = 0;

                            elseif Modele_new == 4                    
                                Bilan_Site_next = C9_corr*V4 + C10*HSO4_stat_next + C11_3_corr*V5_3 - 1;
                                Bilan_SO4_next = C20* HSO4_mob_next + C21*V6 + C22*HSO4_stat_next + C23*V4 - 1;
                                Bilan_H_next = C25*(H_mob_next+HSO4_mob_next+A3H_mob_next) + C26*HSO4_stat_next + C27_3*V12_3 + C28_3*V13_3 + C29_3*V11_3 - 1;             
                                Bilan_A1_next = 0;               
                                Bilan_A2_next = 0;
                                Bilan_A3_next = C31_3*A3H_mob_next + C32_3*V14_3 + C33_3*V5_3 + C34_3*V12_3 + C35_3*V13_3 + C36_3*V11_3 - 1;

                            elseif Modele_new == 5                    
                                Bilan_Site_next = C9_corr*V4 + C10*HSO4_stat_next + C11_2_corr*V5_2 - 1;
                                Bilan_SO4_next = C20* HSO4_mob_next + C21*V6 + C22*HSO4_stat_next + C23*V4 - 1;
                                Bilan_H_next = C25*(H_mob_next+HSO4_mob_next+A2H_mob_next) + C26*HSO4_stat_next + C27_2*V12_2 + C28_2*V13_2 + C29_2*V11_2 - 1;             
                                Bilan_A1_next = 0;               
                                Bilan_A2_next = C31_2*A2H_mob_next + C32_2*V14_2 + C33_2*V5_2 + C34_2*V12_2 + C35_2*V13_2 + C36_2*V11_2 - 1;
                                Bilan_A3_next = 0;

                            elseif Modele_new == 6                    
                                Bilan_Site_next = C9_corr*V4 + C10*HSO4_stat_next + C11_1_corr*V5_1 - 1;
                                Bilan_SO4_next = C20* HSO4_mob_next + C21*V6 + C22*HSO4_stat_next + C23*V4 - 1;
                                Bilan_H_next = C25*(H_mob_next+HSO4_mob_next+A1H_mob_next) + C26*HSO4_stat_next + C27_1*V12_1 + C28_1*V13_1 + C29_1*V11_1 - 1;             
                                Bilan_A1_next = C31_1*A1H_mob_next + C32_1*V14_1 + C33_1*V5_1 + C34_1*V12_1 + C35_1*V13_1 + C36_1*V11_1 - 1;               
                                Bilan_A2_next = 0;
                                Bilan_A3_next = 0;

                            else
                                Bilan_Site_next = C9_corr*V4 + C10*HSO4_stat_next - 1;
                                Bilan_SO4_next = C20* HSO4_mob_next + C21*V6 + C22*HSO4_stat_next + C23*V4 - 1;
                                Bilan_H_next = C25*(H_mob_next+HSO4_mob_next) + C26*HSO4_stat_next - 1;
                                Bilan_A1_next = 0;               
                                Bilan_A2_next = 0;
                                Bilan_A3_next = 0;
                            end

                            Bilans_next = [Bilan_Site_next; Bilan_SO4_next; Bilan_H_next; Bilan_A1_next; Bilan_A2_next; Bilan_A3_next];
                            N_Bilans_next = norm(Bilans_next);
                            Residu_rel_max_new = max(abs(Bilans_next));
                            Residu_abs_max_new = max(abs(Bilans_next.*[qmax_ei; C_SO4; C_H; C_A1H; C_A2H; C_A3H]));

                            % Relaxation de la méthode tant que la nouvelle solution est moins bonne 
                            if N_Bilans_next >= N_Bilans_bef
                                N_iter = 1;
                                Variation_negligeable = 0;
                                while N_Bilans_next >= N_Bilans_bef && ((Residu_abs_max_new >= Precision_Bilan_abs && Residu_rel_max_new >= Precision_Bilan_rel) || Variation_negligeable == 0) && N_iter <= N_iter_max_F
                                    N_iter = N_iter + 1;
                                    N_iter_loc_2_new = N_iter_loc_2_new + 1;
                                    N_iter_loc_2_tot_new  = N_iter_loc_2_tot_new + 1;
                                    Relaxation = Relaxation_ini * Facteur_freinage ^ (N_iter - 1);

                                    % Calcul nouvelle solution relaxée
                                    HSO4_stat_next = HSO4_stat_bef + HSO4_stat_variation_ini / Relaxation;
                                    HSO4_mob_next = HSO4_mob_bef + HSO4_mob_variation_ini / Relaxation;
                                    H_mob_next = H_mob_bef + H_mob_variation_ini / Relaxation;
                                    A1H_mob_next = A1H_mob_bef + A1H_mob_variation_ini / Relaxation;
                                    A2H_mob_next = A2H_mob_bef + A2H_mob_variation_ini / Relaxation;
                                    A3H_mob_next = A3H_mob_bef + A3H_mob_variation_ini / Relaxation;

                                    Variation_negligeable = max(abs([HSO4_stat_variation_ini,HSO4_mob_variation_ini,H_mob_variation_ini,A1H_mob_variation_ini,A2H_mob_variation_ini,A3H_mob_variation_ini]/ Relaxation)) < Precision_Conc_abs;

                                    % Calcul du nouveau vecteur "bilans de matière"
                                    V1 = HSO4_mob_next*H_mob_next;
                                    V2 = HSO4_stat_next^2;
                                    V4 = V2/V1;
                                    V6 = HSO4_mob_next/H_mob_next;

                                    if Modele_new < 7
                                        V7 = N_SO4_min + C6_1*A1H_mob_next + C6_2*A2H_mob_next + C6_3*A3H_mob_next;
                                        V8 = N_HSO4_min + C7_1*A1H_mob_next + C7_2*A2H_mob_next + C7_3*A3H_mob_next;
                                        V9 = N_mat_min + C38_1*A1H_mob_next + C38_2*A2H_mob_next + C38_3*A3H_mob_next;

                                        if  Modele_new == 0 || Modele_new == 2 || Modele_new == 3 || Modele_new == 6;
                                            V3_1 = HSO4_stat_next*A1H_mob_next;
                                            V5_1 = V3_1/V1;
                                            V10_1 = A1H_mob_next/V7;
                                            V11_1 = A1H_mob_next/V9;
                                            V12_1 = V4*V10_1;
                                            V13_1 = V3_1/V8;
                                            V14_1 = A1H_mob_next/H_mob_next;
                                        end

                                        if  Modele_new == 0 || Modele_new == 1 || Modele_new == 3 || Modele_new == 5;
                                            V3_2 = HSO4_stat_next*A2H_mob_next;
                                            V5_2 = V3_2/V1;
                                            V10_2 = A2H_mob_next/V7;
                                            V11_2 = A2H_mob_next/V9;
                                            V12_2 = V4*V10_2;
                                            V13_2 = V3_2/V8;
                                            V14_2 = A2H_mob_next/H_mob_next;
                                        end

                                        if  Modele_new == 0 || Modele_new == 1 || Modele_new == 2 || Modele_new == 4;
                                            V3_3 = HSO4_stat_next*A3H_mob_next;
                                            V5_3 = V3_3/V1;
                                            V10_3 = A3H_mob_next/V7;
                                            V11_3 = A3H_mob_next/V9;
                                            V12_3 = V4*V10_3;
                                            V13_3 = V3_3/V8;
                                            V14_3 = A3H_mob_next/H_mob_next;
                                        end
                                    end

                                    if Modele_new == 0                    
                                        Bilan_Site_next = C9_corr*V4 + C10*HSO4_stat_next + C11_1_corr*V5_1 + C11_2_corr*V5_2 + C11_3_corr*V5_3 - 1;
                                        Bilan_SO4_next = C20* HSO4_mob_next + C21*V6 + C22*HSO4_stat_next + C23*V4 - 1;
                                        Bilan_H_next = C25*(H_mob_next+HSO4_mob_next+A1H_mob_next+A2H_mob_next+A3H_mob_next) + C26*HSO4_stat_next + C27_1*V12_1 + C27_2*V12_2 + C27_3*V12_3 + C28_1*V13_1 + C28_2*V13_2 + C28_3*V13_3 + C29_1*V11_1 + C29_2*V11_2 + C29_3*V11_3 - 1;             
                                        Bilan_A1_next = C31_1*A1H_mob_next + C32_1*V14_1 + C33_1*V5_1 + C34_1*V12_1 + C35_1*V13_1 + C36_1*V11_1 - 1;               
                                        Bilan_A2_next = C31_2*A2H_mob_next + C32_2*V14_2 + C33_2*V5_2 + C34_2*V12_2 + C35_2*V13_2 + C36_2*V11_2 - 1;
                                        Bilan_A3_next = C31_3*A3H_mob_next + C32_3*V14_3 + C33_3*V5_3 + C34_3*V12_3 + C35_3*V13_3 + C36_3*V11_3 - 1;

                                    elseif Modele_new == 1                    
                                        Bilan_Site_next = C9_corr*V4 + C10*HSO4_stat_next + C11_2_corr*V5_2 + C11_3_corr*V5_3 - 1;
                                        Bilan_SO4_next = C20* HSO4_mob_next + C21*V6 + C22*HSO4_stat_next + C23*V4 - 1;
                                        Bilan_H_next = C25*(H_mob_next+HSO4_mob_next+A2H_mob_next+A3H_mob_next) + C26*HSO4_stat_next + C27_2*V12_2 + C27_3*V12_3 + C28_2*V13_2 + C28_3*V13_3 + C29_2*V11_2 + C29_3*V11_3 - 1;             
                                        Bilan_A1_next = 0;               
                                        Bilan_A2_next = C31_2*A2H_mob_next + C32_2*V14_2 + C33_2*V5_2 + C34_2*V12_2 + C35_2*V13_2 + C36_2*V11_2 - 1;
                                        Bilan_A3_next = C31_3*A3H_mob_next + C32_3*V14_3 + C33_3*V5_3 + C34_3*V12_3 + C35_3*V13_3 + C36_3*V11_3 - 1;

                                    elseif Modele_new == 2                    
                                        Bilan_Site_next = C9_corr*V4 + C10*HSO4_stat_next + C11_1_corr*V5_1 + C11_3_corr*V5_3 - 1;
                                        Bilan_SO4_next = C20* HSO4_mob_next + C21*V6 + C22*HSO4_stat_next + C23*V4 - 1;
                                        Bilan_H_next = C25*(H_mob_next+HSO4_mob_next+A1H_mob_next+A3H_mob_next) + C26*HSO4_stat_next + C27_1*V12_1 + C27_3*V12_3 + C28_1*V13_1 + C28_3*V13_3 + C29_1*V11_1 + C29_3*V11_3 - 1;             
                                        Bilan_A1_next = C31_1*A1H_mob_next + C32_1*V14_1 + C33_1*V5_1 + C34_1*V12_1 + C35_1*V13_1 + C36_1*V11_1 - 1;               
                                        Bilan_A2_next = 0;
                                        Bilan_A3_next = C31_3*A3H_mob_next + C32_3*V14_3 + C33_3*V5_3 + C34_3*V12_3 + C35_3*V13_3 + C36_3*V11_3 - 1;

                                    elseif Modele_new == 3                    
                                        Bilan_Site_next = C9_corr*V4 + C10*HSO4_stat_next + C11_1_corr*V5_1 + C11_2_corr*V5_2 - 1;
                                        Bilan_SO4_next = C20* HSO4_mob_next + C21*V6 + C22*HSO4_stat_next + C23*V4 - 1;
                                        Bilan_H_next = C25*(H_mob_next+HSO4_mob_next+A1H_mob_next+A2H_mob_next) + C26*HSO4_stat_next + C27_1*V12_1 + C27_2*V12_2 + C28_1*V13_1 + C28_2*V13_2 + C29_1*V11_1 + C29_2*V11_2 - 1;             
                                        Bilan_A1_next = C31_1*A1H_mob_next + C32_1*V14_1 + C33_1*V5_1 + C34_1*V12_1 + C35_1*V13_1 + C36_1*V11_1 - 1;               
                                        Bilan_A2_next = C31_2*A2H_mob_next + C32_2*V14_2 + C33_2*V5_2 + C34_2*V12_2 + C35_2*V13_2 + C36_2*V11_2 - 1;
                                        Bilan_A3_next = 0;

                                    elseif Modele_new == 4                    
                                        Bilan_Site_next = C9_corr*V4 + C10*HSO4_stat_next + C11_3_corr*V5_3 - 1;
                                        Bilan_SO4_next = C20* HSO4_mob_next + C21*V6 + C22*HSO4_stat_next + C23*V4 - 1;
                                        Bilan_H_next = C25*(H_mob_next+HSO4_mob_next+A3H_mob_next) + C26*HSO4_stat_next + C27_3*V12_3 + C28_3*V13_3 + C29_3*V11_3 - 1;             
                                        Bilan_A1_next = 0;               
                                        Bilan_A2_next = 0;
                                        Bilan_A3_next = C31_3*A3H_mob_next + C32_3*V14_3 + C33_3*V5_3 + C34_3*V12_3 + C35_3*V13_3 + C36_3*V11_3 - 1;

                                    elseif Modele_new == 5                    
                                        Bilan_Site_next = C9_corr*V4 + C10*HSO4_stat_next + C11_2_corr*V5_2 - 1;
                                        Bilan_SO4_next = C20* HSO4_mob_next + C21*V6 + C22*HSO4_stat_next + C23*V4 - 1;
                                        Bilan_H_next = C25*(H_mob_next+HSO4_mob_next+A2H_mob_next) + C26*HSO4_stat_next + C27_2*V12_2 + C28_2*V13_2 + C29_2*V11_2 - 1;             
                                        Bilan_A1_next = 0;               
                                        Bilan_A2_next = C31_2*A2H_mob_next + C32_2*V14_2 + C33_2*V5_2 + C34_2*V12_2 + C35_2*V13_2 + C36_2*V11_2 - 1;
                                        Bilan_A3_next = 0;

                                    elseif Modele_new == 6                    
                                        Bilan_Site_next = C9_corr*V4 + C10*HSO4_stat_next + C11_1_corr*V5_1 - 1;
                                        Bilan_SO4_next = C20* HSO4_mob_next + C21*V6 + C22*HSO4_stat_next + C23*V4 - 1;
                                        Bilan_H_next = C25*(H_mob_next+HSO4_mob_next+A1H_mob_next) + C26*HSO4_stat_next + C27_1*V12_1 + C28_1*V13_1 + C29_1*V11_1 - 1;             
                                        Bilan_A1_next = C31_1*A1H_mob_next + C32_1*V14_1 + C33_1*V5_1 + C34_1*V12_1 + C35_1*V13_1 + C36_1*V11_1 - 1;               
                                        Bilan_A2_next = 0;
                                        Bilan_A3_next = 0;

                                    else
                                        Bilan_Site_next = C9_corr*V4 + C10*HSO4_stat_next - 1;
                                        Bilan_SO4_next = C20* HSO4_mob_next + C21*V6 + C22*HSO4_stat_next + C23*V4 - 1;
                                        Bilan_H_next = C25*(H_mob_next+HSO4_mob_next) + C26*HSO4_stat_next - 1;
                                        Bilan_A1_next = 0;               
                                        Bilan_A2_next = 0;
                                        Bilan_A3_next = 0;
                                    end

                                    Bilans_next = [Bilan_Site_next; Bilan_SO4_next; Bilan_H_next; Bilan_A1_next; Bilan_A2_next; Bilan_A3_next];
                                    N_Bilans_next = norm(Bilans_next);
                                    Residu_rel_max_new = max(abs(Bilans_next));
                                    Residu_abs_max_new = max(abs(Bilans_next.*[qmax_ei; C_SO4; C_H; C_A1H; C_A2H; C_A3H]));
                                end

                            % Accélération si possible de la méthode tant que la nouvelle solution est meilleure que la précédente
                            elseif Acceleration_max >= Facteur_acceleration
                                N_iter = 1;
                                Acceleration = Facteur_acceleration;
                                N_Bilans_next_next = N_Bilans_next;
                                while N_Bilans_next_next < N_Bilans_next && Acceleration <= Acceleration_max && N_iter <= N_iter_max_F
                                    N_iter = N_iter + 1;
                                    N_iter_loc_2_new = N_iter_loc_2_new + 1;
                                    N_iter_loc_2_tot_new  = N_iter_loc_2_tot_new + 1;

                                    % Calcul nouvelle solution accélérée
                                    HSO4_stat_next_next = HSO4_stat_bef + HSO4_stat_variation_ini * Acceleration;
                                    HSO4_mob_next_next = HSO4_mob_bef + HSO4_mob_variation_ini * Acceleration;
                                    H_mob_next_next = H_mob_bef + H_mob_variation_ini * Acceleration;
                                    A1H_mob_next_next = A1H_mob_bef + A1H_mob_variation_ini * Acceleration;
                                    A2H_mob_next_next = A2H_mob_bef + A2H_mob_variation_ini * Acceleration;
                                    A3H_mob_next_next = A3H_mob_bef + A3H_mob_variation_ini * Acceleration;

                                    % Calcul du nouveau vecteur "bilans de matière"
                                    V1 = HSO4_mob_next_next*H_mob_next_next;
                                    V2 = HSO4_stat_next_next^2;
                                    V4 = V2/V1;
                                    V6 = HSO4_mob_next_next/H_mob_next_next;

                                    if Modele_new < 7
                                        V7 = N_SO4_min + C6_1*A1H_mob_next_next + C6_2*A2H_mob_next_next + C6_3*A3H_mob_next_next;
                                        V8 = N_HSO4_min + C7_1*A1H_mob_next_next + C7_2*A2H_mob_next_next + C7_3*A3H_mob_next_next;
                                        V9 = N_mat_min + C38_1*A1H_mob_next_next + C38_2*A2H_mob_next_next + C38_3*A3H_mob_next_next;

                                        if  Modele_new == 0 || Modele_new == 2 || Modele_new == 3 || Modele_new == 6;
                                            V3_1 = HSO4_stat_next_next*A1H_mob_next_next;
                                            V5_1 = V3_1/V1;
                                            V10_1 = A1H_mob_next_next/V7;
                                            V11_1 = A1H_mob_next_next/V9;
                                            V12_1 = V4*V10_1;
                                            V13_1 = V3_1/V8;
                                            V14_1 = A1H_mob_next_next/H_mob_next_next;
                                        end

                                        if  Modele_new == 0 || Modele_new == 1 || Modele_new == 3 || Modele_new == 5;
                                            V3_2 = HSO4_stat_next_next*A2H_mob_next_next;
                                            V5_2 = V3_2/V1;
                                            V10_2 = A2H_mob_next_next/V7;
                                            V11_2 = A2H_mob_next_next/V9;
                                            V12_2 = V4*V10_2;
                                            V13_2 = V3_2/V8;
                                            V14_2 = A2H_mob_next_next/H_mob_next_next;
                                        end

                                        if  Modele_new == 0 || Modele_new == 1 || Modele_new == 2 || Modele_new == 4;
                                            V3_3 = HSO4_stat_next_next*A3H_mob_next_next;
                                            V5_3 = V3_3/V1;
                                            V10_3 = A3H_mob_next_next/V7;
                                            V11_3 = A3H_mob_next_next/V9;
                                            V12_3 = V4*V10_3;
                                            V13_3 = V3_3/V8;
                                            V14_3 = A3H_mob_next_next/H_mob_next_next;
                                        end
                                    end

                                    if Modele_new == 0                    
                                        Bilan_Site_next_next = C9_corr*V4 + C10*HSO4_stat_next_next + C11_1_corr*V5_1 + C11_2_corr*V5_2 + C11_3_corr*V5_3 - 1;
                                        Bilan_SO4_next_next = C20* HSO4_mob_next_next + C21*V6 + C22*HSO4_stat_next_next + C23*V4 - 1;
                                        Bilan_H_next_next = C25*(H_mob_next_next+HSO4_mob_next_next+A1H_mob_next_next+A2H_mob_next_next+A3H_mob_next_next) + C26*HSO4_stat_next_next + C27_1*V12_1 + C27_2*V12_2 + C27_3*V12_3 + C28_1*V13_1 + C28_2*V13_2 + C28_3*V13_3 + C29_1*V11_1 + C29_2*V11_2 + C29_3*V11_3 - 1;             
                                        Bilan_A1_next_next = C31_1*A1H_mob_next_next + C32_1*V14_1 + C33_1*V5_1 + C34_1*V12_1 + C35_1*V13_1 + C36_1*V11_1 - 1;               
                                        Bilan_A2_next_next = C31_2*A2H_mob_next_next + C32_2*V14_2 + C33_2*V5_2 + C34_2*V12_2 + C35_2*V13_2 + C36_2*V11_2 - 1;
                                        Bilan_A3_next_next = C31_3*A3H_mob_next_next + C32_3*V14_3 + C33_3*V5_3 + C34_3*V12_3 + C35_3*V13_3 + C36_3*V11_3 - 1;

                                    elseif Modele_new == 1                    
                                        Bilan_Site_next_next = C9_corr*V4 + C10*HSO4_stat_next_next + C11_2_corr*V5_2 + C11_3_corr*V5_3 - 1;
                                        Bilan_SO4_next_next = C20* HSO4_mob_next_next + C21*V6 + C22*HSO4_stat_next_next + C23*V4 - 1;
                                        Bilan_H_next_next = C25*(H_mob_next_next+HSO4_mob_next_next+A2H_mob_next_next+A3H_mob_next_next) + C26*HSO4_stat_next_next + C27_2*V12_2 + C27_3*V12_3 + C28_2*V13_2 + C28_3*V13_3 + C29_2*V11_2 + C29_3*V11_3 - 1;             
                                        Bilan_A1_next_next = 0;               
                                        Bilan_A2_next_next = C31_2*A2H_mob_next_next + C32_2*V14_2 + C33_2*V5_2 + C34_2*V12_2 + C35_2*V13_2 + C36_2*V11_2 - 1;
                                        Bilan_A3_next_next = C31_3*A3H_mob_next_next + C32_3*V14_3 + C33_3*V5_3 + C34_3*V12_3 + C35_3*V13_3 + C36_3*V11_3 - 1;

                                    elseif Modele_new == 2                    
                                        Bilan_Site_next_next = C9_corr*V4 + C10*HSO4_stat_next_next + C11_1_corr*V5_1 + C11_3_corr*V5_3 - 1;
                                        Bilan_SO4_next_next = C20* HSO4_mob_next_next + C21*V6 + C22*HSO4_stat_next_next + C23*V4 - 1;
                                        Bilan_H_next_next = C25*(H_mob_next_next+HSO4_mob_next_next+A1H_mob_next_next+A3H_mob_next_next) + C26*HSO4_stat_next_next + C27_1*V12_1 + C27_3*V12_3 + C28_1*V13_1 + C28_3*V13_3 + C29_1*V11_1 + C29_3*V11_3 - 1;             
                                        Bilan_A1_next_next = C31_1*A1H_mob_next_next + C32_1*V14_1 + C33_1*V5_1 + C34_1*V12_1 + C35_1*V13_1 + C36_1*V11_1 - 1;               
                                        Bilan_A2_next_next = 0;
                                        Bilan_A3_next_next = C31_3*A3H_mob_next_next + C32_3*V14_3 + C33_3*V5_3 + C34_3*V12_3 + C35_3*V13_3 + C36_3*V11_3 - 1;

                                    elseif Modele_new == 3                    
                                        Bilan_Site_next_next = C9_corr*V4 + C10*HSO4_stat_next_next + C11_1_corr*V5_1 + C11_2_corr*V5_2 - 1;
                                        Bilan_SO4_next_next = C20* HSO4_mob_next_next + C21*V6 + C22*HSO4_stat_next_next + C23*V4 - 1;
                                        Bilan_H_next_next = C25*(H_mob_next_next+HSO4_mob_next_next+A1H_mob_next_next+A2H_mob_next_next) + C26*HSO4_stat_next_next + C27_1*V12_1 + C27_2*V12_2 + C28_1*V13_1 + C28_2*V13_2 + C29_1*V11_1 + C29_2*V11_2 - 1;             
                                        Bilan_A1_next_next = C31_1*A1H_mob_next_next + C32_1*V14_1 + C33_1*V5_1 + C34_1*V12_1 + C35_1*V13_1 + C36_1*V11_1 - 1;               
                                        Bilan_A2_next_next = C31_2*A2H_mob_next_next + C32_2*V14_2 + C33_2*V5_2 + C34_2*V12_2 + C35_2*V13_2 + C36_2*V11_2 - 1;
                                        Bilan_A3_next_next = 0;

                                    elseif Modele_new == 4                    
                                        Bilan_Site_next_next = C9_corr*V4 + C10*HSO4_stat_next_next + C11_3_corr*V5_3 - 1;
                                        Bilan_SO4_next_next = C20* HSO4_mob_next_next + C21*V6 + C22*HSO4_stat_next_next + C23*V4 - 1;
                                        Bilan_H_next_next = C25*(H_mob_next_next+HSO4_mob_next_next+A3H_mob_next_next) + C26*HSO4_stat_next_next + C27_3*V12_3 + C28_3*V13_3 + C29_3*V11_3 - 1;             
                                        Bilan_A1_next_next = 0;               
                                        Bilan_A2_next_next = 0;
                                        Bilan_A3_next_next = C31_3*A3H_mob_next_next + C32_3*V14_3 + C33_3*V5_3 + C34_3*V12_3 + C35_3*V13_3 + C36_3*V11_3 - 1;

                                    elseif Modele_new == 5                    
                                        Bilan_Site_next_next = C9_corr*V4 + C10*HSO4_stat_next_next + C11_2_corr*V5_2 - 1;
                                        Bilan_SO4_next_next = C20* HSO4_mob_next_next + C21*V6 + C22*HSO4_stat_next_next + C23*V4 - 1;
                                        Bilan_H_next_next = C25*(H_mob_next_next+HSO4_mob_next_next+A2H_mob_next_next) + C26*HSO4_stat_next_next + C27_2*V12_2 + C28_2*V13_2 + C29_2*V11_2 - 1;             
                                        Bilan_A1_next_next = 0;               
                                        Bilan_A2_next_next = C31_2*A2H_mob_next_next + C32_2*V14_2 + C33_2*V5_2 + C34_2*V12_2 + C35_2*V13_2 + C36_2*V11_2 - 1;
                                        Bilan_A3_next_next = 0;

                                    elseif Modele_new == 6                    
                                        Bilan_Site_next_next = C9_corr*V4 + C10*HSO4_stat_next_next + C11_1_corr*V5_1 - 1;
                                        Bilan_SO4_next_next = C20* HSO4_mob_next_next + C21*V6 + C22*HSO4_stat_next_next + C23*V4 - 1;
                                        Bilan_H_next_next = C25*(H_mob_next_next+HSO4_mob_next_next+A1H_mob_next_next) + C26*HSO4_stat_next_next + C27_1*V12_1 + C28_1*V13_1 + C29_1*V11_1 - 1;             
                                        Bilan_A1_next_next = C31_1*A1H_mob_next_next + C32_1*V14_1 + C33_1*V5_1 + C34_1*V12_1 + C35_1*V13_1 + C36_1*V11_1 - 1;               
                                        Bilan_A2_next_next = 0;
                                        Bilan_A3_next_next = 0;

                                    else
                                        Bilan_Site_next_next = C9_corr*V4 + C10*HSO4_stat_next_next - 1;
                                        Bilan_SO4_next_next = C20* HSO4_mob_next_next + C21*V6 + C22*HSO4_stat_next_next + C23*V4 - 1;
                                        Bilan_H_next_next = C25*(H_mob_next_next+HSO4_mob_next_next) + C26*HSO4_stat_next_next - 1;
                                        Bilan_A1_next_next = 0;               
                                        Bilan_A2_next_next = 0;
                                        Bilan_A3_next_next = 0;
                                    end

                                    Bilans_next_next = [Bilan_Site_next_next; Bilan_SO4_next_next; Bilan_H_next_next; Bilan_A1_next_next; Bilan_A2_next_next; Bilan_A3_next_next];
                                    N_Bilans_next_next = norm(Bilans_next_next);                     

                                    if N_Bilans_next_next < N_Bilans_next
                                        Acceleration = Facteur_acceleration ^ N_iter; 
                                        HSO4_stat_next = HSO4_stat_next_next;
                                        H_mob_next = H_mob_next_next;
                                        HSO4_mob_next = HSO4_mob_next_next;
                                        A1H_mob_next = A1H_mob_next_next;
                                        A2H_mob_next = A2H_mob_next_next;
                                        A3H_mob_next = A3H_mob_next_next;
                                        Bilans_next = Bilans_next_next;                           
                                        N_Bilans_next = N_Bilans_next_next;
                                        Residu_rel_max_new = max(abs(Bilans_next));
                                        Residu_abs_max_new = max(abs(Bilans_next.*[qmax_ei; C_SO4; C_H; C_A1H; C_A2H; C_A3H]));                           
                                    end
                                end

                                if N_Bilans_next_next > N_Bilans_next
                                    V1 = HSO4_mob_next*H_mob_next;
                                    V2 = HSO4_stat_next^2;
                                    V4 = V2/V1;
                                    V6 = HSO4_mob_next/H_mob_next;

                                    if Modele_new < 7
                                        V7 = N_SO4_min + C6_1*A1H_mob_next + C6_2*A2H_mob_next + C6_3*A3H_mob_next;
                                        V8 = N_HSO4_min + C7_1*A1H_mob_next + C7_2*A2H_mob_next + C7_3*A3H_mob_next;
                                        V9 = N_mat_min + C38_1*A1H_mob_next + C38_2*A2H_mob_next + C38_3*A3H_mob_next;

                                        if  Modele_new == 0 || Modele_new == 2 || Modele_new == 3 || Modele_new == 6;
                                            V3_1 = HSO4_stat_next*A1H_mob_next;
                                            V5_1 = V3_1/V1;
                                            V10_1 = A1H_mob_next/V7;
                                            V11_1 = A1H_mob_next/V9;
                                            V12_1 = V4*V10_1;
                                            V13_1 = V3_1/V8;
                                            V14_1 = A1H_mob_next/H_mob_next;
                                        end

                                        if  Modele_new == 0 || Modele_new == 1 || Modele_new == 3 || Modele_new == 5;
                                            V3_2 = HSO4_stat_next*A2H_mob_next;
                                            V5_2 = V3_2/V1;
                                            V10_2 = A2H_mob_next/V7;
                                            V11_2 = A2H_mob_next/V9;
                                            V12_2 = V4*V10_2;
                                            V13_2 = V3_2/V8;
                                            V14_2 = A2H_mob_next/H_mob_next;
                                        end

                                        if  Modele_new == 0 || Modele_new == 1 || Modele_new == 2 || Modele_new == 4;
                                            V3_3 = HSO4_stat_next*A3H_mob_next;
                                            V5_3 = V3_3/V1;
                                            V10_3 = A3H_mob_next/V7;
                                            V11_3 = A3H_mob_next/V9;
                                            V12_3 = V4*V10_3;
                                            V13_3 = V3_3/V8;
                                            V14_3 = A3H_mob_next/H_mob_next;
                                        end
                                    end
                                end
                            end

                        % Calcul des critères de convergence de la méthode numérique                
                            % Variation absolue maximum des paramètres
                            Ecart_abs_max_new = max(abs([HSO4_stat_next-HSO4_stat_bef,HSO4_mob_next-HSO4_mob_bef,H_mob_next-H_mob_bef,A1H_mob_next-A1H_mob_bef,A2H_mob_next-A2H_mob_bef,A3H_mob_next-A3H_mob_bef]));

                            % Variation relative maximum des paramètres
                            if HSO4_stat_bef > 0
                                HSO4_stat_variation_rel = (HSO4_stat_next-HSO4_stat_bef) / HSO4_stat_bef;
                            elseif HSO4_stat_next > 0
                                HSO4_stat_variation_rel = (HSO4_stat_next-HSO4_stat_bef) / HSO4_stat_next;
                            else
                                HSO4_stat_variation_rel = 0;
                            end
                            if HSO4_mob_bef > 0
                                HSO4_mob_variation_rel = (HSO4_mob_next-HSO4_mob_bef) / HSO4_mob_bef;
                            elseif HSO4_mob_next > 0
                                HSO4_mob_variation_rel = (HSO4_mob_next-HSO4_mob_bef) / HSO4_mob_next;
                            else
                                HSO4_mob_variation_rel = 0;
                            end
                            if H_mob_bef > 0
                                H_mob_variation_rel = (H_mob_next-H_mob_bef) / H_mob_bef;
                            elseif H_mob_next > 0
                                H_mob_variation_rel = (H_mob_next-H_mob_bef) / H_mob_next;
                            else
                                H_mob_variation_rel = 0;
                            end
                            if A1H_mob_bef > 0
                                A1H_mob_variation_rel = (A1H_mob_next-A1H_mob_bef) / A1H_mob_bef;
                            elseif A1H_mob_next > 0
                                A1H_mob_variation_rel = (A1H_mob_next-A1H_mob_bef) / A1H_mob_next;
                            else
                                A1H_mob_variation_rel = 0;
                            end
                            if A2H_mob_bef > 0
                                A2H_mob_variation_rel = (A2H_mob_next-A2H_mob_bef) / A2H_mob_bef;
                            elseif A2H_mob_next > 0
                                A2H_mob_variation_rel = (A2H_mob_next-A2H_mob_bef) / A2H_mob_next;
                            else
                                A2H_mob_variation_rel = 0;
                            end                    
                            if A3H_mob_bef > 0
                                A3H_mob_variation_rel = (A3H_mob_next-A3H_mob_bef) / A3H_mob_bef;
                            elseif A3H_mob_next > 0
                                A3H_mob_variation_rel = (A3H_mob_next-A3H_mob_bef) / A3H_mob_next;
                            else
                                A3H_mob_variation_rel = 0;
                            end               
                            Ecart_rel_max_new = max(abs([HSO4_stat_variation_rel; HSO4_mob_variation_rel; H_mob_variation_rel; A1H_mob_variation_rel; A2H_mob_variation_rel; A3H_mob_variation_rel]));

                        % Contrôle des résultats (mode débug)
                        Resultats = [HSO4_stat_next HSO4_mob_next H_mob_next A1H_mob_next A2H_mob_next A3H_mob_next;
                            HSO4_stat_bef HSO4_mob_bef H_mob_bef A1H_mob_bef A2H_mob_bef A3H_mob_bef;
                            HSO4_stat_old(n) HSO4_mob_old(n) H_mob_old(n) A1H_mob_old(n) A2H_mob_old(n) A3H_mob_old(n);
                            HSO4_stat_next-HSO4_stat_bef HSO4_mob_next-HSO4_mob_bef H_mob_next-H_mob_bef A1H_mob_next-A1H_mob_bef A2H_mob_next-A2H_mob_bef A3H_mob_next-A3H_mob_bef;
                            HSO4_stat_variation_rel HSO4_mob_variation_rel H_mob_variation_rel A1H_mob_variation_rel A2H_mob_variation_rel A3H_mob_variation_rel;
                            Bilans_next'.*[qmax_ei; C_SO4; C_H; C_A1H; C_A2H; C_A3H]';
                            Bilans_bef'.*[qmax_ei; C_SO4; C_H; C_A1H; C_A2H; C_A3H]';
                            Bilans_next';
                            Bilans_bef'];

                        HSO4_stat_bef = HSO4_stat_next;
                        H_mob_bef = H_mob_next;
                        HSO4_mob_bef = HSO4_mob_next;
                        A1H_mob_bef = A1H_mob_next;
                        A2H_mob_bef = A2H_mob_next;
                        A3H_mob_bef = A3H_mob_next;

                        Bilans_bef = Bilans_next;              
                        N_Bilans_bef = N_Bilans_next;
                    end

                    % Contrôle de la convergence de la méthode de Newton-Raphson
                    if (Ecart_abs_max_new < Precision_Conc_abs || Ecart_rel_max_new < Precision_Conc_rel) && (Residu_abs_max_new < Precision_Bilan_abs || Residu_rel_max_new < Precision_Bilan_rel)
                        Convergence_new = 1;
                    end

                    % Enregistrement de la nouvelle composition de l'étage théorique n
                    HSO4_stat_new(n) = HSO4_stat_next;
                    H_mob_new(n) = H_mob_next;
                    HSO4_mob_new(n) = HSO4_mob_next;
                    A1H_mob_new(n) = A1H_mob_next;
                    A2H_mob_new(n) = A2H_mob_next;
                    A3H_mob_new(n) = A3H_mob_next;                   
                    
                    % Calcul A- en solution
                    A1_mob_new(n) = A1H_mob_new(n) * Ka_A1H_0 * Activite_AH_new(n) / H_mob_new(n);
                    A2_mob_new(n) = A2H_mob_new(n) * Ka_A2H_0 * Activite_AH_new(n) / H_mob_new(n);
                    A3_mob_new(n) = A3H_mob_new(n) * Ka_A3H_0 * Activite_AH_new(n) / H_mob_new(n);

                    % Calcul Atot en solution
                    A1tot_mob_new(n) = A1H_mob_new(n) + A1_mob_new(n);
                    A2tot_mob_new(n) = A2H_mob_new(n) + A2_mob_new(n);
                    A3tot_mob_new(n) = A3H_mob_new(n) + A3_mob_new(n);

                    % Calcul SO4 2- en solution
                    SO4_mob_new(n) = HSO4_mob_new(n) * Ka_HSO4_0 * Activite_HSO4_new(n) / H_mob_new(n);

                    % Calcul SO4tot en solution
                    SO4tot_mob_new(n) = HSO4_mob_new(n) + SO4_mob_new(n);        

                    % Calcul SO4 2- sur résine
                    SO4_stat_new(n) = SO4_mob_new(n) * HSO4_stat_new(n)^2 / (HSO4_mob_new(n)^2 * Keq_HSO4_SO4);

                    % Calcul SO4tot sur résine
                    SO4tot_stat_new(n) = HSO4_stat_new(n) + SO4_stat_new(n);

                    % Calcul AH-HSO4 sur résine
                    A1H_HSO4_stat_new(n) = Ks_A1H_HSO4 * A1H_mob_new(n) * HSO4_stat_new(n) / V8;
                    A2H_HSO4_stat_new(n) = Ks_A2H_HSO4 * A2H_mob_new(n) * HSO4_stat_new(n) / V8;
                    A3H_HSO4_stat_new(n) = Ks_A3H_HSO4 * A3H_mob_new(n) * HSO4_stat_new(n) / V8;

                    % Calcul AH-SO4 sur résine
                    A1H_SO4_stat_new(n) = Ks_A1H_SO4 * A1H_mob_new(n) * SO4_stat_new(n) / V7;
                    A2H_SO4_stat_new(n) = Ks_A2H_SO4 * A2H_mob_new(n) * SO4_stat_new(n) / V7;
                    A3H_SO4_stat_new(n) = Ks_A3H_SO4 * A3H_mob_new(n) * SO4_stat_new(n) / V7;

                    % Calcul AH-mat sur résine
                    A1H_mat_stat_new(n) = C8_1 * A1H_mob_new(n) / V9;
                    A2H_mat_stat_new(n) = C8_2 * A2H_mob_new(n) / V9;
                    A3H_mat_stat_new(n) = C8_3 * A3H_mob_new(n) / V9;

                    % Calcul A- sur résine
                    A1_stat_new(n) = Keq_A1_HSO4 * A1_mob_new(n) * HSO4_stat_new(n) / HSO4_mob_new(n);
                    A2_stat_new(n) = Keq_A2_HSO4 * A2_mob_new(n) * HSO4_stat_new(n) / HSO4_mob_new(n);
                    A3_stat_new(n) = Keq_A3_HSO4 * A3_mob_new(n) * HSO4_stat_new(n) / HSO4_mob_new(n);

                    % Calcul Atot sur résine
                    A1tot_stat_new(n) = A1_stat_new(n) + A1H_SO4_stat_new(n) + A1H_HSO4_stat_new(n) + A1H_mat_stat_new(n);
                    A2tot_stat_new(n) = A2_stat_new(n) + A2H_SO4_stat_new(n) + A2H_HSO4_stat_new(n) + A2H_mat_stat_new(n);
                    A3tot_stat_new(n) = A3_stat_new(n) + A3H_SO4_stat_new(n) + A3H_HSO4_stat_new(n) + A3H_mat_stat_new(n);
                    
                    if n == NETtot
                        n = 1;
                    else
                        n = n+1;
                    end
                end

                % Mise à jour des volumes 
                V_inj_etape_new = V_inj_etape_new + Pas_V;
                V_inj_cycle_new = V_inj_cycle_new + Pas_V;
                V_inj_total_new = V_inj_total_new + Pas_V;

                % Mise à jour du graphique
                if (V_inj_sous_etape_new-V_inj_sous_etape_old) >= Pas_V_figure || V_inj_sous_etape_new == Pas_V || V_inj_sous_etape_new == V_alim
                    figure(1)
                    V_etape_fig = num2str(V_inj_etape_new,3);
                    V_sous_etape_fig = num2str(V_inj_sous_etape_new,3);
                    Cycle_fig = num2str(Cycle_new,3);
                    Etape_fig = num2str(Etape_new,3);
                    Sous_Etape_fig = num2str(Sous_Etape_new,3);

                    % H+, HSO4- et SO4 2- en phase mobile
                    subplot(4,1,1)
                    plot(Vecteur_new2,H_mob_new,'LineStyle','none','LineStyle','none','Marker','o','MarkerSize',2,'MarkerSize',2,'Color',[0.8 0.8 0.8]);
                    hold on
                    plot(Vecteur_new2,HSO4_mob_new,'LineStyle','none','LineStyle','none','Marker','o','MarkerSize',2,'MarkerSize',2,'Color',[0.8 0.8 0]);
                    plot(Vecteur_new2,SO4_mob_new,'LineStyle','none','LineStyle','none','Marker','o','MarkerSize',2,'MarkerSize',2,'Color',[0.8 0 0.8]);

                    Titre = ['Cycle = ',Cycle_fig,', Etape = ',Etape_fig,', Sous-Etape = ',Sous_Etape_fig,', Vétape = ',V_etape_fig,' BV, Vsous-étape = ',V_sous_etape_fig,' BV'];
                    title(Titre);
                    xlabel('Etage théorique');
                    ylabel('Concentration (mol/L)');
                    xlim([0 NETtot*1.09]);
                    Ymax_mob_1 = max(max([H_mob_new;HSO4_mob_new;SO4_mob_new]));
                    ylim([0 Ymax_mob_1]);

                    plot(Entree_Z1(Etape_new)*ones(3,1),linspace(0,Ymax_mob_1,3),'r--','LineWidth',2);
                    plot(Entree_Z2(Etape_new)*ones(3,1),linspace(0,Ymax_mob_1,3),'g--','LineWidth',2);
                    plot(Entree_Z3(Etape_new)*ones(3,1),linspace(0,Ymax_mob_1,3),'c--','LineWidth',2);
                    plot(Entree_Z4(Etape_new)*ones(3,1),linspace(0,Ymax_mob_1,3),'b--','LineWidth',2);
                    plot(Sortie_Z1(Etape_new)*ones(3,1),linspace(0,Ymax_mob_1,3),'r-.','MarkerSize',2);
                    plot(Sortie_Z2(Etape_new)*ones(3,1),linspace(0,Ymax_mob_1,3),'g-.','MarkerSize',2);
                    plot(Sortie_Z3(Etape_new)*ones(3,1),linspace(0,Ymax_mob_1,3),'c-.','MarkerSize',2);
                    plot(Sortie_Z4(Etape_new)*ones(3,1),linspace(0,Ymax_mob_1,3),'b-.','MarkerSize',2);

                    hold off
                    hleg1 = legend('H+','HSO4-','SO4 2-');
                    set(hleg1,'FontSize',9);
                    text(Milieu_Z1(Etape_new),Ymax_mob_1/2,'Z1','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','red')
                    text(Milieu_Z2(Etape_new),Ymax_mob_1/2,'Z2','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','green')
                    text(Milieu_Z3(Etape_new),Ymax_mob_1/2,'Z3','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','cyan')
                    text(Milieu_Z4(Etape_new),Ymax_mob_1/2,'Z4','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','blue')
 
                    if Sous_Etape_new == 1
                        text(Entree_Z3(Etape_new),0,'Produit \rightarrow','Rotation',90,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',9,'Color','cyan')
                        text(Sortie_Z3(Etape_new),0,'\rightarrow Raffinat','Rotation',-90,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',9,'Color','cyan')
                    elseif Sous_Etape_new == 2
                        text(Entree_Z1(Etape_new),0,'Eluant \rightarrow','Rotation',90,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',9,'Color','red')
                        text(Sortie_Z1(Etape_new),0,'\rightarrow Extrait','Rotation',-90,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',9,'Color','red')
                    elseif Sous_Etape_new == 3
                        text(Entree_Z1(Etape_new),0,'Eluant \rightarrow','Rotation',90,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',9,'Color','red')
                        text(Sortie_Z3(Etape_new),0,'\rightarrow Raffinat','Rotation',-90,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',9,'Color','cyan')
                    end
                    
                    % A1tot, A2tot, A3tot en phase mobile
                    subplot(4,1,2)
                    plot(Vecteur_new2,A1tot_mob_new,'LineStyle','none','Marker','o','MarkerSize',2,'Color',[1 0 0]);
                    hold on
                    plot(Vecteur_new2,A2tot_mob_new,'LineStyle','none','Marker','o','MarkerSize',2,'Color',[0 1 0]);
                    plot(Vecteur_new2,A3tot_mob_new,'LineStyle','none','Marker','o','MarkerSize',2,'Color',[0 0 1]);

                    xlabel('Etage théorique');
                    ylabel('Concentration (mol/L)');
                    xlim([0 NETtot*1.09]);
                    Ymax_mob_2 = max(max([A1tot_mob_new;A2tot_mob_new;A3tot_mob_new]));
                    if Ymax_mob_2>0
                        ylim([0 Ymax_mob_2]);
                    else
                        Ymax_mob_2 = Ymax_mob_1;
                        ylim([0 Ymax_mob_1]);
                    end

                    plot(Entree_Z1(Etape_new)*ones(3,1),linspace(0,Ymax_mob_2,3),'r--','LineWidth',2);
                    plot(Entree_Z2(Etape_new)*ones(3,1),linspace(0,Ymax_mob_2,3),'g--','LineWidth',2);
                    plot(Entree_Z3(Etape_new)*ones(3,1),linspace(0,Ymax_mob_2,3),'c--','LineWidth',2);
                    plot(Entree_Z4(Etape_new)*ones(3,1),linspace(0,Ymax_mob_2,3),'b--','LineWidth',2);
                    plot(Sortie_Z1(Etape_new)*ones(3,1),linspace(0,Ymax_mob_2,3),'r-.','MarkerSize',2);
                    plot(Sortie_Z2(Etape_new)*ones(3,1),linspace(0,Ymax_mob_2,3),'g-.','MarkerSize',2);
                    plot(Sortie_Z3(Etape_new)*ones(3,1),linspace(0,Ymax_mob_2,3),'c-.','MarkerSize',2);
                    plot(Sortie_Z4(Etape_new)*ones(3,1),linspace(0,Ymax_mob_2,3),'b-.','MarkerSize',2);

                    hold off
                    hleg2 = legend('A1 tot','A2 tot','A3 tot');
                    set(hleg2,'FontSize',9);
                    text(Milieu_Z1(Etape_new),Ymax_mob_2/2,'Z1','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','red')
                    text(Milieu_Z2(Etape_new),Ymax_mob_2/2,'Z2','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','green')
                    text(Milieu_Z3(Etape_new),Ymax_mob_2/2,'Z3','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','cyan')
                    text(Milieu_Z4(Etape_new),Ymax_mob_2/2,'Z4','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','blue')
    
                    if Sous_Etape_new == 1
                        text(Entree_Z3(Etape_new),0,'Produit \rightarrow','Rotation',90,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',9,'Color','cyan')
                        text(Sortie_Z3(Etape_new),0,'\rightarrow Raffinat','Rotation',-90,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',9,'Color','cyan')
                    elseif Sous_Etape_new == 2
                        text(Entree_Z1(Etape_new),0,'Eluant \rightarrow','Rotation',90,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',9,'Color','red')
                        text(Sortie_Z1(Etape_new),0,'\rightarrow Extrait','Rotation',-90,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',9,'Color','red')
                    elseif Sous_Etape_new == 3
                        text(Entree_Z1(Etape_new),0,'Eluant \rightarrow','Rotation',90,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',9,'Color','red')
                        text(Sortie_Z3(Etape_new),0,'\rightarrow Raffinat','Rotation',-90,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',9,'Color','cyan')                    
                    end           

                    % HSO4- et SO4 2- en phase stationnaire
                    subplot(4,1,3)
                    plot(Vecteur_new2,HSO4_stat_new,'LineStyle','none','Marker','o','MarkerSize',2,'Color',[0.8 0.8 0]);
                    hold on    
                    plot(Vecteur_new2,SO4_stat_new,'LineStyle','none','Marker','o','MarkerSize',2,'Color',[0.8 0 0.8]);     

                    xlabel('Etage théorique');
                    ylabel('Concentration (mol/Lrésine)');
                    xlim([0 NETtot*1.09]);
                    Ymax_stat_1 = max(max([HSO4_stat_new;SO4_stat_new]));
                    ylim([0 Ymax_stat_1]);

                    plot(Entree_Z1(Etape_new)*ones(3,1),linspace(0,Ymax_stat_1,3),'r--','LineWidth',2);
                    plot(Entree_Z2(Etape_new)*ones(3,1),linspace(0,Ymax_stat_1,3),'g--','LineWidth',2);
                    plot(Entree_Z3(Etape_new)*ones(3,1),linspace(0,Ymax_stat_1,3),'c--','LineWidth',2);
                    plot(Entree_Z4(Etape_new)*ones(3,1),linspace(0,Ymax_stat_1,3),'b--','LineWidth',2);
                    plot(Sortie_Z1(Etape_new)*ones(3,1),linspace(0,Ymax_stat_1,3),'r-.','MarkerSize',2);
                    plot(Sortie_Z2(Etape_new)*ones(3,1),linspace(0,Ymax_stat_1,3),'g-.','MarkerSize',2);
                    plot(Sortie_Z3(Etape_new)*ones(3,1),linspace(0,Ymax_stat_1,3),'c-.','MarkerSize',2);
                    plot(Sortie_Z4(Etape_new)*ones(3,1),linspace(0,Ymax_stat_1,3),'b-.','MarkerSize',2); 

                    hold off
                    hleg3 = legend('HSO4-','SO4 2-');
                    set(hleg3,'FontSize',9);
                    text(Milieu_Z1(Etape_new),Ymax_stat_1/2,'Z1','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','red')
                    text(Milieu_Z2(Etape_new),Ymax_stat_1/2,'Z2','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','green')
                    text(Milieu_Z3(Etape_new),Ymax_stat_1/2,'Z3','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','cyan')
                    text(Milieu_Z4(Etape_new),Ymax_stat_1/2,'Z4','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','blue')
    
                    if Sous_Etape_new == 1
                        text(Entree_Z3(Etape_new),0,'Produit \rightarrow','Rotation',90,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',9,'Color','cyan')
                        text(Sortie_Z3(Etape_new),0,'\rightarrow Raffinat','Rotation',-90,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',9,'Color','cyan')
                    elseif Sous_Etape_new == 2
                        text(Entree_Z1(Etape_new),0,'Eluant \rightarrow','Rotation',90,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',9,'Color','red')
                        text(Sortie_Z1(Etape_new),0,'\rightarrow Extrait','Rotation',-90,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',9,'Color','red')
                    elseif Sous_Etape_new == 3
                        text(Entree_Z1(Etape_new),0,'Eluant \rightarrow','Rotation',90,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',9,'Color','red')
                        text(Sortie_Z3(Etape_new),0,'\rightarrow Raffinat','Rotation',-90,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',9,'Color','cyan')                    
                    end            

                    % A1adsorbé, A1échangé, A2adsorbé, A2échangé, A3adsorbé et A3échangé en phase stationnaire
                    subplot(4,1,4)
                    A1H_stat=A1tot_stat_new-A1_stat_new;
                    A2H_stat=A2tot_stat_new-A2_stat_new;
                    A3H_stat=A3tot_stat_new-A3_stat_new;
                    plot(Vecteur_new2,A1H_stat,'LineStyle','none','Marker','o','MarkerSize',2,'Color',[1 0 0]);
                    hold on
                    plot(Vecteur_new2,A1_stat_new,'LineStyle','none','Marker','o','MarkerSize',2,'Color',[0.6 0 0]);
                    plot(Vecteur_new2,A2H_stat,'LineStyle','none','Marker','o','MarkerSize',2,'Color',[0 1 0]);
                    plot(Vecteur_new2,A2_stat_new,'LineStyle','none','Marker','o','MarkerSize',2,'Color',[0 0.6 0]);
                    plot(Vecteur_new2,A3H_stat,'LineStyle','none','Marker','o','MarkerSize',2,'Color',[0 0 1]);
                    plot(Vecteur_new2,A3_stat_new,'LineStyle','none','Marker','o','MarkerSize',2,'Color',[0 0 0.6]);

                    xlabel('Etage théorique');
                    ylabel('Concentration (mol/Lrésine)');
                    xlim([0 NETtot*1.09]);
                    Ymax_stat_2 = max(max([A1H_stat;A1_stat_new;A2H_stat;A2_stat_new;A3H_stat;A3_stat_new]));
                    if Ymax_stat_2>0
                        ylim([0 Ymax_stat_2]);
                    else
                        Ymax_stat_2 = Ymax_stat_1;
                        ylim([0 Ymax_stat_1]);
                    end

                    plot(Entree_Z1(Etape_new)*ones(3,1),linspace(0,Ymax_stat_2,3),'r--','LineWidth',2);
                    plot(Entree_Z2(Etape_new)*ones(3,1),linspace(0,Ymax_stat_2,3),'g--','LineWidth',2);
                    plot(Entree_Z3(Etape_new)*ones(3,1),linspace(0,Ymax_stat_2,3),'c--','LineWidth',2);
                    plot(Entree_Z4(Etape_new)*ones(3,1),linspace(0,Ymax_stat_2,3),'b--','LineWidth',2);
                    plot(Sortie_Z1(Etape_new)*ones(3,1),linspace(0,Ymax_stat_2,3),'r-.','MarkerSize',2);
                    plot(Sortie_Z2(Etape_new)*ones(3,1),linspace(0,Ymax_stat_2,3),'g-.','MarkerSize',2);
                    plot(Sortie_Z3(Etape_new)*ones(3,1),linspace(0,Ymax_stat_2,3),'c-.','MarkerSize',2);
                    plot(Sortie_Z4(Etape_new)*ones(3,1),linspace(0,Ymax_stat_2,3),'b-.','MarkerSize',2);

                    hold off
                    hleg4 = legend('A1H','A1-','A2H','A2-','A3H','A3-');
                    set(hleg4,'FontSize',9);
                    text(Milieu_Z1(Etape_new),Ymax_stat_2/2,'Z1','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','red')
                    text(Milieu_Z2(Etape_new),Ymax_stat_2/2,'Z2','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','green')
                    text(Milieu_Z3(Etape_new),Ymax_stat_2/2,'Z3','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','cyan')
                    text(Milieu_Z4(Etape_new),Ymax_stat_2/2,'Z4','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18,'Color','blue')
    
                    if Sous_Etape_new == 1
                        text(Entree_Z3(Etape_new),0,'Produit \rightarrow','Rotation',90,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',9,'Color','cyan')
                        text(Sortie_Z3(Etape_new),0,'\rightarrow Raffinat','Rotation',-90,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',9,'Color','cyan')
                    elseif Sous_Etape_new == 2
                        text(Entree_Z1(Etape_new),0,'Eluant \rightarrow','Rotation',90,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',9,'Color','red')
                        text(Sortie_Z1(Etape_new),0,'\rightarrow Extrait','Rotation',-90,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',9,'Color','red')
                    elseif Sous_Etape_new == 3
                        text(Entree_Z1(Etape_new),0,'Eluant \rightarrow','Rotation',90,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',9,'Color','red')
                        text(Sortie_Z3(Etape_new),0,'\rightarrow Raffinat','Rotation',-90,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',9,'Color','cyan')
                    end
                    
                    V_inj_sous_etape_old = V_inj_sous_etape_new;
                
                    pause(0.01)
                end

                % Mise à jour des résultats
                if Sous_Etape_new == 1
                    V_ech_total_new = V_ech_total_new + Pas_V;
                    V_ech_cycle_new = V_ech_cycle_new + Pas_V;
                    V_ech_etape_new = V_ech_etape_new + Pas_V;
                    V_ech_sous_etape_new = V_ech_sous_etape_new + Pas_V;
                    
                    H_raff_total_new = (H_raff_total_new * V_raff_total_new + Pas_V * H_mob_new(n_fin))/(V_raff_total_new + Pas_V);
                    H_raff_cycle_new = (H_raff_cycle_new * V_raff_cycle_new + Pas_V * H_mob_new(n_fin))/(V_raff_cycle_new + Pas_V);
                    H_raff_etape_new = (H_raff_etape_new * V_raff_etape_new + Pas_V * H_mob_new(n_fin))/(V_raff_etape_new + Pas_V);
                    H_raff_sous_etape_new = (H_raff_sous_etape_new * V_raff_sous_etape_new + Pas_V * H_mob_new(n_fin))/(V_raff_sous_etape_new + Pas_V);                    
                    HSO4_raff_total_new = (HSO4_raff_total_new * V_raff_total_new + Pas_V * HSO4_mob_new(n_fin))/(V_raff_total_new + Pas_V);
                    HSO4_raff_cycle_new = (HSO4_raff_cycle_new * V_raff_cycle_new + Pas_V * HSO4_mob_new(n_fin))/(V_raff_cycle_new + Pas_V);
                    HSO4_raff_etape_new = (HSO4_raff_etape_new * V_raff_etape_new + Pas_V * HSO4_mob_new(n_fin))/(V_raff_etape_new + Pas_V);
                    HSO4_raff_sous_etape_new = (HSO4_raff_sous_etape_new * V_raff_sous_etape_new + Pas_V * HSO4_mob_new(n_fin))/(V_raff_sous_etape_new + Pas_V);
                    SO4_raff_total_new = (SO4_raff_total_new * V_raff_total_new + Pas_V * SO4_mob_new(n_fin))/(V_raff_total_new + Pas_V);
                    SO4_raff_cycle_new = (SO4_raff_cycle_new * V_raff_cycle_new + Pas_V * SO4_mob_new(n_fin))/(V_raff_cycle_new + Pas_V);
                    SO4_raff_etape_new = (SO4_raff_etape_new * V_raff_etape_new + Pas_V * SO4_mob_new(n_fin))/(V_raff_etape_new + Pas_V);
                    SO4_raff_sous_etape_new = (SO4_raff_sous_etape_new * V_raff_sous_etape_new + Pas_V * SO4_mob_new(n_fin))/(V_raff_sous_etape_new + Pas_V);
                    A1H_raff_total_new = (A1H_raff_total_new * V_raff_total_new + Pas_V * A1H_mob_new(n_fin))/(V_raff_total_new + Pas_V);
                    A1H_raff_cycle_new = (A1H_raff_cycle_new * V_raff_cycle_new + Pas_V * A1H_mob_new(n_fin))/(V_raff_cycle_new + Pas_V);
                    A1H_raff_etape_new = (A1H_raff_etape_new * V_raff_etape_new + Pas_V * A1H_mob_new(n_fin))/(V_raff_etape_new + Pas_V);
                    A1H_raff_sous_etape_new = (A1H_raff_sous_etape_new * V_raff_sous_etape_new + Pas_V * A1H_mob_new(n_fin))/(V_raff_sous_etape_new + Pas_V);
                    A1_raff_total_new = (A1_raff_total_new * V_raff_total_new + Pas_V * A1_mob_new(n_fin))/(V_raff_total_new + Pas_V);
                    A1_raff_cycle_new =  (A1_raff_cycle_new * V_raff_cycle_new + Pas_V * A1_mob_new(n_fin))/(V_raff_cycle_new + Pas_V);
                    A1_raff_etape_new = (A1_raff_etape_new * V_raff_etape_new + Pas_V * A1_mob_new(n_fin))/(V_raff_etape_new + Pas_V);
                    A1_raff_sous_etape_new = (A1_raff_sous_etape_new * V_raff_sous_etape_new + Pas_V * A1_mob_new(n_fin))/(V_raff_sous_etape_new + Pas_V);
                    A2H_raff_total_new = (A2H_raff_total_new * V_raff_total_new + Pas_V * A2H_mob_new(n_fin))/(V_raff_total_new + Pas_V);
                    A2H_raff_cycle_new = (A2H_raff_cycle_new * V_raff_cycle_new + Pas_V * A2H_mob_new(n_fin))/(V_raff_cycle_new + Pas_V);
                    A2H_raff_etape_new = (A2H_raff_etape_new * V_raff_etape_new + Pas_V * A2H_mob_new(n_fin))/(V_raff_etape_new + Pas_V);
                    A2H_raff_sous_etape_new = (A2H_raff_sous_etape_new * V_raff_sous_etape_new + Pas_V * A2H_mob_new(n_fin))/(V_raff_sous_etape_new + Pas_V);
                    A2_raff_total_new = (A2_raff_total_new * V_raff_total_new + Pas_V * A2_mob_new(n_fin))/(V_raff_total_new + Pas_V);
                    A2_raff_cycle_new = (A2_raff_cycle_new * V_raff_cycle_new + Pas_V * A2_mob_new(n_fin))/(V_raff_cycle_new + Pas_V);
                    A2_raff_etape_new = (A2_raff_etape_new * V_raff_etape_new + Pas_V * A2_mob_new(n_fin))/(V_raff_etape_new + Pas_V);
                    A2_raff_sous_etape_new = (A2_raff_sous_etape_new * V_raff_sous_etape_new + Pas_V * A2_mob_new(n_fin))/(V_raff_sous_etape_new + Pas_V);
                    A3H_raff_total_new = (A3H_raff_total_new * V_raff_total_new + Pas_V * A3H_mob_new(n_fin))/(V_raff_total_new + Pas_V);
                    A3H_raff_cycle_new = (A3H_raff_cycle_new * V_raff_cycle_new + Pas_V * A3H_mob_new(n_fin))/(V_raff_cycle_new + Pas_V);
                    A3H_raff_etape_new = (A3H_raff_etape_new * V_raff_etape_new + Pas_V * A3H_mob_new(n_fin))/(V_raff_etape_new + Pas_V);
                    A3H_raff_sous_etape_new = (A3H_raff_sous_etape_new * V_raff_sous_etape_new + Pas_V * A3H_mob_new(n_fin))/(V_raff_sous_etape_new + Pas_V);
                    A3_raff_total_new = (A3_raff_total_new * V_raff_total_new + Pas_V * A3_mob_new(n_fin))/(V_raff_total_new + Pas_V);
                    A3_raff_cycle_new = (A3_raff_cycle_new * V_raff_cycle_new + Pas_V * A3_mob_new(n_fin))/(V_raff_cycle_new + Pas_V);
                    A3_raff_etape_new = (A3_raff_etape_new * V_raff_etape_new + Pas_V * A3_mob_new(n_fin))/(V_raff_etape_new + Pas_V);
                    A3_raff_sous_etape_new = (A3_raff_sous_etape_new * V_raff_sous_etape_new + Pas_V * A3_mob_new(n_fin))/(V_raff_sous_etape_new + Pas_V);                            
                    
                    V_raff_total_new = V_raff_total_new + Pas_V;
                    V_raff_cycle_new = V_raff_cycle_new + Pas_V;
                    V_raff_etape_new = V_raff_etape_new + Pas_V;
                    V_raff_sous_etape_new = V_raff_sous_etape_new + Pas_V;
                    
                elseif Sous_Etape_new == 2
                    V_elu_total_new = V_elu_total_new + Pas_V;
                    V_elu_cycle_new = V_elu_cycle_new + Pas_V;
                    V_elu_etape_new = V_elu_etape_new + Pas_V;
                    V_elu_sous_etape_new = V_elu_sous_etape_new + Pas_V;
                    
                    H_ext_total_new = (H_ext_total_new * V_ext_total_new + Pas_V * H_mob_new(n_fin))/(V_ext_total_new + Pas_V);
                    H_ext_cycle_new = (H_ext_cycle_new * V_ext_cycle_new + Pas_V * H_mob_new(n_fin))/(V_ext_cycle_new + Pas_V);
                    H_ext_etape_new = (H_ext_etape_new * V_ext_etape_new + Pas_V * H_mob_new(n_fin))/(V_ext_etape_new + Pas_V);
                    H_ext_sous_etape_new = (H_ext_sous_etape_new * V_ext_sous_etape_new + Pas_V * H_mob_new(n_fin))/(V_ext_sous_etape_new + Pas_V);                      
                    HSO4_ext_total_new = (HSO4_ext_total_new * V_ext_total_new + Pas_V * HSO4_mob_new(n_fin))/(V_ext_total_new + Pas_V);
                    HSO4_ext_cycle_new = (HSO4_ext_cycle_new * V_ext_cycle_new + Pas_V * HSO4_mob_new(n_fin))/(V_ext_cycle_new + Pas_V);
                    HSO4_ext_etape_new = (HSO4_ext_etape_new * V_ext_etape_new + Pas_V * HSO4_mob_new(n_fin))/(V_ext_etape_new + Pas_V);
                    HSO4_ext_sous_etape_new = (HSO4_ext_sous_etape_new * V_ext_sous_etape_new + Pas_V * HSO4_mob_new(n_fin))/(V_ext_sous_etape_new + Pas_V);
                    SO4_ext_total_new = (SO4_ext_total_new * V_ext_total_new + Pas_V * SO4_mob_new(n_fin))/(V_ext_total_new + Pas_V);
                    SO4_ext_cycle_new = (SO4_ext_cycle_new * V_ext_cycle_new + Pas_V * SO4_mob_new(n_fin))/(V_ext_cycle_new + Pas_V);
                    SO4_ext_etape_new = (SO4_ext_etape_new * V_ext_etape_new + Pas_V * SO4_mob_new(n_fin))/(V_ext_etape_new + Pas_V);
                    SO4_ext_sous_etape_new = (SO4_ext_sous_etape_new * V_ext_sous_etape_new + Pas_V * SO4_mob_new(n_fin))/(V_ext_sous_etape_new + Pas_V);
                    A1H_ext_total_new = (A1H_ext_total_new * V_ext_total_new + Pas_V * A1H_mob_new(n_fin))/(V_ext_total_new + Pas_V);
                    A1H_ext_cycle_new = (A1H_ext_cycle_new * V_ext_cycle_new + Pas_V * A1H_mob_new(n_fin))/(V_ext_cycle_new + Pas_V);
                    A1H_ext_etape_new = (A1H_ext_etape_new * V_ext_etape_new + Pas_V * A1H_mob_new(n_fin))/(V_ext_etape_new + Pas_V);
                    A1H_ext_sous_etape_new = (A1H_ext_sous_etape_new * V_ext_sous_etape_new + Pas_V * A1H_mob_new(n_fin))/(V_ext_sous_etape_new + Pas_V);
                    A1_ext_total_new = (A1_ext_total_new * V_ext_total_new + Pas_V * A1_mob_new(n_fin))/(V_ext_total_new + Pas_V);
                    A1_ext_cycle_new =  (A1_ext_cycle_new * V_ext_cycle_new + Pas_V * A1_mob_new(n_fin))/(V_ext_cycle_new + Pas_V);
                    A1_ext_etape_new = (A1_ext_etape_new * V_ext_etape_new + Pas_V * A1_mob_new(n_fin))/(V_ext_etape_new + Pas_V);
                    A1_ext_sous_etape_new = (A1_ext_sous_etape_new * V_ext_sous_etape_new + Pas_V * A1_mob_new(n_fin))/(V_ext_sous_etape_new + Pas_V);
                    A2H_ext_total_new = (A2H_ext_total_new * V_ext_total_new + Pas_V * A2H_mob_new(n_fin))/(V_ext_total_new + Pas_V);
                    A2H_ext_cycle_new = (A2H_ext_cycle_new * V_ext_cycle_new + Pas_V * A2H_mob_new(n_fin))/(V_ext_cycle_new + Pas_V);
                    A2H_ext_etape_new = (A2H_ext_etape_new * V_ext_etape_new + Pas_V * A2H_mob_new(n_fin))/(V_ext_etape_new + Pas_V);
                    A2H_ext_sous_etape_new = (A2H_ext_sous_etape_new * V_ext_sous_etape_new + Pas_V * A2H_mob_new(n_fin))/(V_ext_sous_etape_new + Pas_V);
                    A2_ext_total_new = (A2_ext_total_new * V_ext_total_new + Pas_V * A2_mob_new(n_fin))/(V_ext_total_new + Pas_V);
                    A2_ext_cycle_new = (A2_ext_cycle_new * V_ext_cycle_new + Pas_V * A2_mob_new(n_fin))/(V_ext_cycle_new + Pas_V);
                    A2_ext_etape_new = (A2_ext_etape_new * V_ext_etape_new + Pas_V * A2_mob_new(n_fin))/(V_ext_etape_new + Pas_V);
                    A2_ext_sous_etape_new = (A2_ext_sous_etape_new * V_ext_sous_etape_new + Pas_V * A2_mob_new(n_fin))/(V_ext_sous_etape_new + Pas_V);
                    A3H_ext_total_new = (A3H_ext_total_new * V_ext_total_new + Pas_V * A3H_mob_new(n_fin))/(V_ext_total_new + Pas_V);
                    A3H_ext_cycle_new = (A3H_ext_cycle_new * V_ext_cycle_new + Pas_V * A3H_mob_new(n_fin))/(V_ext_cycle_new + Pas_V);
                    A3H_ext_etape_new = (A3H_ext_etape_new * V_ext_etape_new + Pas_V * A3H_mob_new(n_fin))/(V_ext_etape_new + Pas_V);
                    A3H_ext_sous_etape_new = (A3H_ext_sous_etape_new * V_ext_sous_etape_new + Pas_V * A3H_mob_new(n_fin))/(V_ext_sous_etape_new + Pas_V);
                    A3_ext_total_new = (A3_ext_total_new * V_ext_total_new + Pas_V * A3_mob_new(n_fin))/(V_ext_total_new + Pas_V);
                    A3_ext_cycle_new = (A3_ext_cycle_new * V_ext_cycle_new + Pas_V * A3_mob_new(n_fin))/(V_ext_cycle_new + Pas_V);
                    A3_ext_etape_new = (A3_ext_etape_new * V_ext_etape_new + Pas_V * A3_mob_new(n_fin))/(V_ext_etape_new + Pas_V);
                    A3_ext_sous_etape_new = (A3_ext_sous_etape_new * V_ext_sous_etape_new + Pas_V * A3_mob_new(n_fin))/(V_ext_sous_etape_new + Pas_V);                    

                    V_ext_total_new = V_ext_total_new + Pas_V;
                    V_ext_cycle_new = V_ext_cycle_new + Pas_V;
                    V_ext_etape_new = V_ext_etape_new + Pas_V;
                    V_ext_sous_etape_new = V_ext_sous_etape_new + Pas_V;
                    
                elseif Sous_Etape_new == 3    
                    V_elu_total_new = V_elu_total_new + Pas_V;
                    V_elu_cycle_new = V_elu_cycle_new + Pas_V;
                    V_elu_etape_new = V_elu_etape_new + Pas_V;
                    V_elu_sous_etape_new = V_elu_sous_etape_new + Pas_V;
                    
                    H_raff_total_new = (H_raff_total_new * V_raff_total_new + Pas_V * H_mob_new(n_fin))/(V_raff_total_new + Pas_V);
                    H_raff_cycle_new = (H_raff_cycle_new * V_raff_cycle_new + Pas_V * H_mob_new(n_fin))/(V_raff_cycle_new + Pas_V);
                    H_raff_etape_new = (H_raff_etape_new * V_raff_etape_new + Pas_V * H_mob_new(n_fin))/(V_raff_etape_new + Pas_V);
                    H_raff_sous_etape_new = (H_raff_sous_etape_new * V_raff_sous_etape_new + Pas_V * H_mob_new(n_fin))/(V_raff_sous_etape_new + Pas_V);  
                    HSO4_raff_total_new = (HSO4_raff_total_new * V_raff_total_new + Pas_V * HSO4_mob_new(n_fin))/(V_raff_total_new + Pas_V);
                    HSO4_raff_cycle_new = (HSO4_raff_cycle_new * V_raff_cycle_new + Pas_V * HSO4_mob_new(n_fin))/(V_raff_cycle_new + Pas_V);
                    HSO4_raff_etape_new = (HSO4_raff_etape_new * V_raff_etape_new + Pas_V * HSO4_mob_new(n_fin))/(V_raff_etape_new + Pas_V);
                    HSO4_raff_sous_etape_new = (HSO4_raff_sous_etape_new * V_raff_sous_etape_new + Pas_V * HSO4_mob_new(n_fin))/(V_raff_sous_etape_new + Pas_V);
                    SO4_raff_total_new = (SO4_raff_total_new * V_raff_total_new + Pas_V * SO4_mob_new(n_fin))/(V_raff_total_new + Pas_V);
                    SO4_raff_cycle_new = (SO4_raff_cycle_new * V_raff_cycle_new + Pas_V * SO4_mob_new(n_fin))/(V_raff_cycle_new + Pas_V);
                    SO4_raff_etape_new = (SO4_raff_etape_new * V_raff_etape_new + Pas_V * SO4_mob_new(n_fin))/(V_raff_etape_new + Pas_V);
                    SO4_raff_sous_etape_new = (SO4_raff_sous_etape_new * V_raff_sous_etape_new + Pas_V * SO4_mob_new(n_fin))/(V_raff_sous_etape_new + Pas_V);
                    A1H_raff_total_new = (A1H_raff_total_new * V_raff_total_new + Pas_V * A1H_mob_new(n_fin))/(V_raff_total_new + Pas_V);
                    A1H_raff_cycle_new = (A1H_raff_cycle_new * V_raff_cycle_new + Pas_V * A1H_mob_new(n_fin))/(V_raff_cycle_new + Pas_V);
                    A1H_raff_etape_new = (A1H_raff_etape_new * V_raff_etape_new + Pas_V * A1H_mob_new(n_fin))/(V_raff_etape_new + Pas_V);
                    A1H_raff_sous_etape_new = (A1H_raff_sous_etape_new * V_raff_sous_etape_new + Pas_V * A1H_mob_new(n_fin))/(V_raff_sous_etape_new + Pas_V);
                    A1_raff_total_new = (A1_raff_total_new * V_raff_total_new + Pas_V * A1_mob_new(n_fin))/(V_raff_total_new + Pas_V);
                    A1_raff_cycle_new =  (A1_raff_cycle_new * V_raff_cycle_new + Pas_V * A1_mob_new(n_fin))/(V_raff_cycle_new + Pas_V);
                    A1_raff_etape_new = (A1_raff_etape_new * V_raff_etape_new + Pas_V * A1_mob_new(n_fin))/(V_raff_etape_new + Pas_V);
                    A1_raff_sous_etape_new = (A1_raff_sous_etape_new * V_raff_sous_etape_new + Pas_V * A1_mob_new(n_fin))/(V_raff_sous_etape_new + Pas_V);
                    A2H_raff_total_new = (A2H_raff_total_new * V_raff_total_new + Pas_V * A2H_mob_new(n_fin))/(V_raff_total_new + Pas_V);
                    A2H_raff_cycle_new = (A2H_raff_cycle_new * V_raff_cycle_new + Pas_V * A2H_mob_new(n_fin))/(V_raff_cycle_new + Pas_V);
                    A2H_raff_etape_new = (A2H_raff_etape_new * V_raff_etape_new + Pas_V * A2H_mob_new(n_fin))/(V_raff_etape_new + Pas_V);
                    A2H_raff_sous_etape_new = (A2H_raff_sous_etape_new * V_raff_sous_etape_new + Pas_V * A2H_mob_new(n_fin))/(V_raff_sous_etape_new + Pas_V);
                    A2_raff_total_new = (A2_raff_total_new * V_raff_total_new + Pas_V * A2_mob_new(n_fin))/(V_raff_total_new + Pas_V);
                    A2_raff_cycle_new = (A2_raff_cycle_new * V_raff_cycle_new + Pas_V * A2_mob_new(n_fin))/(V_raff_cycle_new + Pas_V);
                    A2_raff_etape_new = (A2_raff_etape_new * V_raff_etape_new + Pas_V * A2_mob_new(n_fin))/(V_raff_etape_new + Pas_V);
                    A2_raff_sous_etape_new = (A2_raff_sous_etape_new * V_raff_sous_etape_new + Pas_V * A2_mob_new(n_fin))/(V_raff_sous_etape_new + Pas_V);
                    A3H_raff_total_new = (A3H_raff_total_new * V_raff_total_new + Pas_V * A3H_mob_new(n_fin))/(V_raff_total_new + Pas_V);
                    A3H_raff_cycle_new = (A3H_raff_cycle_new * V_raff_cycle_new + Pas_V * A3H_mob_new(n_fin))/(V_raff_cycle_new + Pas_V);
                    A3H_raff_etape_new = (A3H_raff_etape_new * V_raff_etape_new + Pas_V * A3H_mob_new(n_fin))/(V_raff_etape_new + Pas_V);
                    A3H_raff_sous_etape_new = (A3H_raff_sous_etape_new * V_raff_sous_etape_new + Pas_V * A3H_mob_new(n_fin))/(V_raff_sous_etape_new + Pas_V);
                    A3_raff_total_new = (A3_raff_total_new * V_raff_total_new + Pas_V * A3_mob_new(n_fin))/(V_raff_total_new + Pas_V);
                    A3_raff_cycle_new = (A3_raff_cycle_new * V_raff_cycle_new + Pas_V * A3_mob_new(n_fin))/(V_raff_cycle_new + Pas_V);
                    A3_raff_etape_new = (A3_raff_etape_new * V_raff_etape_new + Pas_V * A3_mob_new(n_fin))/(V_raff_etape_new + Pas_V);
                    A3_raff_sous_etape_new = (A3_raff_sous_etape_new * V_raff_sous_etape_new + Pas_V * A3_mob_new(n_fin))/(V_raff_sous_etape_new + Pas_V);      
                    
                    V_raff_total_new = V_raff_total_new + Pas_V;
                    V_raff_cycle_new = V_raff_cycle_new + Pas_V;
                    V_raff_etape_new = V_raff_etape_new + Pas_V;
                    V_raff_sous_etape_new = V_raff_sous_etape_new + Pas_V;              
                end              
                    
                % Mémorisation définitive partielle des résultats
                if (V_inj_sous_etape_new-V_inj_sous_etape(i)) >= Pas_V_memoire || V_inj_sous_etape_new == Pas_V || V_inj_sous_etape_new == V_alim
                    i = i+1;
                    
                    V_inj_total(i) = V_inj_total_new;
                    V_inj_cycle(i) = V_inj_cycle_new;
                    V_inj_etape(i) = V_inj_etape_new;
                    V_inj_sous_etape(i) = V_inj_sous_etape_new;
                    V_ech_total(i) = V_ech_total_new;
                    V_ech_cycle(i) = V_ech_cycle_new;
                    V_ech_etape(i) = V_ech_etape_new;
                    V_ech_sous_etape(i) = V_ech_sous_etape_new;
                    V_elu_total(i) = V_elu_total_new;
                    V_elu_cycle(i) = V_elu_cycle_new;
                    V_elu_etape(i) = V_elu_etape_new;
                    V_elu_sous_etape(i) = V_elu_sous_etape_new;
                    V_ext_total(i) = V_ext_total_new;
                    V_ext_cycle(i) = V_ext_cycle_new;
                    V_ext_etape(i) = V_ext_etape_new;
                    V_ext_sous_etape(i) = V_ext_sous_etape_new;
                    V_raff_total(i) = V_raff_total_new;
                    V_raff_cycle(i) = V_raff_cycle_new;
                    V_raff_etape(i) = V_raff_etape_new;
                    V_raff_sous_etape(i) = V_raff_sous_etape_new;
                    
                    Cycle(i) = Cycle_new;
                    Etape(i) = Etape_new;
                    Sous_Etape(i) = Sous_Etape_new;

                    H_mob(i,:) = H_mob_new(1:Pas_n_memoire:NETtot);
                    HSO4_mob(i,:) = HSO4_mob_new(1:Pas_n_memoire:NETtot);
                    SO4_mob(i,:) = SO4_mob_new(1:Pas_n_memoire:NETtot);
                    SO4tot_mob(i,:) = SO4tot_mob_new(1:Pas_n_memoire:NETtot);
                    A1H_mob(i,:) = A1H_mob_new(1:Pas_n_memoire:NETtot);
                    A2H_mob(i,:) = A2H_mob_new(1:Pas_n_memoire:NETtot);
                    A3H_mob(i,:) = A3H_mob_new(1:Pas_n_memoire:NETtot);
                    A1_mob(i,:) = A1_mob_new(1:Pas_n_memoire:NETtot);
                    A2_mob(i,:) = A2_mob_new(1:Pas_n_memoire:NETtot);
                    A3_mob(i,:) = A3_mob_new(1:Pas_n_memoire:NETtot);
                    A1tot_mob(i,:) = A1tot_mob_new(1:Pas_n_memoire:NETtot);
                    A2tot_mob(i,:) = A2tot_mob_new(1:Pas_n_memoire:NETtot);
                    A3tot_mob(i,:) = A3tot_mob_new(1:Pas_n_memoire:NETtot);

                    HSO4_stat(i,:) = HSO4_stat_new(1:Pas_n_memoire:NETtot);
                    SO4_stat(i,:) = SO4_stat_new(1:Pas_n_memoire:NETtot);
                    SO4tot_stat(i,:) = SO4tot_stat_new(1:Pas_n_memoire:NETtot);
                    A1H_HSO4_stat(i,:) = A1H_HSO4_stat_new(1:Pas_n_memoire:NETtot);
                    A2H_HSO4_stat(i,:) = A2H_HSO4_stat_new(1:Pas_n_memoire:NETtot);
                    A3H_HSO4_stat(i,:) = A3H_HSO4_stat_new(1:Pas_n_memoire:NETtot);
                    A1H_SO4_stat(i,:) = A1H_SO4_stat_new(1:Pas_n_memoire:NETtot);
                    A2H_SO4_stat(i,:) = A2H_SO4_stat_new(1:Pas_n_memoire:NETtot);
                    A3H_SO4_stat(i,:) = A3H_SO4_stat_new(1:Pas_n_memoire:NETtot);
                    A1H_mat_stat(i,:) = A1H_mat_stat_new(1:Pas_n_memoire:NETtot);
                    A2H_mat_stat(i,:) = A2H_mat_stat_new(1:Pas_n_memoire:NETtot);
                    A3H_mat_stat(i,:) = A3H_mat_stat_new(1:Pas_n_memoire:NETtot);
                    A1_stat(i,:) = A1_stat_new(1:Pas_n_memoire:NETtot);
                    A2_stat(i,:) = A2_stat_new(1:Pas_n_memoire:NETtot);
                    A3_stat(i,:) = A3_stat_new(1:Pas_n_memoire:NETtot);
                    A1tot_stat(i,:) = A1tot_stat_new(1:Pas_n_memoire:NETtot);
                    A2tot_stat(i,:) = A2tot_stat_new(1:Pas_n_memoire:NETtot);
                    A3tot_stat(i,:) = A3tot_stat_new(1:Pas_n_memoire:NETtot);

                    F_ionique(i,:) = F_ionique_new(1:Pas_n_memoire:NETtot);
                    Activite_AH(i,:) = Activite_AH_new(1:Pas_n_memoire:NETtot);
                    Activite_HSO4(i,:) = Activite_HSO4_new(1:Pas_n_memoire:NETtot);

                    H_ext_total(i) = H_ext_total_new;
                    H_ext_cycle(i) = H_ext_cycle_new;
                    H_ext_etape(i) = H_ext_etape_new;
                    H_ext_sous_etape(i) = H_ext_sous_etape_new;                    
                    HSO4_ext_total(i) = HSO4_ext_total_new;
                    HSO4_ext_cycle(i) = HSO4_ext_cycle_new;
                    HSO4_ext_etape(i) = HSO4_ext_etape_new;
                    HSO4_ext_sous_etape(i) = HSO4_ext_sous_etape_new;
                    SO4_ext_total(i) = SO4_ext_total_new;
                    SO4_ext_cycle(i) = SO4_ext_cycle_new;
                    SO4_ext_etape(i) = SO4_ext_etape_new;
                    SO4_ext_sous_etape(i) = SO4_ext_sous_etape_new;
                    A1H_ext_total(i) = A1H_ext_total_new;
                    A1H_ext_cycle(i) = A1H_ext_cycle_new;
                    A1H_ext_etape(i) = A1H_ext_etape_new;
                    A1H_ext_sous_etape(i) = A1H_ext_sous_etape_new;
                    A1_ext_total(i) = A1_ext_total_new;
                    A1_ext_cycle(i) = A1_ext_cycle_new;
                    A1_ext_etape(i) = A1_ext_etape_new;
                    A1_ext_sous_etape(i) = A1_ext_sous_etape_new;
                    A2H_ext_total(i) = A2H_ext_total_new;
                    A2H_ext_cycle(i) = A2H_ext_cycle_new;
                    A2H_ext_etape(i) = A2H_ext_etape_new;
                    A2H_ext_sous_etape(i) = A2H_ext_sous_etape_new;
                    A2_ext_total(i) = A2_ext_total_new;
                    A2_ext_cycle(i) = A2_ext_cycle_new;
                    A2_ext_etape(i) = A2_ext_etape_new;
                    A2_ext_sous_etape(i) = A2_ext_sous_etape_new;
                    A3H_ext_total(i) = A3H_ext_total_new;
                    A3H_ext_cycle(i) = A3H_ext_cycle_new;
                    A3H_ext_etape(i) = A3H_ext_etape_new;
                    A3H_ext_sous_etape(i) = A3H_ext_sous_etape_new;
                    A3_ext_total(i) = A3_ext_total_new;
                    A3_ext_cycle(i) = A3_ext_cycle_new;
                    A3_ext_etape(i) = A3_ext_etape_new;
                    A3_ext_sous_etape(i) = A3_ext_sous_etape_new;

                    H_raff_total(i) = H_raff_total_new;
                    H_raff_cycle(i) = H_raff_cycle_new;
                    H_raff_etape(i) = H_raff_etape_new;
                    H_raff_sous_etape(i) = H_raff_sous_etape_new;                    
                    HSO4_raff_total(i) = HSO4_raff_total_new;
                    HSO4_raff_cycle(i) = HSO4_raff_cycle_new;
                    HSO4_raff_etape(i) = HSO4_raff_etape_new;
                    HSO4_raff_sous_etape(i) = HSO4_raff_sous_etape_new;
                    SO4_raff_total(i) = SO4_raff_total_new;
                    SO4_raff_cycle(i) = SO4_raff_cycle_new;
                    SO4_raff_etape(i) = SO4_raff_etape_new;
                    SO4_raff_sous_etape(i) = SO4_raff_sous_etape_new;
                    A1H_raff_total(i) = A1H_raff_total_new;
                    A1H_raff_cycle(i) = A1H_raff_cycle_new;
                    A1H_raff_etape(i) = A1H_raff_etape_new;
                    A1H_raff_sous_etape(i) = A1H_raff_sous_etape_new;
                    A1_raff_total(i) = A1_raff_total_new;
                    A1_raff_cycle(i) = A1_raff_cycle_new;
                    A1_raff_etape(i) = A1_raff_etape_new;
                    A1_raff_sous_etape(i) = A1_raff_sous_etape_new;
                    A2H_raff_total(i) = A2H_raff_total_new;
                    A2H_raff_cycle(i) = A2H_raff_cycle_new;
                    A2H_raff_etape(i) = A2H_raff_etape_new;
                    A2H_raff_sous_etape(i) = A2H_raff_sous_etape_new;
                    A2_raff_total(i) = A2_raff_total_new;
                    A2_raff_cycle(i) = A2_raff_cycle_new;
                    A2_raff_etape(i) = A2_raff_etape_new;
                    A2_raff_sous_etape(i) = A2_raff_sous_etape_new;
                    A3H_raff_total(i) = A3H_raff_total_new;
                    A3H_raff_cycle(i) = A3H_raff_cycle_new;
                    A3H_raff_etape(i) = A3H_raff_etape_new;
                    A3H_raff_sous_etape(i) = A3H_raff_sous_etape_new;
                    A3_raff_total(i) = A3_raff_total_new;
                    A3_raff_cycle(i) = A3_raff_cycle_new;
                    A3_raff_etape(i) = A3_raff_etape_new;
                    A3_raff_sous_etape(i) = A3_raff_sous_etape_new;
                  
                    N_iter_loc_1(i,:) = N_iter_loc_1_new;
                    N_iter_loc_1_tot(i) = N_iter_loc_1_tot_new;
                    N_iter_loc_2(i,:) = N_iter_loc_2_new;
                    N_iter_loc_2_tot(i) = N_iter_loc_2_tot_new;
                    N_iter_glob(i) = N_iter_glob_new;
                    Ecart_abs_max(i,:) = Ecart_abs_max_new;
                    Ecart_rel_max(i,:) = Ecart_rel_max_new;
                    Residu_abs_max(i,:) = Residu_abs_max_new;
                    Residu_rel_max(i,:) = Residu_rel_max_new;
                    Convergence(i,:) = Convergence_new;
                    Modele(i) = Modele_new;
                end
         
                % Mise à jour de l'ancienne composition des étages théoriques des différentes colonnes
                H_mob_old = H_mob_new;
                HSO4_mob_old = HSO4_mob_new;
                SO4_mob_old = SO4_mob_new;
                SO4tot_mob_old = SO4tot_mob_new;
                A1H_mob_old = A1H_mob_new;
                A2H_mob_old = A2H_mob_new;
                A3H_mob_old = A3H_mob_new;
                A1_mob_old = A1_mob_new;
                A2_mob_old = A2_mob_new;
                A3_mob_old = A3_mob_new;
                A1tot_mob_old = A1tot_mob_new;
                A2tot_mob_old = A2tot_mob_new;
                A3tot_mob_old = A3tot_mob_new; 
                
                HSO4_stat_old = HSO4_stat_new;
                SO4_stat_old = SO4_stat_new;
                SO4tot_stat_old = SO4tot_stat_new;
                A1H_HSO4_stat_old = A1H_HSO4_stat_new;
                A2H_HSO4_stat_old = A2H_HSO4_stat_new;
                A3H_HSO4_stat_old = A3H_HSO4_stat_new;
                A1H_SO4_stat_old = A1H_SO4_stat_new;
                A2H_SO4_stat_old = A2H_SO4_stat_new;
                A3H_SO4_stat_old = A3H_SO4_stat_new;
                A1H_mat_stat_old = A1H_mat_stat_new;
                A2H_mat_stat_old = A2H_mat_stat_new;
                A3H_mat_stat_old = A3H_mat_stat_new;
                A1_stat_old = A1_stat_new;
                A2_stat_old = A2_stat_new;
                A3_stat_old = A3_stat_new;
                A1tot_stat_old = A1tot_stat_new;
                A2tot_stat_old = A2tot_stat_new;
                A3tot_stat_old = A3tot_stat_new;             
            end
            
            Sous_Etape_new = Sous_Etape_new + 1;
        end
        
        Etape_new = Etape_new + 1;
    end

    % Optimisation des volumes d'alimentation lors des sous-étapes du SSMB
    if Cycle_new >= Ncycles_avant_opti && Etape_new > Ncol_tot
        n_fin_Z1 = Sortie_Z1(Ncol_tot,1);            
        n_deb_Z2 = Entree_Z2(Ncol_tot,1);
        n_fin_Z2 = Sortie_Z2(Ncol_tot,1);
        
        V_Alim = F_V_ech(Cycle_new)*V_ech_est;
        V_Ext = F_V_elu1(Cycle_new)*V_elu1_est;
        V_Raff = V_Alim+F_V_elu2(Cycle_new)*V_elu2_est;      
        
        n_penetration_lent = max([n_fin_Z1-NETcol+1;n_fin_Z1-ceil(Taux_penetration_Z4_lent_cible*(NETcol-1)/100)]);
        n_penetration_rapide_gauche = min([n_deb_Z2+NETcol-1;n_deb_Z2+ceil(Taux_penetration_Z1_gauche_rapide_cible*(NETcol-1)/100)]);
        n_penetration_rapide_droite = max([n_fin_Z2-NETcol+1;n_fin_Z2-ceil(Taux_penetration_Z1_droite_rapide_cible*(NETcol-1)/100)]);
        
        Quantite_lent_Alim = 0;
        Quantite_lent_Z4 = 0;            
        Quantite_lent_Ext = 0;
        Quantite_lent_Raff = 0;
        Correction_sels_elu = SO4tot_elu*(V_Raff+V_Ext-V_Alim)/(V_Raff+V_Ext);

        for Lent=1:length(Lents)
            if Lents(Lent) == 1
                Quantite_lent_Alim = Quantite_lent_Alim + V_Alim*(A1H_ech+A1_ech);
                Quantite_lent_Z4 = Quantite_lent_Z4 + (Porosite*sum(A1tot_mob_new(n_penetration_lent:n_fin_Z1))+(1-Porosite)*sum(A1tot_stat_new(n_penetration_lent:n_fin_Z1)))/(n_fin_Z1-n_penetration_lent+1);            
                Quantite_lent_Raff = Quantite_lent_Raff + V_Raff*(A1H_raff_etape_new+A1_raff_etape_new);
            elseif Lents(Lent) == 2
                Quantite_lent_Alim = Quantite_lent_Alim + V_Alim*(A2H_ech+A2_ech);
                Quantite_lent_Z4 = Quantite_lent_Z4 + (Porosite*sum(A2tot_mob_new(n_penetration_lent:n_fin_Z1))+(1-Porosite)*sum(A2tot_stat_new(n_penetration_lent:n_fin_Z1)))/(n_fin_Z1-n_penetration_lent+1);                    
                Quantite_lent_Raff = Quantite_lent_Raff + V_Raff*(A2H_raff_etape_new+A2_raff_etape_new);
            elseif Lents(Lent) == 3
                Quantite_lent_Alim = Quantite_lent_Alim + V_Alim*(A3H_ech+A3_ech);
                Quantite_lent_Z4 = Quantite_lent_Z4 + (Porosite*sum(A3tot_mob_new(n_penetration_lent:n_fin_Z1))+(1-Porosite)*sum(A3tot_stat_new(n_penetration_lent:n_fin_Z1)))/(n_fin_Z1-n_penetration_lent+1);                    
                Quantite_lent_Raff = Quantite_lent_Raff + V_Raff*(A3H_raff_etape_new+A3_raff_etape_new);
            elseif Lents(Lent) == 4
                Quantite_lent_Alim = Quantite_lent_Alim + V_Alim*(HSO4_ech+SO4_ech);
                Quantite_lent_Z4 = Quantite_lent_Z4 + Porosite*max([0,sum(SO4tot_mob_new(n_penetration_lent:n_fin_Z1))/(n_fin_Z1-n_penetration_lent+1)-Correction_sels_elu]);                      
                Quantite_lent_Raff = Quantite_lent_Raff + V_Raff*max([0,HSO4_raff_etape_new+SO4_raff_etape_new-Correction_sels_elu]);    
            end 
        end
        
        Quantite_rapide_Alim = 0;
        Quantite_rapide_Z1_gauche = 0;            
        Quantite_rapide_Z1_droite = 0;
        
        for Rapide=1:length(Rapides)
            if Rapides(Rapide) == 1
                Quantite_rapide_Alim = Quantite_rapide_Alim + V_Alim*(A1H_ech+A1_ech);
                Quantite_rapide_Z1_gauche = Quantite_rapide_Z1_gauche + (Porosite*sum(A1tot_mob_new(n_deb_Z2:n_penetration_rapide_gauche))+(1-Porosite)*sum(A1tot_stat_new(n_deb_Z2:n_penetration_rapide_gauche)))/(n_penetration_rapide_gauche-n_deb_Z2+1);
                Quantite_rapide_Z1_droite = Quantite_rapide_Z1_droite + (Porosite*sum(A1tot_mob_new(n_penetration_rapide_droite:n_fin_Z2))+(1-Porosite)*sum(A1tot_stat_new(n_penetration_rapide_droite:n_fin_Z2)))/(n_fin_Z2-n_penetration_rapide_droite+1); 
            elseif Rapides(Rapide) == 2
                Quantite_rapide_Alim = Quantite_rapide_Alim + V_Alim*(A2H_ech+A2_ech);
                Quantite_rapide_Z1_gauche = Quantite_rapide_Z1_gauche + (Porosite*sum(A2tot_mob_new(n_deb_Z2:n_penetration_rapide_gauche))+(1-Porosite)*sum(A2tot_stat_new(n_deb_Z2:n_penetration_rapide_gauche)))/(n_penetration_rapide_gauche-n_deb_Z2+1);
                Quantite_rapide_Z1_droite = Quantite_rapide_Z1_droite + (Porosite*sum(A2tot_mob_new(n_penetration_rapide_droite:n_fin_Z2))+(1-Porosite)*sum(A2tot_stat_new(n_penetration_rapide_droite:n_fin_Z2)))/(n_fin_Z2-n_penetration_rapide_droite+1);
            elseif Rapides(Rapide) == 3
                Quantite_rapide_Alim = Quantite_rapide_Alim + V_Alim*(A3H_ech+A3_ech);
                Quantite_rapide_Z1_gauche = Quantite_rapide_Z1_gauche + (Porosite*sum(A3tot_mob_new(n_deb_Z2:n_penetration_rapide_gauche))+(1-Porosite)*sum(A3tot_stat_new(n_deb_Z2:n_penetration_rapide_gauche)))/(n_penetration_rapide_gauche-n_deb_Z2+1);
                Quantite_rapide_Z1_droite = Quantite_rapide_Z1_droite + (Porosite*sum(A3tot_mob_new(n_penetration_rapide_droite:n_fin_Z2))+(1-Porosite)*sum(A3tot_stat_new(n_penetration_rapide_droite:n_fin_Z2)))/(n_fin_Z2-n_penetration_rapide_droite+1);
            elseif Rapides(Rapide) == 4
                Quantite_rapide_Alim = Quantite_rapide_Alim + V_Alim*(HSO4_ech+SO4_ech);
                Quantite_rapide_Z1_gauche = Quantite_rapide_Z1_gauche + Porosite*max([0,sum(SO4tot_mob_new(n_deb_Z2:n_penetration_rapide_gauche))/(n_penetration_rapide_gauche-n_deb_Z2+1)-Correction_sels_elu]);   
                Quantite_rapide_Z1_droite = Quantite_rapide_Z1_droite + Porosite*max([0,sum(SO4tot_mob_new(n_penetration_rapide_droite:n_fin_Z2))/(n_fin_Z2-n_penetration_rapide_droite+1)-Correction_sels_elu]); 
            end 
        end
        
        % Calcul du taux de pénétration des espèces lentes dans la future Z4
        Perte_lent_Z4_droite = 100*Quantite_lent_Z4/Quantite_lent_Alim;
        
        % Calcul de la perte totale d'espèces lentes dans le raffinat
        Perte_lent = 100*Quantite_lent_Raff/Quantite_lent_Alim;
        if Perte_lent > Perte_lent_Z4_droite
            Perte_lent_Z3_droite = Perte_lent-Perte_lent_Z4_droite;
        else
            Perte_lent_Z3_droite = Perte_lent_Z3_droite_cible;
        end
        
        % Calcul du taux de pénétration des espèces rapides à l'entrée de la future Z1
        Perte_rapide_Z1_gauche = 100*Quantite_rapide_Z1_gauche/Quantite_rapide_Alim;
        
        % Calcul du taux de pénétration des espèces rapides à la sortie de la future Z1
        Perte_rapide_Z1_droite = 100*Quantite_rapide_Z1_droite/Quantite_rapide_Alim;
        
        % Calcul des volumes d'alimentation pour la prochaine étape
        BV2_old = F_V_elu2(Cycle_new)*V_elu2_est + F_V_boucle(Cycle_new)*V_boucle_est;
        BV1_old = F_V_elu1(Cycle_new)*V_elu1_est + BV2_old;
        BV3_old = F_V_ech(Cycle_new)*V_ech_est + BV2_old;
        BV4_old = F_V_boucle(Cycle_new)*V_boucle_est;
        
        BV1_new = BV1_old*max(min((2-Perte_lent_Z4_droite_cible/Perte_lent_Z4_droite),1+DeltaBV_max/100),1-DeltaBV_max/100);
        BV3_new = BV3_old*max(min((Perte_lent_Z3_droite_cible/Perte_lent_Z3_droite),1+DeltaBV_max/100),1-DeltaBV_max/100);
        BV2_new = BV2_old*max(min((2-Perte_rapide_Z1_droite_cible/Perte_rapide_Z1_droite),1+DeltaBV_max/100),1-DeltaBV_max/100);
        BV4_new = BV4_old*max(min((Perte_rapide_Z1_gauche_cible/Perte_rapide_Z1_gauche),1+DeltaBV_max/100),1-DeltaBV_max/100);         

        % Calcul de l'éventuelle correction pour que les volumes soient positifs
        Correction_Vech = max([0,((BV3_old-BV2_old)-(BV3_new-BV2_new))/(BV3_old-BV2_old)]); 
        Correction_Velu1 = max([0,((BV1_old-BV2_old)-(BV1_new-BV2_new))/(BV1_old-BV2_old)]); 
        Correction_Velu2 = max([0,((BV2_old-BV4_old)-(BV2_new-BV4_new))/(BV2_old-BV4_old)]);
        Correction_Vboucle = max([0,(BV4_old-BV4_new)/BV4_old]);
        
        Correction = max([Correction_Vech,Correction_Velu1,Correction_Velu2,Correction_Vboucle]);
        if Correction > 0
            Correction = Correction/(1+Correction);
            
        else
            if Perte_lent_Z4_droite < Perte_lent_Z4_droite_limite
                BV1_new = BV1_old*max(min((2-Perte_lent_Z4_droite_limite/Perte_lent_Z4_droite),1+DeltaBV_max/100),1-DeltaBV_max/100);
            end
            if Perte_lent_Z3_droite < Perte_lent_Z3_droite_limite
                BV3_new = BV3_old*max(min((Perte_lent_Z3_droite_limite/Perte_lent_Z3_droite),1+DeltaBV_max/100),1-DeltaBV_max/100);
            end
            if Perte_rapide_Z1_droite < Perte_rapide_Z1_droite_limite
                BV2_new = BV2_old*max(min((2-Perte_rapide_Z1_droite_limite/Perte_rapide_Z1_droite),1+DeltaBV_max/100),1-DeltaBV_max/100);
            end
            if Perte_rapide_Z1_gauche < Perte_rapide_Z1_gauche_limite
                BV4_new = BV4_old*max(min((Perte_rapide_Z1_gauche_limite/Perte_rapide_Z1_gauche),1+DeltaBV_max/100),1-DeltaBV_max/100);
            end
            
            Correction_Vech = max([0,((BV3_old-BV2_old)-(BV3_new-BV2_new))/(BV3_old-BV2_old)]); 
            Correction_Velu1 = max([0,((BV1_old-BV2_old)-(BV1_new-BV2_new))/(BV1_old-BV2_old)]); 
            Correction_Velu2 = max([0,((BV2_old-BV4_old)-(BV2_new-BV4_new))/(BV2_old-BV4_old)]);
            Correction_Vboucle = max([0,(BV4_old-BV4_new)/BV4_old]);
            
            Correction = max([Correction_Vech,Correction_Velu1,Correction_Velu2,Correction_Vboucle]);           
            if Correction > 0
                Correction = Correction/(1+Correction);
            end
        end
        
        BV1_new = Correction*BV1_new + (1-Correction)*BV1_old;
        BV3_new = Correction*BV3_new + (1-Correction)*BV3_old;
        BV2_new = Correction*BV2_new + (1-Correction)*BV2_old;
        BV4_new = Correction*BV4_new + (1-Correction)*BV4_old;         
        
        F_V_ech_new = (BV3_new - BV2_new)/V_ech_est;
        F_V_elu1_new = (BV1_new - BV2_new)/V_elu1_est;
        F_V_elu2_new = (BV2_new - BV4_new)/V_elu2_est;
        F_V_boucle_new = BV4_new/V_boucle_est;

        % Mise à jour des volumes d'alimentation pour la prochaine étape
        if Cycle_new < Ncycles
            Correction_F(Cycle_new+1) = Correction;
            F_V_ech(Cycle_new+1) = F_V_ech_new;
            F_V_elu1(Cycle_new+1) = F_V_elu1_new;
            F_V_elu2(Cycle_new+1) = F_V_elu2_new; 
            F_V_boucle(Cycle_new+1) = F_V_boucle_new;                
        end
    end
        
    Cycle_new = Cycle_new + 1;      
end

% Calcul des performances de la SSMB
i_fin = i;
Purete_A1tot_ext_total = zeros(i_fin,1);
Purete_A2tot_ext_total = zeros(i_fin,1);
Purete_A3tot_ext_total = zeros(i_fin,1);
Purete_SO4tot_ext_total = zeros(i_fin,1);
Purete_A1tot_ext_cycle = zeros(i_fin,1);
Purete_A2tot_ext_cycle = zeros(i_fin,1);
Purete_A3tot_ext_cycle = zeros(i_fin,1);
Purete_SO4tot_ext_cycle = zeros(i_fin,1);
Purete_A1tot_ext_etape = zeros(i_fin,1);
Purete_A2tot_ext_etape = zeros(i_fin,1);
Purete_A3tot_ext_etape = zeros(i_fin,1);        
Purete_SO4tot_ext_etape = zeros(i_fin,1);
Purete_A1tot_ext_sous_etape = zeros(i_fin,1);
Purete_A2tot_ext_sous_etape = zeros(i_fin,1);
Purete_A3tot_ext_sous_etape = zeros(i_fin,1);
Purete_SO4tot_ext_sous_etape = zeros(i_fin,1);
Purete_A1tot_raff_total = zeros(i_fin,1);
Purete_A2tot_raff_total = zeros(i_fin,1);
Purete_A3tot_raff_total = zeros(i_fin,1);
Purete_SO4tot_raff_total = zeros(i_fin,1);
Purete_A1tot_raff_cycle = zeros(i_fin,1);
Purete_A2tot_raff_cycle = zeros(i_fin,1);
Purete_A3tot_raff_cycle = zeros(i_fin,1);
Purete_SO4tot_raff_cycle = zeros(i_fin,1);
Purete_A1tot_raff_etape = zeros(i_fin,1);
Purete_A2tot_raff_etape = zeros(i_fin,1);
Purete_A3tot_raff_etape = zeros(i_fin,1);        
Purete_SO4tot_raff_etape = zeros(i_fin,1);
Purete_A1tot_raff_sous_etape = zeros(i_fin,1);
Purete_A2tot_raff_sous_etape = zeros(i_fin,1);
Purete_A3tot_raff_sous_etape = zeros(i_fin,1);
Purete_SO4tot_raff_sous_etape = zeros(i_fin,1);

Recuperation_A1tot_ext_total = zeros(i_fin,1);
Recuperation_A2tot_ext_total = zeros(i_fin,1);
Recuperation_A3tot_ext_total = zeros(i_fin,1);
Recuperation_SO4tot_ext_total = zeros(i_fin,1);
Recuperation_A1tot_ext_cycle = zeros(i_fin,1);
Recuperation_A2tot_ext_cycle = zeros(i_fin,1);
Recuperation_A3tot_ext_cycle = zeros(i_fin,1);
Recuperation_SO4tot_ext_cycle = zeros(i_fin,1);
Recuperation_A1tot_ext_etape = zeros(i_fin,1);
Recuperation_A2tot_ext_etape = zeros(i_fin,1);
Recuperation_A3tot_ext_etape = zeros(i_fin,1);        
Recuperation_SO4tot_ext_etape = zeros(i_fin,1);
Recuperation_A1tot_ext_sous_etape = zeros(i_fin,1);
Recuperation_A2tot_ext_sous_etape = zeros(i_fin,1);
Recuperation_A3tot_ext_sous_etape = zeros(i_fin,1);
Recuperation_SO4tot_ext_sous_etape = zeros(i_fin,1);
Recuperation_A1tot_raff_total = zeros(i_fin,1);
Recuperation_A2tot_raff_total = zeros(i_fin,1);
Recuperation_A3tot_raff_total = zeros(i_fin,1);
Recuperation_SO4tot_raff_total = zeros(i_fin,1);
Recuperation_A1tot_raff_cycle = zeros(i_fin,1);
Recuperation_A2tot_raff_cycle = zeros(i_fin,1);
Recuperation_A3tot_raff_cycle = zeros(i_fin,1);
Recuperation_SO4tot_raff_cycle = zeros(i_fin,1);
Recuperation_A1tot_raff_etape = zeros(i_fin,1);
Recuperation_A2tot_raff_etape = zeros(i_fin,1);
Recuperation_A3tot_raff_etape = zeros(i_fin,1);        
Recuperation_SO4tot_raff_etape = zeros(i_fin,1);
Recuperation_A1tot_raff_sous_etape = zeros(i_fin,1);
Recuperation_A2tot_raff_sous_etape = zeros(i_fin,1);
Recuperation_A3tot_raff_sous_etape = zeros(i_fin,1);
Recuperation_SO4tot_raff_sous_etape = zeros(i_fin,1);

for i=2:i_fin
    
    % Pureté des composés dans l'extrait
    A1tot_ext_total = A1H_ext_total(i)+A1_ext_total(i);
    A2tot_ext_total = A2H_ext_total(i)+A2_ext_total(i);
    A3tot_ext_total = A3H_ext_total(i)+A3_ext_total(i);
    SO4tot_ext_total = HSO4_ext_total(i)+SO4_ext_total(i);
    Mtot_ext_total = A1tot_ext_total+A2tot_ext_total+A3tot_ext_total+SO4tot_ext_total;
    
    if Mtot_ext_total > 0
        Purete_A1tot_ext_total(i) = 100*A1tot_ext_total/Mtot_ext_total;
        Purete_A2tot_ext_total(i) = 100*A2tot_ext_total/Mtot_ext_total;
        Purete_A3tot_ext_total(i) = 100*A3tot_ext_total/Mtot_ext_total;
        Purete_SO4tot_ext_total(i) = 100*SO4tot_ext_total/Mtot_ext_total;
    else
        Purete_A1tot_ext_total(i) = 0;
        Purete_A2tot_ext_total(i) = 0;
        Purete_A3tot_ext_total(i) = 0;
        Purete_SO4tot_ext_total(i) = 0;
    end

    A1tot_ext_cycle = A1H_ext_cycle(i)+A1_ext_cycle(i);
    A2tot_ext_cycle = A2H_ext_cycle(i)+A2_ext_cycle(i);
    A3tot_ext_cycle = A3H_ext_cycle(i)+A3_ext_cycle(i);
    SO4tot_ext_cycle = HSO4_ext_cycle(i)+SO4_ext_cycle(i);
    Mtot_ext_cycle = A1tot_ext_cycle+A2tot_ext_cycle+A3tot_ext_cycle+SO4tot_ext_cycle;
    
    if Mtot_ext_cycle > 0
        Purete_A1tot_ext_cycle(i) = 100*A1tot_ext_cycle/Mtot_ext_cycle;
        Purete_A2tot_ext_cycle(i) = 100*A2tot_ext_cycle/Mtot_ext_cycle;
        Purete_A3tot_ext_cycle(i) = 100*A3tot_ext_cycle/Mtot_ext_cycle;
        Purete_SO4tot_ext_cycle(i) = 100*SO4tot_ext_cycle/Mtot_ext_cycle;
    else
        Purete_A1tot_ext_cycle(i) = 0;
        Purete_A2tot_ext_cycle(i) = 0;
        Purete_A3tot_ext_cycle(i) = 0;
        Purete_SO4tot_ext_cycle(i) = 0;
    end

    A1tot_ext_etape = A1H_ext_etape(i)+A1_ext_etape(i);
    A2tot_ext_etape = A2H_ext_etape(i)+A2_ext_etape(i);
    A3tot_ext_etape = A3H_ext_etape(i)+A3_ext_etape(i);
    SO4tot_ext_etape = HSO4_ext_etape(i)+SO4_ext_etape(i);
    Mtot_ext_etape = A1tot_ext_etape+A2tot_ext_etape+A3tot_ext_etape+SO4tot_ext_etape;
    
    if Mtot_ext_etape > 0
        Purete_A1tot_ext_etape(i) = 100*A1tot_ext_etape/Mtot_ext_etape;
        Purete_A2tot_ext_etape(i) = 100*A2tot_ext_etape/Mtot_ext_etape;
        Purete_A3tot_ext_etape(i) = 100*A3tot_ext_etape/Mtot_ext_etape;
        Purete_SO4tot_ext_etape(i) = 100*SO4tot_ext_etape/Mtot_ext_etape;
    else
        Purete_A1tot_ext_etape(i) = 0;
        Purete_A2tot_ext_etape(i) = 0;
        Purete_A3tot_ext_etape(i) = 0;
        Purete_SO4tot_ext_etape(i) = 0;
    end

    A1tot_ext_sous_etape = A1H_ext_sous_etape(i)+A1_ext_sous_etape(i);
    A2tot_ext_sous_etape = A2H_ext_sous_etape(i)+A2_ext_sous_etape(i);
    A3tot_ext_sous_etape = A3H_ext_sous_etape(i)+A3_ext_sous_etape(i);
    SO4tot_ext_sous_etape = HSO4_ext_sous_etape(i)+SO4_ext_sous_etape(i);
    Mtot_ext_sous_etape = A1tot_ext_sous_etape+A2tot_ext_sous_etape+A3tot_ext_sous_etape+SO4tot_ext_sous_etape;
    
    if Mtot_ext_sous_etape > 0
        Purete_A1tot_ext_sous_etape(i) = 100*A1tot_ext_sous_etape/Mtot_ext_sous_etape;
        Purete_A2tot_ext_sous_etape(i) = 100*A2tot_ext_sous_etape/Mtot_ext_sous_etape;
        Purete_A3tot_ext_sous_etape(i) = 100*A3tot_ext_sous_etape/Mtot_ext_sous_etape;
        Purete_SO4tot_ext_sous_etape(i) = 100*SO4tot_ext_sous_etape/Mtot_ext_sous_etape;
    else
        Purete_A1tot_ext_sous_etape(i) = 0;
        Purete_A2tot_ext_sous_etape(i) = 0;
        Purete_A3tot_ext_sous_etape(i) = 0;
        Purete_SO4tot_ext_sous_etape(i) = 0;
    end

    % Pureté des composés dans le raffinat    
    A1tot_raff_total = A1H_raff_total(i)+A1_raff_total(i);
    A2tot_raff_total = A2H_raff_total(i)+A2_raff_total(i);
    A3tot_raff_total = A3H_raff_total(i)+A3_raff_total(i);
    SO4tot_raff_total = HSO4_raff_total(i)+SO4_raff_total(i);
    Mtot_raff_total = A1tot_raff_total+A2tot_raff_total+A3tot_raff_total+SO4tot_raff_total;
    
    if Mtot_raff_total > 0
        Purete_A1tot_raff_total(i) = 100*A1tot_raff_total/Mtot_raff_total;
        Purete_A2tot_raff_total(i) = 100*A2tot_raff_total/Mtot_raff_total;
        Purete_A3tot_raff_total(i) = 100*A3tot_raff_total/Mtot_raff_total;
        Purete_SO4tot_raff_total(i) = 100*SO4tot_raff_total/Mtot_raff_total;
    else
        Purete_A1tot_raff_total(i) = 0;
        Purete_A2tot_raff_total(i) = 0;
        Purete_A3tot_raff_total(i) = 0;
        Purete_SO4tot_raff_total(i) = 0;
    end

    A1tot_raff_cycle = A1H_raff_cycle(i)+A1_raff_cycle(i);
    A2tot_raff_cycle = A2H_raff_cycle(i)+A2_raff_cycle(i);
    A3tot_raff_cycle = A3H_raff_cycle(i)+A3_raff_cycle(i);
    SO4tot_raff_cycle = HSO4_raff_cycle(i)+SO4_raff_cycle(i);
    Mtot_raff_cycle = A1tot_raff_cycle+A2tot_raff_cycle+A3tot_raff_cycle+SO4tot_raff_cycle;
    
    if Mtot_raff_cycle > 0
        Purete_A1tot_raff_cycle(i) = 100*A1tot_raff_cycle/Mtot_raff_cycle;
        Purete_A2tot_raff_cycle(i) = 100*A2tot_raff_cycle/Mtot_raff_cycle;
        Purete_A3tot_raff_cycle(i) = 100*A3tot_raff_cycle/Mtot_raff_cycle;
        Purete_SO4tot_raff_cycle(i) = 100*SO4tot_raff_cycle/Mtot_raff_cycle;
    else
        Purete_A1tot_raff_cycle(i) = 0;
        Purete_A2tot_raff_cycle(i) = 0;
        Purete_A3tot_raff_cycle(i) = 0;
        Purete_SO4tot_raff_cycle(i) = 0;
    end

    A1tot_raff_etape = A1H_raff_etape(i)+A1_raff_etape(i);
    A2tot_raff_etape = A2H_raff_etape(i)+A2_raff_etape(i);
    A3tot_raff_etape = A3H_raff_etape(i)+A3_raff_etape(i);
    SO4tot_raff_etape = HSO4_raff_etape(i)+SO4_raff_etape(i);
    Mtot_raff_etape = A1tot_raff_etape+A2tot_raff_etape+A3tot_raff_etape+SO4tot_raff_etape;
    
    if Mtot_raff_etape > 0
        Purete_A1tot_raff_etape(i) = 100*A1tot_raff_etape/Mtot_raff_etape;
        Purete_A2tot_raff_etape(i) = 100*A2tot_raff_etape/Mtot_raff_etape;
        Purete_A3tot_raff_etape(i) = 100*A3tot_raff_etape/Mtot_raff_etape;
        Purete_SO4tot_raff_etape(i) = 100*SO4tot_raff_etape/Mtot_raff_etape;
    else
        Purete_A1tot_raff_etape(i) = 0;
        Purete_A2tot_raff_etape(i) = 0;
        Purete_A3tot_raff_etape(i) = 0;
        Purete_SO4tot_raff_etape(i) = 0;
    end

    A1tot_raff_sous_etape = A1H_raff_sous_etape(i)+A1_raff_sous_etape(i);
    A2tot_raff_sous_etape = A2H_raff_sous_etape(i)+A2_raff_sous_etape(i);
    A3tot_raff_sous_etape = A3H_raff_sous_etape(i)+A3_raff_sous_etape(i);
    SO4tot_raff_sous_etape = HSO4_raff_sous_etape(i)+SO4_raff_sous_etape(i);
    Mtot_raff_sous_etape = A1tot_raff_sous_etape+A2tot_raff_sous_etape+A3tot_raff_sous_etape+SO4tot_raff_sous_etape;
    
    if Mtot_raff_sous_etape > 0
        Purete_A1tot_raff_sous_etape(i) = 100*A1tot_raff_sous_etape/Mtot_raff_sous_etape;
        Purete_A2tot_raff_sous_etape(i) = 100*A2tot_raff_sous_etape/Mtot_raff_sous_etape;
        Purete_A3tot_raff_sous_etape(i) = 100*A3tot_raff_sous_etape/Mtot_raff_sous_etape;
        Purete_SO4tot_raff_sous_etape(i) = 100*SO4tot_raff_sous_etape/Mtot_raff_sous_etape;
    else
        Purete_A1tot_raff_sous_etape(i) = 0;
        Purete_A2tot_raff_sous_etape(i) = 0;
        Purete_A3tot_raff_sous_etape(i) = 0;
        Purete_SO4tot_raff_sous_etape(i) = 0;
    end
    
    % Taux de récupération des composés dans l'extrait et dans le raffinat
    V_Ech = F_V_ech(Cycle(i))*V_ech_est;
    V_Ext = V_Ech + F_V_elu1(Cycle(i))*V_elu1_est;
    V_Raff = F_V_elu2(Cycle(i))*V_elu2_est;
    Correction_sels_elu = SO4tot_elu*(V_Raff+V_Ext-V_Ech)/(V_Raff+V_Ext);    
    
    A1tot_alim_total = V_ech_total(i)*A1tot_ech;
    if A1tot_alim_total > 0
        Recuperation_A1tot_ext_total(i) = 100*V_ext_total(i)*(A1H_ext_total(i)+A1_ext_total(i))/A1tot_alim_total;
        Recuperation_A1tot_raff_total(i) = 100* V_raff_total(i)*(A1H_raff_total(i)+A1_raff_total(i))/A1tot_alim_total;
    else
        Recuperation_A1tot_ext_total(i) = 0;
        Recuperation_A1tot_raff_total(i) = 0;
    end
    
    A2tot_alim_total = V_ech_total(i)*A2tot_ech;
    if A2tot_alim_total > 0
        Recuperation_A2tot_ext_total(i) = 100*V_ext_total(i)*(A2H_ext_total(i)+A2_ext_total(i))/A2tot_alim_total;
        Recuperation_A2tot_raff_total(i) = 100* V_raff_total(i)*(A2H_raff_total(i)+A2_raff_total(i))/A2tot_alim_total;
    else
        Recuperation_A2tot_ext_total(i) = 0;
        Recuperation_A2tot_raff_total(i) = 0;
    end
    
    A3tot_alim_total = V_ech_total(i)*A3tot_ech;
    if A3tot_alim_total > 0
        Recuperation_A3tot_ext_total(i) = 100*V_ext_total(i)*(A3H_ext_total(i)+A3_ext_total(i))/A3tot_alim_total;
        Recuperation_A3tot_raff_total(i) = 100*V_raff_total(i)*(A3H_raff_total(i)+A3_raff_total(i))/A3tot_alim_total;
    else
        Recuperation_A3tot_ext_total(i) = 0;
        Recuperation_A3tot_raff_total(i) = 0;
    end

    SO4tot_alim_total = V_ech_total(i)*SO4tot_ech;
    if SO4tot_alim_total > 0
        Recuperation_SO4tot_ext_total(i) = 100*V_ext_total(i)*(HSO4_ext_total(i)+SO4_ext_total(i)-Correction_sels_elu)/SO4tot_alim_total;
        Recuperation_SO4tot_raff_total(i) = 100*(V_raff_total(i)*(HSO4_raff_total(i)+SO4_raff_total(i))-(V_raff_total(i)-V_ech_total(i))*Correction_sels_elu)/SO4tot_alim_total;
    else
        Recuperation_SO4tot_ext_total(i) = 0;
        Recuperation_SO4tot_raff_total(i) = 0;
    end    
    
    A1tot_alim_cycle = V_ech_cycle(i)*A1tot_ech;
    if A1tot_alim_cycle > 0
        Recuperation_A1tot_ext_cycle(i) = 100*V_ext_cycle(i)*(A1H_ext_cycle(i)+A1_ext_cycle(i))/A1tot_alim_cycle;
        Recuperation_A1tot_raff_cycle(i) = 100* V_raff_cycle(i)*(A1H_raff_cycle(i)+A1_raff_cycle(i))/A1tot_alim_cycle;
    else
        Recuperation_A1tot_ext_cycle(i) = 0;
        Recuperation_A1tot_raff_cycle(i) = 0;
    end
    
    A2tot_alim_cycle = V_ech_cycle(i)*A2tot_ech;
    if A2tot_alim_cycle > 0
        Recuperation_A2tot_ext_cycle(i) = 100*V_ext_cycle(i)*(A2H_ext_cycle(i)+A2_ext_cycle(i))/A2tot_alim_cycle;
        Recuperation_A2tot_raff_cycle(i) = 100* V_raff_cycle(i)*(A2H_raff_cycle(i)+A2_raff_cycle(i))/A2tot_alim_cycle;
    else
        Recuperation_A2tot_ext_cycle(i) = 0;
        Recuperation_A2tot_raff_cycle(i) = 0;
    end
    
    A3tot_alim_cycle = V_ech_cycle(i)*A3tot_ech;
    if A3tot_alim_cycle > 0
        Recuperation_A3tot_ext_cycle(i) = 100*V_ext_cycle(i)*(A3H_ext_cycle(i)+A3_ext_cycle(i))/A3tot_alim_cycle;
        Recuperation_A3tot_raff_cycle(i) = 100*V_raff_cycle(i)*(A3H_raff_cycle(i)+A3_raff_cycle(i))/A3tot_alim_cycle;
    else
        Recuperation_A3tot_ext_cycle(i) = 0;
        Recuperation_A3tot_raff_cycle(i) = 0;
    end

    SO4tot_alim_cycle = V_ech_cycle(i)*SO4tot_ech;
    if SO4tot_alim_cycle > 0
        Recuperation_SO4tot_ext_cycle(i) = 100*V_ext_cycle(i)*(HSO4_ext_cycle(i)+SO4_ext_cycle(i)-Correction_sels_elu)/SO4tot_alim_cycle;
        Recuperation_SO4tot_raff_cycle(i) = 100*(V_raff_cycle(i)*(HSO4_raff_cycle(i)+SO4_raff_cycle(i))-(V_raff_cycle(i)-V_ech_cycle(i))*Correction_sels_elu)/SO4tot_alim_cycle;
    else
        Recuperation_SO4tot_ext_cycle(i) = 0;
        Recuperation_SO4tot_raff_cycle(i) = 0;
    end
    
    A1tot_alim_etape = V_ech_etape(i)*A1tot_ech;
    if A1tot_alim_etape > 0
        Recuperation_A1tot_ext_etape(i) = 100*V_ext_etape(i)*(A1H_ext_etape(i)+A1_ext_etape(i))/A1tot_alim_etape;
        Recuperation_A1tot_raff_etape(i) = 100* V_raff_etape(i)*(A1H_raff_etape(i)+A1_raff_etape(i))/A1tot_alim_etape;
    else
        Recuperation_A1tot_ext_etape(i) = 0;
        Recuperation_A1tot_raff_etape(i) = 0;
    end
    
    A2tot_alim_etape = V_ech_etape(i)*A2tot_ech;
    if A2tot_alim_etape > 0
        Recuperation_A2tot_ext_etape(i) = 100*V_ext_etape(i)*(A2H_ext_etape(i)+A2_ext_etape(i))/A2tot_alim_etape;
        Recuperation_A2tot_raff_etape(i) = 100* V_raff_etape(i)*(A2H_raff_etape(i)+A2_raff_etape(i))/A2tot_alim_etape;
    else
        Recuperation_A2tot_ext_etape(i) = 0;
        Recuperation_A2tot_raff_etape(i) = 0;
    end
    
    A3tot_alim_etape = V_ech_etape(i)*A3tot_ech;
    if A3tot_alim_etape > 0
        Recuperation_A3tot_ext_etape(i) = 100*V_ext_etape(i)*(A3H_ext_etape(i)+A3_ext_etape(i))/A3tot_alim_etape;
        Recuperation_A3tot_raff_etape(i) = 100*V_raff_etape(i)*(A3H_raff_etape(i)+A3_raff_etape(i))/A3tot_alim_etape;
    else
        Recuperation_A3tot_ext_etape(i) = 0;
        Recuperation_A3tot_raff_etape(i) = 0;
    end

    SO4tot_alim_etape = V_ech_etape(i)*SO4tot_ech;
    if SO4tot_alim_etape > 0
        Recuperation_SO4tot_ext_etape(i) = 100*V_ext_etape(i)*(HSO4_ext_etape(i)+SO4_ext_etape(i)-Correction_sels_elu)/SO4tot_alim_etape;
        Recuperation_SO4tot_raff_etape(i) = 100*(V_raff_etape(i)*(HSO4_raff_etape(i)+SO4_raff_etape(i))-(V_raff_etape(i)-V_ech_etape(i))*Correction_sels_elu)/SO4tot_alim_etape;
    else
        Recuperation_SO4tot_ext_etape(i) = 0;
        Recuperation_SO4tot_raff_etape(i) = 0;
    end
    
    A1tot_alim_sous_etape = V_ech_sous_etape(i)*A1tot_ech;
    if A1tot_alim_sous_etape > 0
        Recuperation_A1tot_ext_sous_etape(i) = 100*V_ext_sous_etape(i)*(A1H_ext_sous_etape(i)+A1_ext_sous_etape(i))/A1tot_alim_sous_etape;
        Recuperation_A1tot_raff_sous_etape(i) = 100* V_raff_sous_etape(i)*(A1H_raff_sous_etape(i)+A1_raff_sous_etape(i))/A1tot_alim_sous_etape;
    else
        Recuperation_A1tot_ext_sous_etape(i) = 0;
        Recuperation_A1tot_raff_sous_etape(i) = 0;
    end
    
    A2tot_alim_sous_etape = V_ech_sous_etape(i)*A2tot_ech;
    if A2tot_alim_sous_etape > 0
        Recuperation_A2tot_ext_sous_etape(i) = 100*V_ext_sous_etape(i)*(A2H_ext_sous_etape(i)+A2_ext_sous_etape(i))/A2tot_alim_sous_etape;
        Recuperation_A2tot_raff_sous_etape(i) = 100* V_raff_sous_etape(i)*(A2H_raff_sous_etape(i)+A2_raff_sous_etape(i))/A2tot_alim_sous_etape;
    else
        Recuperation_A2tot_ext_sous_etape(i) = 0;
        Recuperation_A2tot_raff_sous_etape(i) = 0;
    end
    
    A3tot_alim_sous_etape = V_ech_sous_etape(i)*A3tot_ech;
    if A3tot_alim_sous_etape > 0
        Recuperation_A3tot_ext_sous_etape(i) = 100*V_ext_sous_etape(i)*(A3H_ext_sous_etape(i)+A3_ext_sous_etape(i))/A3tot_alim_sous_etape;
        Recuperation_A3tot_raff_sous_etape(i) = 100*V_raff_sous_etape(i)*(A3H_raff_sous_etape(i)+A3_raff_sous_etape(i))/A3tot_alim_sous_etape;
    else
        Recuperation_A3tot_ext_sous_etape(i) = 0;
        Recuperation_A3tot_raff_sous_etape(i) = 0;
    end

    SO4tot_alim_sous_etape = V_ech_sous_etape(i)*SO4tot_ech;
    if SO4tot_alim_sous_etape > 0
        Recuperation_SO4tot_ext_sous_etape(i) = 100*V_ext_sous_etape(i)*(HSO4_ext_sous_etape(i)+SO4_ext_sous_etape(i)-Correction_sels_elu)/SO4tot_alim_sous_etape;
        Recuperation_SO4tot_raff_sous_etape(i) = 100*(V_raff_sous_etape(i)*(HSO4_raff_sous_etape(i)+SO4_raff_sous_etape(i))-(V_raff_sous_etape(i)-V_ech_sous_etape(i))*Correction_sels_elu)/SO4tot_alim_sous_etape;
    else
        Recuperation_SO4tot_ext_sous_etape(i) = 0;
        Recuperation_SO4tot_raff_sous_etape(i) = 0;
    end
end                      

Bilan_A1_total = Recuperation_A1tot_ext_total+Recuperation_A1tot_raff_total;
Bilan_A2_total = Recuperation_A2tot_ext_total+Recuperation_A2tot_raff_total;
Bilan_A3_total = Recuperation_A3tot_ext_total+Recuperation_A3tot_raff_total;
Bilan_A1_cycle = Recuperation_A1tot_ext_cycle+Recuperation_A1tot_raff_cycle;
Bilan_A2_cycle = Recuperation_A2tot_ext_cycle+Recuperation_A2tot_raff_cycle;
Bilan_A3_cycle = Recuperation_A3tot_ext_cycle+Recuperation_A3tot_raff_cycle;
Bilan_A1_etape = Recuperation_A1tot_ext_etape+Recuperation_A1tot_raff_etape;
Bilan_A2_etape = Recuperation_A2tot_ext_etape+Recuperation_A2tot_raff_etape;
Bilan_A3_etape = Recuperation_A3tot_ext_etape+Recuperation_A3tot_raff_etape;
Bilan_A1_sous_etape = Recuperation_A1tot_ext_sous_etape+Recuperation_A1tot_raff_sous_etape;
Bilan_A2_sous_etape = Recuperation_A2tot_ext_sous_etape+Recuperation_A2tot_raff_sous_etape;
Bilan_A3_sous_etape = Recuperation_A3tot_ext_sous_etape+Recuperation_A3tot_raff_sous_etape;

% Création graphique : évolution des performances de la SSMB
figure(2)

    % Pureté des acides organiques dans l'extrait en fonction du volume total injecté
    subplot(2,1,1)
    plot(V_inj_total(1:i_fin),Purete_A1tot_ext_total(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[1 0 0]);
    hold on
    plot(V_inj_total(1:i_fin),Purete_A2tot_ext_total(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 1 0]);
    plot(V_inj_total(1:i_fin),Purete_A3tot_ext_total(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0 1]);
    plot(V_inj_total(1:i_fin),Purete_SO4tot_ext_total(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 1 1]);
    plot(V_inj_total(1:i_fin),Purete_A1tot_ext_cycle(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0.7 0 0]);
    plot(V_inj_total(1:i_fin),Purete_A2tot_ext_cycle(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0.7 0]);
    plot(V_inj_total(1:i_fin),Purete_A3tot_ext_cycle(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0 0.7]);
    plot(V_inj_total(1:i_fin),Purete_SO4tot_ext_cycle(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0.7 0.7]);
    
    if Ncycles*Ncol_tot <= 12
        plot(V_inj_total(1:i_fin),Purete_A1tot_ext_etape(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0.4 0 0]);
        plot(V_inj_total(1:i_fin),Purete_A2tot_ext_etape(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0.4 0]);
        plot(V_inj_total(1:i_fin),Purete_A3tot_ext_etape(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0 0.4]);
        plot(V_inj_total(1:i_fin),Purete_SO4tot_ext_etape(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0.4 0.4]);
    end
    
    for i=2:i_fin
        V_ech = F_V_ech(Cycle(i))*V_ech_est;
        V_elu1 = F_V_elu1(Cycle(i))*V_elu1_est;
        V_elu2 = F_V_elu2(Cycle(i))*V_elu2_est;
        V_boucle = F_V_boucle(Cycle(i))*V_boucle_est;
        
        if Ncycles == 1
            if V_inj_sous_etape(i) == V_ech && Sous_Etape(i) == 1
                plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k--','LineWidth',1);
                text(V_inj_total(i)-V_ech/2,100,'Feed','HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','red');

            elseif V_inj_sous_etape(i) == V_elu1 && Sous_Etape(i) == 2 
                plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k--','LineWidth',1);
                text(V_inj_total(i)-V_elu1/2,100,'Extrait','HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','green');

            elseif V_inj_sous_etape(i) == V_elu2 && Sous_Etape(i) == 3
                plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k--','LineWidth',1);
                text(V_inj_total(i)-V_elu2/2,100,'Raffinat','HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','cyan');
                
            elseif V_inj_sous_etape(i) == V_boucle && Sous_Etape(i) == 4 
                if Etape(i) < Ncol_tot
                    plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k--','LineWidth',2);
                else
                    plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k-','LineWidth',2);
                end
                text(V_inj_total(i)-V_boucle/2,100,'Boucle','HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','magenta');
            end
        
        elseif Ncycles*Ncol_tot <= 12
             if V_inj_sous_etape(i) == V_boucle && Sous_Etape(i) == 4 
                if Etape(i) < Ncol_tot
                    plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k--','LineWidth',2);
                else
                    plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k-','LineWidth',2);
                end
                Etape_fig2 = num2str(Etape(i),'%i');
                text(V_inj_total(i)-V_inj_etape(i)/2,100,Etape_fig2,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','black');
             end
             
        else
            if V_inj_sous_etape(i) == V_boucle && Sous_Etape(i) == 4 && Etape(i) == Ncol_tot
                plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k-','LineWidth',2);
                Cycle_fig2 = num2str(Cycle(i),'%i');
                text(V_inj_total(i)-V_inj_cycle(i)/2,100,Cycle_fig2,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','black');
            end
        end
    end
    hold off    
    
    Titre = 'Pureté des composés dans l extrait';
    title(Titre);
    xlabel('Volume injecté (BV)');
    ylabel('Pureté (%)');
    xlim([0 max(V_inj_total(1:i_fin))]);
    ylim([0 100]);

    if Ncycles*Ncol_tot <= 12
        hleg5 = legend('A1 global','A2 global','A3 global','SO4 global','A1 cycle','A2 cycle','A3 cycle','SO4 cycle','A1 etape','A2 etape','A3 etape','SO4 etape','Location','East');
    else    
        hleg5 = legend('A1 global','A2 global','A3 global','SO4 global','A1 cycle','A2 cycle','A3 cycle','SO4 cycle','Location','East');
    end
    set(hleg5,'FontSize',9);
    
    % Pureté des acides organiques dans le raffinat en fonction du volume total injecté
    subplot(2,1,2)
    plot(V_inj_total(1:i_fin),Purete_A1tot_raff_total(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[1 0 0]);
    hold on
    plot(V_inj_total(1:i_fin),Purete_A2tot_raff_total(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 1 0]);
    plot(V_inj_total(1:i_fin),Purete_A3tot_raff_total(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0 1]);
    plot(V_inj_total(1:i_fin),Purete_SO4tot_raff_total(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 1 1]);
    plot(V_inj_total(1:i_fin),Purete_A1tot_raff_cycle(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0.7 0 0]);
    plot(V_inj_total(1:i_fin),Purete_A2tot_raff_cycle(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0.7 0]);
    plot(V_inj_total(1:i_fin),Purete_A3tot_raff_cycle(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0 0.7]);
    plot(V_inj_total(1:i_fin),Purete_SO4tot_raff_cycle(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0.7 0.7]);
    
    if Ncycles*Ncol_tot <= 12
        plot(V_inj_total(1:i_fin),Purete_A1tot_raff_etape(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0.4 0 0]);
        plot(V_inj_total(1:i_fin),Purete_A2tot_raff_etape(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0.4 0]);
        plot(V_inj_total(1:i_fin),Purete_A3tot_raff_etape(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0 0.4]);
        plot(V_inj_total(1:i_fin),Purete_SO4tot_raff_etape(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0.4 0.4]);
    end
    
    for i=2:i_fin
        V_ech = F_V_ech(Cycle(i))*V_ech_est;
        V_elu1 = F_V_elu1(Cycle(i))*V_elu1_est;
        V_elu2 = F_V_elu2(Cycle(i))*V_elu2_est;
        V_boucle = F_V_boucle(Cycle(i))*V_boucle_est;
        
        if Ncycles == 1
            if V_inj_sous_etape(i) == V_ech && Sous_Etape(i) == 1
                plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k--','LineWidth',1);
                text(V_inj_total(i)-V_ech/2,100,'Feed','HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','red');

            elseif V_inj_sous_etape(i) == V_elu1 && Sous_Etape(i) == 2 
                plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k--','LineWidth',1);
                text(V_inj_total(i)-V_elu1/2,100,'Extrait','HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','green');

            elseif V_inj_sous_etape(i) == V_elu2 && Sous_Etape(i) == 3
                plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k--','LineWidth',1);
                text(V_inj_total(i)-V_elu2/2,100,'Raffinat','HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','cyan');
                
            elseif V_inj_sous_etape(i) == V_boucle && Sous_Etape(i) == 4 
                if Etape(i) < Ncol_tot
                    plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k--','LineWidth',2);
                else
                    plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k-','LineWidth',2);
                end
                text(V_inj_total(i)-V_boucle/2,100,'Boucle','HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','magenta');
            end
        
        elseif Ncycles*Ncol_tot <= 12
             if V_inj_sous_etape(i) == V_boucle && Sous_Etape(i) == 4 
                if Etape(i) < Ncol_tot
                    plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k--','LineWidth',2);
                else
                    plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k-','LineWidth',2);
                end
                Etape_fig2 = num2str(Etape(i),'%i');
                text(V_inj_total(i)-V_inj_etape(i)/2,100,Etape_fig2,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','black');
             end
             
        else
            if V_inj_sous_etape(i) == V_boucle && Sous_Etape(i) == 4 && Etape(i) == Ncol_tot
                plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k-','LineWidth',2);
                Cycle_fig2 = num2str(Cycle(i),'%i');
                text(V_inj_total(i)-V_inj_cycle(i)/2,100,Cycle_fig2,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','black');
            end
        end
    end
    hold off    
    
    Titre = 'Pureté des composés dans le raffinat';
    title(Titre);
    xlabel('Volume injecté (BV)');
    ylabel('Pureté (%)');
    xlim([0 max(V_inj_total(1:i_fin))]);
    ylim([0 100]);

    if Ncycles*Ncol_tot <= 12
        hleg6 = legend('A1 global','A2 global','A3 global','SO4 global','A1 cycle','A2 cycle','A3 cycle','SO4 cycle','A1 etape','A2 etape','A3 etape','SO4 etape','Location','East');
    else    
        hleg6 = legend('A1 global','A2 global','A3 global','SO4 global','A1 cycle','A2 cycle','A3 cycle','SO4 cycle','Location','East');
    end
    set(hleg6,'FontSize',9);

figure(3)

    % Taux de récupération des composés dans l'extrait en fonction du volume total injecté
    subplot(2,1,1)
    plot(V_inj_total(1:i_fin),Recuperation_A1tot_ext_total(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[1 0 0]);
    hold on
    plot(V_inj_total(1:i_fin),Recuperation_A2tot_ext_total(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 1 0]);
    plot(V_inj_total(1:i_fin),Recuperation_A3tot_ext_total(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0 1]);
    plot(V_inj_total(1:i_fin),Recuperation_SO4tot_ext_total(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 1 1]);
    plot(V_inj_total(1:i_fin),Recuperation_A1tot_ext_cycle(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0.7 0 0]);
    plot(V_inj_total(1:i_fin),Recuperation_A2tot_ext_cycle(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0.7 0]);
    plot(V_inj_total(1:i_fin),Recuperation_A3tot_ext_cycle(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0 0.7]);
    plot(V_inj_total(1:i_fin),Recuperation_SO4tot_ext_cycle(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0.7 0.7]);
    
    if Ncycles*Ncol_tot <= 12
        plot(V_inj_total(1:i_fin),Recuperation_A1tot_ext_etape(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0.4 0 0]);
        plot(V_inj_total(1:i_fin),Recuperation_A2tot_ext_etape(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0.4 0]);
        plot(V_inj_total(1:i_fin),Recuperation_A3tot_ext_etape(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0 0.4]);
        plot(V_inj_total(1:i_fin),Recuperation_SO4tot_ext_etape(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0.4 0.4]);
    end
    
    for i=2:i_fin
        V_ech = F_V_ech(Cycle(i))*V_ech_est;
        V_elu1 = F_V_elu1(Cycle(i))*V_elu1_est;
        V_elu2 = F_V_elu2(Cycle(i))*V_elu2_est;
        V_boucle = F_V_boucle(Cycle(i))*V_boucle_est;
        
        if Ncycles == 1
            if V_inj_sous_etape(i) == V_ech && Sous_Etape(i) == 1
                plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k--','LineWidth',1);
                text(V_inj_total(i)-V_ech/2,100,'Feed','HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','red');

            elseif V_inj_sous_etape(i) == V_elu1 && Sous_Etape(i) == 2 
                plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k--','LineWidth',1);
                text(V_inj_total(i)-V_elu1/2,100,'Extrait','HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','green');

            elseif V_inj_sous_etape(i) == V_elu2 && Sous_Etape(i) == 3
                plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k--','LineWidth',1);
                text(V_inj_total(i)-V_elu2/2,100,'Raffinat','HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','cyan');
                
            elseif V_inj_sous_etape(i) == V_boucle && Sous_Etape(i) == 4
                if Etape(i) < Ncol_tot
                    plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k--','LineWidth',2);
                else
                    plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k-','LineWidth',2);
                end
                text(V_inj_total(i)-V_boucle/2,100,'Boucle','HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','magenta');
            end
        
        elseif Ncycles*Ncol_tot <= 12
             if V_inj_sous_etape(i) == V_boucle && Sous_Etape(i) == 4 
                if Etape(i) < Ncol_tot
                    plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k--','LineWidth',2);
                else
                    plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k-','LineWidth',2);
                end
                Etape_fig2 = num2str(Etape(i),'%i');
                text(V_inj_total(i)-V_inj_etape(i)/2,100,Etape_fig2,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','black');
             end
             
        else
            if V_inj_sous_etape(i) == V_boucle && Sous_Etape(i) == 4 && Etape(i) == Ncol_tot
                plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k-','LineWidth',2);
                Cycle_fig2 = num2str(Cycle(i),'%i');
                text(V_inj_total(i)-V_inj_cycle(i)/2,100,Cycle_fig2,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','black');
            end
        end
    end
    hold off    
    
    Titre = 'Taux de récupération des composés dans l extrait';
    title(Titre);
    xlabel('Volume injecté (BV)');
    ylabel('Taux de récupération (%)');
    xlim([0 max(V_inj_total(1:i_fin))]);
    ylim([0 100])
    
    if Ncycles*Ncol_tot <= 12
        hleg8 = legend('A1 global','A2 global','A3 global','SO4 global','A1 cycle','A2 cycle','A3 cycle','SO4 cycle','A1 etape','A2 etape','A3 etape','SO4 etape','Location','East');
    else    
        hleg8 = legend('A1 global','A2 global','A3 global','SO4 global','A1 cycle','A2 cycle','A3 cycle','SO4 cycle','Location','East');
    end
    set(hleg8,'FontSize',9);
    
    % Taux de récupération des composés dans le raffinat en fonction du volume total injecté
    subplot(2,1,2)
    plot(V_inj_total(1:i_fin),Recuperation_A1tot_raff_total(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[1 0 0]);
    hold on
    plot(V_inj_total(1:i_fin),Recuperation_A2tot_raff_total(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 1 0]);
    plot(V_inj_total(1:i_fin),Recuperation_A3tot_raff_total(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0 1]);
    plot(V_inj_total(1:i_fin),Recuperation_SO4tot_raff_total(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 1 1]);
    plot(V_inj_total(1:i_fin),Recuperation_A1tot_raff_cycle(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0.7 0 0]);
    plot(V_inj_total(1:i_fin),Recuperation_A2tot_raff_cycle(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0.7 0]);
    plot(V_inj_total(1:i_fin),Recuperation_A3tot_raff_cycle(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0 0.7]);
    plot(V_inj_total(1:i_fin),Recuperation_SO4tot_raff_cycle(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0.7 0.7]);
    
    if Ncycles*Ncol_tot <= 12
        plot(V_inj_total(1:i_fin),Recuperation_A1tot_raff_etape(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0.4 0 0]);
        plot(V_inj_total(1:i_fin),Recuperation_A2tot_raff_etape(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0.4 0]);
        plot(V_inj_total(1:i_fin),Recuperation_A3tot_raff_etape(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0 0.4]);
        plot(V_inj_total(1:i_fin),Recuperation_SO4tot_raff_etape(:),'LineStyle','-','LineWidth',1,'Marker','o','MarkerSize',2,'Color',[0 0.4 0.4]);
    end
    
    for i=2:i_fin
        V_ech = F_V_ech(Cycle(i))*V_ech_est;
        V_elu1 = F_V_elu1(Cycle(i))*V_elu1_est;
        V_elu2 = F_V_elu2(Cycle(i))*V_elu2_est;
        V_boucle = F_V_boucle(Cycle(i))*V_boucle_est;
        
        if Ncycles == 1
            if V_inj_sous_etape(i) == V_ech && Sous_Etape(i) == 1
                plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k--','LineWidth',1);
                text(V_inj_total(i)-V_ech/2,100,'Feed','HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','red');

            elseif V_inj_sous_etape(i) == V_elu1 && Sous_Etape(i) == 2 
                plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k--','LineWidth',1);
                text(V_inj_total(i)-V_elu1/2,100,'Extrait','HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','green');

            elseif V_inj_sous_etape(i) == V_elu2 && Sous_Etape(i) == 3
                plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k--','LineWidth',1);
                text(V_inj_total(i)-V_elu2/2,100,'Raffinat','HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','cyan');

            elseif V_inj_sous_etape(i) == V_boucle && Sous_Etape(i) == 4 
                if Etape(i) < Ncol_tot
                    plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k--','LineWidth',2);
                else
                    plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k-','LineWidth',2);
                end
                text(V_inj_total(i)-V_boucle/2,100,'Boucle','HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','magenta');
            end
        
        elseif Ncycles*Ncol_tot <= 12
             if V_inj_sous_etape(i) == V_boucle && Sous_Etape(i) == 4 
                if Etape(i) < Ncol_tot
                    plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k--','LineWidth',2);
                else
                    plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k-','LineWidth',2);
                end
                Etape_fig2 = num2str(Etape(i),'%i');
                text(V_inj_total(i)-V_inj_etape(i)/2,100,Etape_fig2,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','black');
             end
             
        else
            if V_inj_sous_etape(i) == V_boucle && Sous_Etape(i) == 4 && Etape(i) == Ncol_tot
                plot(V_inj_total(i)*ones(3,1),linspace(0,100,3),'k-','LineWidth',2);
                Cycle_fig2 = num2str(Cycle(i),'%i');
                text(V_inj_total(i)-V_inj_cycle(i)/2,100,Cycle_fig2,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','black');
            end
        end
    end
    hold off    
    
    Titre = 'Taux de récupération des composés dans le raffinat';
    title(Titre);
    xlabel('Volume injecté (BV)');
    ylabel('Taux de récupération (%)');
    xlim([0 max(V_inj_total(1:i_fin))]);
    ylim([0 100])
    
    if Ncycles*Ncol_tot <= 12
        hleg9 = legend('A1 global','A2 global','A3 global','SO4 global','A1 cycle','A2 cycle','A3 cycle','SO4 cycle','A1 etape','A2 etape','A3 etape','SO4 etape','Location','East');
    else    
        hleg9 = legend('A1 global','A2 global','A3 global','SO4 global','A1 cycle','A2 cycle','A3 cycle','SO4 cycle','Location','East');
    end
    set(hleg9,'FontSize',9);
    
% Création graphique : évolution des conditions opératoires de la SSMB
figure(4)

    % Volumes injectés lors des différentes sous-étapes
    V_ech_fig3 = zeros(Ncol_tot*Ncycles,1);
    V_elu1_fig3 = zeros(Ncol_tot*Ncycles,1);
    V_elu2_fig3 = zeros(Ncol_tot*Ncycles,1);
    V_boucle_fig3 = zeros(Ncol_tot*Ncycles,1);

    for i=1:Ncycles
        for j=1:Ncol_tot
            V_ech_fig3(Ncol_tot*(i-1)+j) = F_V_ech(i)*V_ech_est;
            V_elu1_fig3(Ncol_tot*(i-1)+j) = F_V_elu1(i)*V_elu1_est;
            V_elu2_fig3(Ncol_tot*(i-1)+j) = F_V_elu2(i)*V_elu2_est;
            V_boucle_fig3(Ncol_tot*(i-1)+j) = F_V_boucle(i)*V_boucle_est;
        end
    end

    subplot(2,1,1)
    plot(1:Ncol_tot*Ncycles,V_ech_fig3(:),'LineStyle','--','LineWidth',2,'Marker','o','MarkerSize',4,'Color','red');
    hold on
    plot(1:Ncol_tot*Ncycles,V_elu1_fig3(:),'LineStyle','--','LineWidth',2,'Marker','o','MarkerSize',4,'Color','green');
    plot(1:Ncol_tot*Ncycles,V_elu2_fig3(:),'LineStyle','--','LineWidth',2,'Marker','o','MarkerSize',4,'Color','cyan');
    plot(1:Ncol_tot*Ncycles,V_boucle_fig3(:),'LineStyle','--','LineWidth',2,'Marker','o','MarkerSize',4,'Color','magenta');
    
    Titre = 'Volumes injectés lors des différentes sous-étapes';
    title(Titre);
    xlabel('Etape');
    ylabel('Volume injecté (BV)');
    legend('1:Feed+Raff','2:Elu+Ext','3:Elu+Raff','4:Boucle');
    
    V_max_fig3 = max(max([V_ech_fig3(:);V_elu1_fig3(:);V_elu2_fig3(:);V_boucle_fig3(:)]));
    xlim([0.5 (Ncol_tot*Ncycles+0.5)]);
    ylim([0 V_max_fig3])
    
    plot(0.5*ones(3,1),linspace(0,V_max_fig3,3),'k-','LineWidth',2);
    for i=1:Ncycles
        plot((Ncol_tot*i+0.5)*ones(3,1),linspace(0,V_max_fig3,3),'k-','LineWidth',2);
        Cycle_fig3 = num2str(i,'%i');
        text(Ncol_tot*i-1.5,V_max_fig3,Cycle_fig3,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','black');
    end
    hold off
    
    % Volumes injectés dans les différentes zones du SSMB
    BV1_fig3 = V_elu1_fig3+V_elu2_fig3+V_boucle_fig3+Volume_tuyau;
    BV2_fig3 = V_elu2_fig3+V_boucle_fig3+Volume_tuyau;
    BV3_fig3 = V_ech_fig3+V_elu2_fig3+V_boucle_fig3+Volume_tuyau;
    BV4_fig3 = V_boucle_fig3+Volume_tuyau;
    
    subplot(2,1,2)
    plot(1:Ncol_tot*Ncycles,BV1_fig3(:),'LineStyle','--','LineWidth',2,'Marker','o','MarkerSize',4,'Color','green');
    hold on
    plot(1:Ncol_tot*Ncycles,BV2_fig3(:),'LineStyle','--','LineWidth',2,'Marker','o','MarkerSize',4,'Color','cyan');
    plot(1:Ncol_tot*Ncycles,BV3_fig3(:),'LineStyle','--','LineWidth',2,'Marker','o','MarkerSize',4,'Color','red');
    plot(1:Ncol_tot*Ncycles,BV4_fig3(:),'LineStyle','--','LineWidth',2,'Marker','o','MarkerSize',4,'Color','blue');
    
    Titre = 'Volumes injectés dans les différentes zones du SSMB';
    title(Titre);
    xlabel('Etape');
    ylabel('Volume injecté (BV)');
    legend('BV1','BV2','BV3','BV4');    
    
    BV_max_fig3 = max(max([BV1_fig3(:);BV2_fig3(:);BV3_fig3(:);BV4_fig3(:)]));
    xlim([0.5 (Ncol_tot*Ncycles+0.5)]);
    ylim([0 BV_max_fig3])
    
    plot(0.5*ones(3,1),linspace(0,BV_max_fig3,3),'k-','LineWidth',2);
    for i=1:Ncycles
        plot((Ncol_tot*i+0.5)*ones(3,1),linspace(0,BV_max_fig3,3),'k-','LineWidth',2);
        Cycle_fig3 = num2str(i,'%i');
        text(Ncol_tot*i-1.5,BV_max_fig3,Cycle_fig3,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9,'Color','black');
    end
    hold off

% Création graphique : Visualisation de la variation des concentrations en un point donné (sortie de la dernière colonne) pendant un cycle
figure(5)

    i = i_fin;
    while Cycle(i) == Ncycles
        i = i-1;                
    end
    i_deb = i+1;

    SO4tot_mob_fig4 = SO4tot_mob(i_deb:i_fin,end);
    A1tot_mob_fig4 = A1tot_mob(i_deb:i_fin,end);
    A2tot_mob_fig4 = A2tot_mob(i_deb:i_fin,end);
    A3tot_mob_fig4 = A3tot_mob(i_deb:i_fin,end);

    n = 0;
    V_inj_fig4 = zeros(i_fin-i_deb+1,1);
    for i=i_deb:i_fin
        n = n+1;
        if (Etape(i)==1 && Sous_Etape(i)==4) || (Etape(i)==2 && (Sous_Etape(i)==1 || Sous_Etape(i)==3 || Sous_Etape(i)==4)) || (Etape(i)==3 && (Sous_Etape(i)==3 || Sous_Etape(i)==4)) || (Etape(i)==4 && (Sous_Etape(i)==2 || Sous_Etape(i)==3 || Sous_Etape(i)==4))
            if i>i_deb
                V_inj_fig4(n) = V_inj_fig4(n-1) + V_inj_total(i) - V_inj_total(i-1);
            else
                V_inj_fig4(n) = V_inj_total(i) - V_inj_total(i-1);
            end
        else
            if i>i_deb
                V_inj_fig4(n) = V_inj_fig4(n-1);
            else
                V_inj_fig4(n) = 0;
            end
        end
    end

    plot(V_inj_fig4(:),SO4tot_mob_fig4(:),'LineStyle','-','LineWidth',2,'Marker','o','MarkerSize',6,'Color','blue');
    hold on
    plot(V_inj_fig4(:),A1tot_mob_fig4(:),'LineStyle','-','LineWidth',2,'Marker','o','MarkerSize',6,'Color','cyan');
    plot(V_inj_fig4(:),A2tot_mob_fig4(:),'LineStyle','-','LineWidth',2,'Marker','o','MarkerSize',6,'Color','red');
    plot(V_inj_fig4(:),A3tot_mob_fig4(:),'LineStyle','-','LineWidth',2,'Marker','o','MarkerSize',6,'Color','green');
    
    if load_exp == true
        Exp_resultats_structure = load(nom_fichier_exp);
        Exp_resultats = Exp_resultats_structure.Exp_resultats;
        Exp_volume_inj_BV(:,1) = Exp_resultats(:,1);
        Exp_A1tot_mob(:,1)= Exp_resultats(:,2);
        Exp_A2tot_mob(:,1)= Exp_resultats(:,3);
        Exp_A3tot_mob(:,1)= Exp_resultats(:,4);
        
        plot(Exp_V_inj_BV(:),Exp_Cond_mob(:),'LineStyle','--','LineWidth',2,'Marker','*','MarkerSize',6,'Color','blue');
        plot(Exp_V_inj_BV(:),Exp_A1tot_mob(:),'LineStyle','--','LineWidth',2,'Marker','*','MarkerSize',6,'Color','cyan');
        plot(Exp_V_inj_BV(:),Exp_A2tot_mob(:),'LineStyle','--','LineWidth',2,'Marker','*','MarkerSize',6,'Color','red');
        plot(Exp_V_inj_BV(:),Exp_A3tot_mob(:),'LineStyle','--','LineWidth',2,'Marker','*','MarkerSize',6,'Color','green');
        
        legend('SO4 sim','A1 sim','A2 sim','A3 sim','Cond exp','A1 exp','A2 exp','A3 exp','Location','East');    
        C_max_fig4 = max(max([max([SO4tot_mob_fig4(:);A1tot_mob_fig4(:);A2tot_mob_fig4(:);A3tot_mob_fig4(:)]);max([Exp_Cond_mob(:);Exp_A1tot_mob(:);Exp_A2tot_mob(:);Exp_A3tot_mob(:)])]));
    
    else
        legend('SO4','A1','A2','A3','Location','East');    
        C_max_fig4 = max(max([SO4tot_mob_fig4(:);A1tot_mob_fig4(:);A2tot_mob_fig4(:);A3tot_mob_fig4(:)]));
    end

    Titre = 'Evolution des concentrations en sortie de colonne pendant un cycle';
    title(Titre);
    xlabel('Volume déplacé (BV)');
    ylabel('Concentration (mol/L)');
    xlim([0 max(V_inj_fig4(:))]);
    ylim([0 1.05*C_max_fig4]);

    plot(0*ones(3,1),linspace(0,1.05*C_max_fig4,3),'k-','LineWidth',5);
    text(V_boucle_fig3(end)/2,1.05*C_max_fig4,'BV4','HorizontalAlignment','center','VerticalAlignment','top','FontSize',9);
    text(V_boucle_fig3(end)/2,C_max_fig4,'Vboucle','HorizontalAlignment','center','VerticalAlignment','top','FontSize',7);
    plot(V_boucle_fig3(end)*ones(3,1),linspace(0,1.05*C_max_fig4,3),'g-','LineWidth',3);
    text(V_boucle_fig3(end)+V_ech_fig3(end)/2,C_max_fig4,'Vech','HorizontalAlignment','center','VerticalAlignment','top','FontSize',7);
    plot((V_boucle_fig3(end)+V_ech_fig3(end))*ones(3,1),linspace(0,C_max_fig4,3),'g--','LineWidth',2);
    text(V_boucle_fig3(end)+(V_ech_fig3(end)+V_elu2_fig3(end)+V_boucle_fig3(end))/2,1.05*C_max_fig4,'BV3','HorizontalAlignment','center','VerticalAlignment','top','FontSize',9);
    text(V_boucle_fig3(end)+V_ech_fig3(end)+V_elu2_fig3(end)/2,C_max_fig4,'Velu2','HorizontalAlignment','center','VerticalAlignment','top','FontSize',7);
    plot((V_boucle_fig3(end)+V_ech_fig3(end)+V_elu2_fig3(end))*ones(3,1),linspace(0,C_max_fig4,3),'g--','LineWidth',2);
    text(V_boucle_fig3(end)+V_ech_fig3(end)+V_elu2_fig3(end)+V_boucle_fig3(end)/2,C_max_fig4,'Vboucle','HorizontalAlignment','center','VerticalAlignment','top','FontSize',7);
    plot((V_boucle_fig3(end)+V_ech_fig3(end)+V_elu2_fig3(end)+V_boucle_fig3(end))*ones(3,1),linspace(0,1.05*C_max_fig4,3),'k-','LineWidth',3);
    text(V_boucle_fig3(end)+BV3_fig3(end)-Volume_tuyau+V_elu2_fig3(end)/2,C_max_fig4,'Velu2','HorizontalAlignment','center','VerticalAlignment','top','FontSize',7);
    plot((V_boucle_fig3(end)+BV3_fig3(end)-Volume_tuyau+V_elu2_fig3(end))*ones(3,1),linspace(0,C_max_fig4,3),'k--','LineWidth',2);
    text(V_boucle_fig3(end)+BV3_fig3(end)-Volume_tuyau+(V_elu2_fig3(end)+V_boucle_fig3(end))/2,1.05*C_max_fig4,'BV2','HorizontalAlignment','center','VerticalAlignment','top','FontSize',9);
    text(V_boucle_fig3(end)+BV3_fig3(end)-Volume_tuyau+V_elu2_fig3(end)+V_boucle_fig3(end)/2,C_max_fig4,'Vboucle','HorizontalAlignment','center','VerticalAlignment','top','FontSize',7);
    plot((V_boucle_fig3(end)+BV3_fig3(end)-Volume_tuyau+BV2_fig3(end)-Volume_tuyau)*ones(3,1),linspace(0,1.05*C_max_fig4,3),'r-','LineWidth',3);
    text(V_boucle_fig3(end)+BV3_fig3(end)-Volume_tuyau+BV2_fig3(end)-Volume_tuyau+V_elu1_fig3(end)/2,C_max_fig4,'Velu1','HorizontalAlignment','center','VerticalAlignment','top','FontSize',7);
    plot((V_boucle_fig3(end)+BV3_fig3(end)-Volume_tuyau+BV2_fig3(end)-Volume_tuyau+V_elu1_fig3(end))*ones(3,1),linspace(0,C_max_fig4,3),'r--','LineWidth',2);
    text(V_boucle_fig3(end)+BV3_fig3(end)-Volume_tuyau+BV2_fig3(end)-Volume_tuyau+(V_elu1_fig3(end)+V_elu2_fig3(end)+V_boucle_fig3(end))/2,1.05*C_max_fig4,'BV1','HorizontalAlignment','center','VerticalAlignment','top','FontSize',9);
    text(V_boucle_fig3(end)+BV3_fig3(end)-Volume_tuyau+BV2_fig3(end)-Volume_tuyau+V_elu1_fig3(end)+V_elu2_fig3(end)/2,C_max_fig4,'Velu2','HorizontalAlignment','center','VerticalAlignment','top','FontSize',7);
    plot((V_boucle_fig3(end)+BV3_fig3(end)-Volume_tuyau+BV2_fig3(end)-Volume_tuyau+V_elu1_fig3(end)+V_elu2_fig3(end))*ones(3,1),linspace(0,C_max_fig4,3),'k--','LineWidth',2);
    text(V_boucle_fig3(end)+BV3_fig3(end)-Volume_tuyau+BV2_fig3(end)-Volume_tuyau+V_elu1_fig3(end)+V_elu2_fig3(end)+V_boucle_fig3(end)/2,C_max_fig4,'Vboucle','HorizontalAlignment','center','VerticalAlignment','top','FontSize',7);
    plot((V_boucle_fig3(end)+BV3_fig3(end)-Volume_tuyau+BV2_fig3(end)-Volume_tuyau+V_elu1_fig3(end)+V_elu2_fig3(end)+V_boucle_fig3(end))*ones(3,1),linspace(0,1.05*C_max_fig4,3),'k-','LineWidth',4);    
    
% Sauvegarde de toutes les données
save (nom_fichier)

% Sauvegarde de la composition finale des colonnes
H_mob_load = H_mob_new;
HSO4_mob_load = HSO4_mob_new;
SO4_mob_load = SO4_mob_new;
SO4tot_mob_load = SO4tot_mob_new;
A1H_mob_load = A1H_mob_new;
A2H_mob_load = A2H_mob_new;
A3H_mob_load = A3H_mob_new;
A1_mob_load = A1_mob_new;
A2_mob_load = A2_mob_new;
A3_mob_load = A3_mob_new;
A1tot_mob_load = A1tot_mob_new;
A2tot_mob_load = A2tot_mob_new;
A3tot_mob_load = A3tot_mob_new;

HSO4_stat_load = HSO4_stat_new;
SO4_stat_load = SO4_stat_new;
SO4tot_stat_load = SO4tot_stat_new;
A1H_HSO4_stat_load = A1H_HSO4_stat_new;
A2H_HSO4_stat_load = A2H_HSO4_stat_new;
A3H_HSO4_stat_load = A3H_HSO4_stat_new;
A1H_SO4_stat_load = A1H_SO4_stat_new;
A2H_SO4_stat_load = A2H_SO4_stat_new;
A3H_SO4_stat_load = A3H_SO4_stat_new;
A1H_mat_stat_load = A1H_mat_stat_new;
A2H_mat_stat_load = A2H_mat_stat_new;
A3H_mat_stat_load = A3H_mat_stat_new;
A1_stat_load = A1_stat_new;
A2_stat_load = A2_stat_new;
A3_stat_load = A3_stat_new;
A1tot_stat_load = A1tot_stat_new;
A2tot_stat_load = A2tot_stat_new;
A3tot_stat_load = A3tot_stat_new;

F_ionique_load = F_ionique_new;
Activite_AH_load = Activite_AH_new;
Activite_HSO4_load = Activite_HSO4_new;

save(nom_fichier_CI_new,'H_mob_load','HSO4_mob_load','SO4_mob_load','SO4tot_mob_load','A1H_mob_load','A2H_mob_load','A3H_mob_load','A1_mob_load','A2_mob_load','A3_mob_load','A1tot_mob_load','A2tot_mob_load','A3tot_mob_load','HSO4_stat_load','SO4_stat_load','SO4tot_stat_load','A1H_HSO4_stat_load','A2H_HSO4_stat_load','A3H_HSO4_stat_load','A1H_SO4_stat_load','A2H_SO4_stat_load','A3H_SO4_stat_load','A1H_mat_stat_load','A2H_mat_stat_load','A3H_mat_stat_load','A1_stat_load','A2_stat_load','A3_stat_load','A1tot_stat_load','A2tot_stat_load','A3tot_stat_load','F_ionique_load','Activite_AH_load','Activite_HSO4_load','BV1_fig3','BV2_fig3','BV3_fig3','BV4_fig3')
