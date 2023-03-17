
% This is a version that has time dynamic ATPi
% Incorporate dstates are from energy metabolites and force generation model
% A96 to A101 and A124 to A138

% put all derivatives first

% experimenting with parameter values to recover ATPi and CrPi recovery
% using parameter values from cellml if there is discrepency

% Summary of changes - all parameters using cellml values
% added outputting normalised force FN_Ca


function [dCam, dADPm, dDpsi, dSCoA, dNADH, dISOC, dAKG, dSucc, dFUM, dMAL, dOaa, dASP, dATPi, dATPi_cyto, dCrPi_mito, dCrPi_cyto, dP0, dP1, dP2, dP3, dN1, dLTRPNCa, dHTRPNCa, force_norm, force, VNO] ...
    = calculateF_mitene_full_v6_cellml_VnaCa_unscaled(Cam, ADPm, Dpsi, SCoA, NADH, ISOC, AKG, Succ, FUM, MAL, Oaa, ASP, ATPi, ATPi_cyto, CrPi_mito, CrPi_cyto, P0, P1, P2, P3, N1, LTRPNCa, HTRPNCa, Cai, Nai, Jup, INaK, IpCa)

    % Elements needed from TorOrd: Cai, Nai, Jup, INaK, IpCa

    % Base conditions
    Cmito = 1.812e-3; % mM/mV Inner membrane capacitance
    inv_Cmito = 1/Cmito; % inverse
    fm = 3.0e-4; % fraction of free Cam
    b = 0.5; % delPsi dependence on NCX 
    kact = 3.8e-4; % mM activation constant
    F = 96.5; % C/mmol Faraday constant
    T = 310; % K absolute temp 
    R = 8.314; % J.mol-1.K-1 Universal gas constant
    
    % Cam Vuni
    Vmuni = 0.0275; % mM.ms-1 Vmax uniporter Ca transport
    ktrans = 0.019; % Kd for translocated Ca
    Vmuni_ktrans = Vmuni / ktrans;

    RT_over_F = R*T/F;
    F_over_RT = 1.0 / RT_over_F;
    FRT2 = 2.0 * F_over_RT;
    FRT2_Dpsi = FRT2 * ( (Dpsi) - 91.0 ); %%% 91 not 81

    inv_ktrans = 1.0 / ktrans;
    Cai_ktrans_plus1 = 1.0 + (Cai) * inv_ktrans;
    Cai_ktrans_plus1_p3 = Cai_ktrans_plus1 * Cai_ktrans_plus1 * Cai_ktrans_plus1;

    inv_kact = 1.0 / kact;
    na = 2.8; % dimensionless Uniporter activation cooperativity
    L = 110.0; % dimensionless? Keq for conformational transitions in uniporter

    Vuni =  Vmuni_ktrans * (Cai) * FRT2_Dpsi * Cai_ktrans_plus1_p3  / ...
        ( ( Cai_ktrans_plus1_p3 * Cai_ktrans_plus1 + L / power( 1.0 + (Cai) * inv_kact, na) ) * ( 1.0 - exp(-FRT2_Dpsi) ) );
    
    % Cam VnaCa
    b_05 = b * 0.5;
    Kna = 9.4; % mM Antiporter Na constant 9.4
    n = 3; % dimensionless, NCX antiporter cooperativity
    Knca = 3.75e-4; % mM antiporter Ca constant
    VmNC = 0.625E-4; % 1e-4 mM.ms-1 vmax of NaCa antiporter, 0.625E-4 in cellml %%%%% EXPERIMENTING - cellml parameter doesn't do much

    VnaCa = VmNC * exp( b_05 * FRT2_Dpsi ) * (Cam) / ( (Cai) * power( ( 1.0 + Kna / (Nai) ) , n) * ( 1.0 + Knca / (Cam) ) );

    % ADPm VANT
    VmDT = 0.015; % mM.ms-1 ANT max rate? 0.025 in paper, 0.015 in cellml
	VmDT_75 = 0.75 * VmDT;
	VmDT_20 = 20.0 * VmDT;
    ADP = 8.0 - (ATPi);
    ATPi_ADP = (ATPi) / ADP;
    Cm = 1.5; % mM as from Cellml - C_A in paper
    ATPm = Cm - (ADPm);
    ADPm_ATPm = (ADPm) / ATPm;
    tenDiv9 = 10.0 / 9.0;
    hm = 0.5; % dimensionless, fraction of delta psi m
    hm_F_over_RT = hm * F_over_RT;

    VANT = ( VmDT_75 - VmDT_20 * ATPi_ADP * ADPm_ATPm * exp(-F_over_RT * (Dpsi)) ) /...
		( ( 1.0 + tenDiv9 * ATPi_ADP * exp( -hm_F_over_RT * (Dpsi) ) ) * (1.0 + 18 * ADPm_ATPm ) );

    % ADPm VATPase
    p1 = 1.346e-8; % dimensionless, sum of products of rate constants
    pa = 1.656e-8; % ms-1 sum of products of rate constants
    pc2 = 4.585e-17; % ms-1 sum of products of rate constants
    pc1 = 9.651e-17; % ms-1 sum of products of rate constants
    Dpsio = 50; % mV phase boundary potential
    exp_3_FRT_Dpsio = exp( 3.0 * Dpsio * F_over_RT); 
    VATPase_C1 = ( 100.0 * pa + pc1 * exp_3_FRT_Dpsio);

    FRT_3 = 3.0 * F_over_RT;
    DpH = -0.6; % pH gradient across inner membrane
    DmuH_Constant = -2.303 * RT_over_F * DpH;
    DmuH = DmuH_Constant + (Dpsi);
    exp_3FRT_DmuH = exp(FRT_3 * DmuH);

    rhoF1 = 0.05; % mM Conc of F1F0-ATPase 1.5 in cortassa, 0.05 in cellml
    p1_exp_3_FRT_Dpsio = p1 * exp_3_FRT_Dpsio;
    kf1 = 1.71e6; % dimensionless equilibrium constant of ATP hydrolysis
    Pi = 2.0; % mM Inorganic phosphate conc
    kf1_Pi = kf1 / Pi;

    AF1 =  kf1_Pi * ATPm / (ADPm);
    p2 = 7.739e-7; % dimensionless sum of products of rate constants 7.739e-4 in cortassa, 7.739e-7 in cellml
    p3 = 6.65e-15; % dimensionless sum of products of rate constants
    denominator_VATPase_Vhu = -rhoF1 / ( exp_3_FRT_Dpsio + p1_exp_3_FRT_Dpsio * AF1 + ( p2 + p3 * AF1 ) * exp_3FRT_DmuH); 
    VATPase = ( (VATPase_C1 + pc2 * exp_3FRT_DmuH) * AF1 - pa * exp_3FRT_DmuH ) * denominator_VATPase_Vhu;

    % ADPm VSL
    kfSL = 0.005; %mM-1.ms-1 forward rate constant of SL 5e-4 in cortassa, 0.005 in cellml
    CoA = 0.02; % mM Coenzyme A conc
    KSLeq = 3.115; % dimensionless, equilibrium constant of SL reaction
    CoA_KSLeq = CoA / KSLeq;
    VSL = kfSL * ( (SCoA) * (ADPm) - CoA_KSLeq * (Succ) * ATPm);

    % DPsi VHNe
    rhoREN = 1.01e-1; % mM Conc of electron carriers 1.01e-1 in cellml
    ra = 6.394e-13; % ms-1 sum of products of rate constants
    rb = 1.762e-16; % ms-1 sum of products of rate constants
    rhoRen_6_ra = 6.0 * rhoREN * ra;

    kres = 1.35e18; % dimensionless equilibrium constant of respiration
    KmIDNAD = 0.923; % mM Michaelis constant for NAD+
    CPN = 10; % mM Total sum of mito pyridine nucleotides
    NAD = CPN - (NADH);
    KmIDNAD_NAD = KmIDNAD / NAD;
    kres_sq_KmIDNAD = kres * kres / KmIDNAD;
    exp_6_FRT_Dpsio = exp( 6.0 * Dpsio * F_over_RT);
    r1 = 2.077e-18; % dimensionless, sum of products of rate constants
    r1_exp_6_FRT_Dpsio = r1 * exp_6_FRT_Dpsio;
    r2 = 1.728e-9; % dimensionless, sum of products of rate constants
    r3 = 1.059e-26; % dimensionless, sum of products of rate 
    g = 0.85; % dimensionless, correction factor for voltage
    FRT_6_g = 6.0 * g * F_over_RT; 
    exp_FRT_6_g_DmuH = exp(FRT_6_g * DmuH); 
    
	AREN =  sqrt ( (NADH) * kres_sq_KmIDNAD * KmIDNAD_NAD) ;
	denominator_VHNe_VNo = 1.0 / ( (exp_6_FRT_Dpsio + r1_exp_6_FRT_Dpsio *  AREN ) + ( r2 + r3 *  AREN ) * exp_FRT_6_g_DmuH );

    rhoRen_6_ra_rb = 6.0 * rhoREN * ( ra + rb );
	
	VHNe = (rhoRen_6_ra * AREN  - rhoRen_6_ra_rb * exp_FRT_6_g_DmuH ) *	denominator_VHNe_VNo;

    % Dpsi VHFe
    rhoREF = 3.75e-4; % mM Conc of electron carriers (II, III, IV) 
    kresf = 5.765e13; % dimensionless, equilibrium constant of FADH2 oxidation
    FADH2 = 1.24; % mM conc of FADH2 reduced
    FAD = 0.01; % mM conc of FAD+ oxidized
    AREF = RT_over_F * log ( kresf * sqrt ( FADH2 / FAD ) );
    exp_AREF_FRT = exp( AREF * F_over_RT );
    VFO_VHFe_C1 = ( 1.0 + r1 * exp_AREF_FRT) * exp_6_FRT_Dpsio;

    r2_r3_exp_AREF_FRT = r2 + r3 * exp_AREF_FRT;

    ra_exp_AREF_FRT = 4.0 * ra * exp_AREF_FRT; % this is where the X4 is - why? who knows
    ra_rb = 4.0 * (ra + rb);
    exp_FRT_6_g_DmuH = exp(FRT_6_g * DmuH);
    
    denominator_VHFe = rhoREF / (VFO_VHFe_C1 + r2_r3_exp_AREF_FRT * exp_FRT_6_g_DmuH);
    VHFe = ( ra_exp_AREF_FRT - ra_rb*exp_FRT_6_g_DmuH ) * denominator_VHFe;

    % Dpsi Vhu
    pb = 3.373e-10; % ms-1 sum of products of rate constants
    pa_300 = 300.0 * pa;
    pa_pb_3 = 3.0 * (pa + pb);
    Vhu = ( pa_300 + pa_300 * AF1  - pa_pb_3 * exp_3FRT_DmuH ) * denominator_VATPase_Vhu;

    % Dspi VANT established in ADPm

    % Dpsi Vhleak
    gh = 2.0e-7; % mM. ms-1. mV-1 Ionic conductance of inner membrane 1e-8 in cortassa, 2.0e-7 in cellml
    Vhleak = gh * DmuH;

    % Dpsi VnaCa established in Cam

    % Dpsi Vuni established in Cam
    
    % NADH VNO
    rc1 = 2.656e-22; % ms-1 sum of products of rate constants
    rc2 = 8.632e-30; % ms-1 sum of products of rate constants 
    ra_rc1_exp_6_FRT_Dpsio = ra + rc1 * exp_6_FRT_Dpsio;
    rhoREN_ra_rc1_exp_6_FRT_Dpsio = 0.5 * rhoREN * ra_rc1_exp_6_FRT_Dpsio;

    rhoREN_rc2 = 0.5 * rhoREN * rc2;
    rhoREN_ra = 0.5 * rhoREN * ra;

    VNO  = ( (rhoREN_ra_rc1_exp_6_FRT_Dpsio + rhoREN_rc2 * exp_FRT_6_g_DmuH) *  AREN  - rhoREN_ra * exp_FRT_6_g_DmuH ) * denominator_VHNe_VNo;

    % NADH VIDH
	KADP = 0.62; % mM activation constant by ADP
    KaCa = 0.0005; % mM activation constant for Ca2+
    KidhNADH = 0.19; % mM inhibition constant by NADH
    inv_KADP = 1.0 / KADP;
	inv_KaCa = 1.0 / KaCa;
    inv_KidhNADH = 1.0 / KidhNADH;
    Fa = 1.0 / (( 1.0 + (ADPm) * inv_KADP) * (1.0 + (Cam) * inv_KaCa));
    Fi = 1.0 + (NADH) * inv_KidhNADH;

    kh_1 = 8.1e-5; % mM ionization constant of IDH
    kh_2 = 5.98e-5; % mM ionization constant of IDH - why separate kh_1, kh_2?
    H = 2.5e-5; % mM Matrix proton concentration
    VIDH_Constant = 1.0 + H / kh_1 + kh_2 / H;

    Kmiso = 1.52; % mM Michaelis constant for isocitrate
    nID = 2; % Dimensionless, cooperativity for isocitrate
    kIDH = 0.05; % mM rate constant of IDH 0.03 in cortassa, 0.05 in cellml
    EtID = 0.109; % mM conc of IDH
    kIDH_EtID = kIDH * EtID;
    
    VIDH = kIDH_EtID / ...
		  (VIDH_Constant + KmIDNAD_NAD * Fi + power( Kmiso/(ISOC), nID) * Fa * (1.0 + KmIDNAD_NAD  * Fi));

    % NADH VKGDH
    Mg = 0.4; % mM Mg2+ conc in mitochondria
    Kmg = 0.0308; % mM activation constant for Mg2+
    Kca = 1.27e-3; % mM activation constant for Ca2+
    Mg_Kmg_1 = Mg / Kmg + 1.0;
    Mg_Kmg_1_Kca = Mg_Kmg_1 / Kca;
    a = ( (Mg_Kmg_1 + Mg_Kmg_1_Kca*(Cam)) );

    kKGDH = 7.5e-2; % ms-1 rate constant of KGDH 0.05 in cortassa, 7.5e-2 in cellml
    EtKG = 0.5; % mM Conc of KGDH
    kKGDH_EtKG = kKGDH * EtKG;

    KmKG = 1.94; % mM Michaelis constant for alphaKG
    nKG = 1.2; % dimensionless, Hill coefficient of KGDH for alphaKG

    KmKGNAD = 38.7; % mM Michaelis constant for NAD
    KmKGNAD_KmIDNAD = KmKGNAD / KmIDNAD;
    
    VKGDH = kKGDH_EtKG * a / ( a + power(KmKG/(AKG),nKG) + KmKGNAD_KmIDNAD * KmIDNAD_NAD );

    % NADH VMDH
    Kh1 = 1.13e-5; % mM ionization constant of MDH
    Kh2 = 26.7; % mM ionization cosnstant of MDH 
    Kh3 = 6.68e-9; % mM ionization constant of MDH
    Kh4 = 5.62e-6; % mM ionization constant of MDH - why all different??
    Koff = 3.99e-2; % dimensionless, pH independent term in pH activation factor of MDH
    kMDH = 0.111; % ms-1 rate constat of MDH
    EtMD = 0.154; % mM total MDH enzyme conc

    kMDH_Fh_EtMD = power( ( 1.0 / ( 1.0 + Kh3 / H + Kh3 * Kh4 / power(H,2) ) ) ,2) * ...
		                ( 1.0 / ( 1.0 + H / Kh1 + power(H,2) / ( Kh1 * Kh2 ) ) + Koff) * ...
						  kMDH * EtMD;
    Kmal = 1.493; % mM Michaelis constant for malate
    Kioaa = 3.1e-3; % mM Inhibition constant for oxalacetate
    Kmal_Kioaa = Kmal / Kioaa;
    KmmNAD = 0.2244; % mM Michaelis constant for NAD+

    VMDH = kMDH_Fh_EtMD * (MAL) * NAD / ...
		( ( (MAL) + Kmal + (Oaa) * Kmal_Kioaa ) * ( KmmNAD + NAD ) );

    % ISOC VACO

    KACOeq = 2.22; % dimensionless, equilibrium constant of ACO
    kfACO = 1.25e-2; % ms-1 forward rate constant of ACO
    CIK = 1; % mM Sum of TCA cycle intermediates - used to back calculate CIT
    one_inv_KACOeq = 1.0 + 1.0 / KACOeq;


    VACO = kfACO * ( CIK - (AKG) - (SCoA) - (Succ) - (FUM) - (MAL) - (Oaa) - (ISOC)* one_inv_KACOeq );

    % ISOC VIDH established in NADH

    % AKG VIDH established in NADH

    % AKG VAAT 
    GLU = 10; % mM glutamate conc
    kfAAT = 6.44e-4; % ms-1 forward rate constant of AAT
    kcnsASP = 1.5e-6; % ms-1 rate constant of aspartate consumption
    KAATeq = 6.6; % dimensionless, equilibrium constant of AAT

    VAAT_Constant = kfAAT * GLU * kcnsASP * KAATeq / kfAAT;
    kcnsASP_KAATeq_kfAAT = kcnsASP * KAATeq / kfAAT;

    VAAT = VAAT_Constant * (Oaa) / ( kcnsASP_KAATeq_kfAAT + (AKG));

    % AKG VKGDH established in NADH

    % SCoA VKGDH established in NADH
    % SCoA VSL established in ADPm

    % Succ VSL established in ADPm
    % Succ VSDH
    
    kSDH = 0.005; % ms-1 rate constant of SDH 3e-3 in cortassa, 0.005 in cellml
    EtSDH = 0.5; % mM SDH enzyme conc
    kSDH_EtSDH = kSDH * EtSDH;
    KmSucc = 0.03; % mM Michaelis constant for succinate
    KiFUM = 1.3; % mM Inhibition constant by fumarate
    KmSucc_KiFUM = KmSucc / KiFUM;
    KiOxaa = 0.15; % mM Inhibition constant by oxalacetate
    inv_KiOxaa = 1.0 / KiOxaa; 

    VSDH = kSDH_EtSDH * (Succ) / ((Succ) + (KmSucc + KmSucc_KiFUM*(FUM)) * (1.0 + inv_KiOxaa*(Oaa)) );

    % FUM VSDH established in Succ
    % FUM VFH
    kfFH = 3.32e-3; % ms-1 forward rate constant for FH
    KFHeq = 1; % dimensionless, equilibrium constnat of FH
    kfFH_KFHeq = kfFH / KFHeq;

    VFH = kfFH * (FUM) - kfFH_KFHeq * (MAL);

    % MAL VFH established in FUM
    % MAL VMDH established in NADH

    % Oaa established in NADH
    % Oaa VCS
    
    KCS = 0.05; % ms-1 catalytic constant of CS
    EtCS = 0.4; % mM concentration of CS
    AcCoA = 1.0; % mM Acetyl CoA conc
    KmAcCoA = 1.26e-2; % mM Michaelis constant for AcCoA
    KmOaa = 6.4e-4; % mM Michaelis constant for OAA

    VCS_C1 = (KCS * EtCS * AcCoA) / (KmAcCoA + AcCoA);

    VCS = VCS_C1 * (Oaa) / ((Oaa) +  KmOaa);

    % Oaa VAAT established in AKG

    % ASP VAAT established in Oaa

    % ATPi V_AM

    V_AM_scaler = 15; % dimensionless, what is this parameter? taken from cellml
    V_AM_max = 0.00048; % mM.ms-1 maximal rate of ATP hydrolysis by myofibrils - in cell ml 0.00048 - in cortassa 7.2e-3
    f_xb = 0.05; % ms-1 transition rate from weak to strong crossbridge
    f_01 = 3.0  * f_xb;
	f_12 = 10.0 * f_xb;
	f_23 = 7.0  * f_xb;
    inv_ATPi = 1.0 / (ATPi);
    KmATP_AM = 0.03; % mM ATP half saturation constant of AM ATPase
    Ki_AM = 0.26; % mM ADP inhibition constant of AM ATPase

    KmATP_AM_Ki_AM = KmATP_AM / Ki_AM;
    V_AM_scaler_max_1_f_01_12_23 = V_AM_scaler * V_AM_max / (f_01 + f_12 + f_23);

    V_AM = V_AM_scaler_max_1_f_01_12_23 * (f_01 * (P0) + f_12 * (P1) + f_23 * (P2)) / ...
		   (1.0 + inv_ATPi * ( KmATP_AM + KmATP_AM_Ki_AM * ADP ));

    % ATPi VCK_mito
    
    kf_3 = 1.33e-6; % ms-1 forward rate constant of mito CK 
    CRT_mito = 25; % mM total conc of creat metabolites - both mito and cyto
    keq = 0.0095; % dimensionless, equilibrium constant of CK 0.0095 in cellml
    inv_keq = 1.0 / keq;

    VCK_mito = kf_3 * ( ( CRT_mito - (CrPi_mito) ) * (ATPi) - (CrPi_mito) * ADP * inv_keq);

    % ATPi_cyto VCK_cyto
    kf_2 = 1.4e-4; % ms-1 forward rate constant of cytoplasmic CK#
    CRT_cyto = 25; % mM total conc of creatinine metabolites - both mito and cyto

    VCK_cyto = kf_2 * ( ( CRT_cyto - (CrPi_cyto) ) * (ATPi_cyto) -  (CrPi_cyto) * ( 8.0 - (ATPi_cyto) ) * inv_keq);

    % ATPi_cyto VATPase_cyto
    VATPase_cyto = 1e-5; % mM.ms-1 constitutive cytosolic ATP consumption rate

    % Crpi_myto VCK_mito established in ATPi

    % Crpi_myto Vt_CRP2

    kt_2 = 2e-3; % ms-1 transfer rate constant of CrP
    Vt_CRP2  = kt_2 * ( (CrPi_mito) - (CrPi_cyto) );

    % P0 kTrop_pn_f_01
    kTrop_pn = 0.04; % ms-1 transition rate from tropomyosin permissive to nonpermissive
    kTrop_pn_f_01 = -kTrop_pn - f_01;

    % P0 kTrop_np

    SL = 2.15; % micrometer, sarcomere length
    LTRPNtot = 0.07; % mM total troponin high affinity sites
    Ntrop = 3.5 * SL - 2.0;
    kltrpn_minus = 4e-2; % ms-1 Ca2+ off rate for troponin low affinity sites
    kltrpn_plus = 100; % mM-1.ms-1 Ca2+ on-rate for troponin low affinity sites

    Ktrop_Ca = kltrpn_minus / kltrpn_plus;
    Ktrop_half = 1.0 / (1.0 + (Ktrop_Ca / (1.7 / 1000.0 + ((0.9 / 1000.0 - 1.7 / 1000.0) / (2.3 - 1.7)) * (SL - 1.7))));
    inv_LTRPNtot_Ktrop_half = 1.0 / ( LTRPNtot * Ktrop_half );

    kTrop_np = kTrop_pn * power( (LTRPNCa) * inv_LTRPNtot_Ktrop_half, Ntrop); 

    % p0 g_01_mod
    gmin_xb = 0.1; % ms-1 minimum transition rate from strong to weak crossbridge
    g0_01 = 1.0 * gmin_xb;
    mod_factor = 1 + (2.3 - SL) / power((2.3 - 1.7),1.6);
    g_01_mod = g0_01 * mod_factor;

    % P1 kTrop_pn_f_12_g_01_mod
    kTrop_pn_f_12_g_01_mod = -(kTrop_pn + f_12 + g_01_mod);

    % P1 kTrop_np and f_01 established in P0

    % P1 g_12_mod
    g0_12 = 2.0 * gmin_xb;
    g_12_mod = g0_12 * mod_factor;

    % P2 f_23_g_12_mod
    f_23_g_12_mod = -(f_23 + g_12_mod);

    % P2 f_12 established in ATPi V_AM, 

    % P2 g_23_mod
    g0_23 = 3.0 * gmin_xb;
    g_23_mod = g0_23 * mod_factor;

    % P3 g_23_mod established in P2 f_23 established in ATPi V_AM

    % N1 kTrop_pn estabished in P0

    % N1 g_01_off_mod
    g_01_off = 30.0 / 1000.0;
    g_01_off_mod = g_01_off * mod_factor;

    % LTRPNCa twoThirds
    twoThirds = 2.0 / 3.0;

    % LTRPNCa FN_Ca
    La = 1; % um from cellml
    Lm_prime = 1.5; % um from cellml
    Lz = 0.1; % um from cellml
    Lb = 0.1; % um
    Lm = Lm_prime-Lb;

    paths = g0_01 * g0_12 * g0_23 + f_01 * g0_12 * g0_23 + f_01 * f_12 * g0_23 + f_01 * f_12 * f_23;
    P1max = (f_01 * (2.0 * gmin_xb) * (3.0 * gmin_xb)) / paths;
    P2max = (f_01 * f_12 * (3.0 * gmin_xb)) / paths;
    P3max = (f_01 * f_12 * f_23) / paths;

    fnormmax2 = P1max+P2max+P3max; % 

    if SL<2.2
       alpha_SL = min(1.0, (SL - 2.0 * La + Lm_prime - Lz) / Lm); % this might take long to execute
    else
       alpha_SL = 1 - (SL - 2.2) / Lm;
    end

    %alpha_SL = 1 - (SL - 2.2) / Lm;
    alpha_SL_fnormmax2 = alpha_SL / fnormmax2;
    P1_N1_P2_P3 = P1 + N1 + P2 + P3;
    FN_Ca = alpha_SL_fnormmax2 * P1_N1_P2_P3; 

    % HTRPNCa constants
    khtrpn_plus = 100; % mM-1.ms-1 Ca2+ on-rate for troponin high affinity sites
    HTRPNtot = 0.14; % mM total troponin high affinity sites
    khtrpn_minus = 3.3e-4; % ms-1 Ca2+ off rate for troponin high affinity sites

    %%%%%%%%% force - not normalised from cortassa 2006 A121 %%%%%%%%
    zeta = 0.1; % N/mm2 conversion factor normalizing to physiological force
    Fmax  = P1max + 2.0 * P2max + 3 * P3max;
    fnormmax = Fmax / 3.0;
    alpha_SL_fnormmax  = alpha_SL / fnormmax / 3.0;
    zeta_alpha_SL_fnormmax = zeta * alpha_SL_fnormmax;

    force_norm = alpha_SL_fnormmax*(P1_N1_P2_P3+P2+P3+P3);

    force = zeta_alpha_SL_fnormmax*(P1_N1_P2_P3+P2+P3+P3);

    % derivatives and outputs
    two_b = b*2; % dimensionless
    
    %CIT = CIK - (ISOC + AKG + SCoA + Succ + FUM + MAL + Oaa); % mM citrate conc

	dCam  = fm * (Vuni - VnaCa);
	dADPm = VANT - VATPase - VSL;
	dDpsi = -(-VHNe - VHFe + Vhu + VANT + Vhleak + two_b * VnaCa + 2.0 * Vuni ) * inv_Cmito; % mito membrane potential - two_b is 1 anyway
	dNADH = -VNO + VIDH + VKGDH + VMDH;
	dISOC = VACO - VIDH;
	dAKG  = VIDH + VAAT - VKGDH;
	dSCoA = VKGDH - VSL;
	dSucc = VSL - VSDH;
	dFUM  = VSDH - VFH;
	dMAL  = VFH - VMDH;
	dOaa  = VMDH - VCS - VAAT;
	dASP  = VAAT - kcnsASP * (ASP); % I think ASP is aspartate? but it's not in the main paper

    % This bit is force generation model

    N0 = 1-(N1+P0+P1+P2+P3);

	dATPi = 0.615 * VANT - V_AM - 0.05*Jup - ( 6.371e-5 * ( INaK + IpCa ) ) - VCK_mito; %% 0.05 scalar to Jup artificially added!
	dATPi_cyto = - VCK_cyto - VATPase_cyto;
    dCrPi_mito = VCK_mito - Vt_CRP2;
	dCrPi_cyto = Vt_CRP2 + VCK_cyto;
    dP0 = kTrop_pn_f_01 * (P0) + kTrop_np * (N0) + g_01_mod * (P1);
	dP1 = kTrop_pn_f_12_g_01_mod * (P1) + kTrop_np * (N1) + f_01 * (P0) + g_12_mod * (P2);
	dP2 = f_23_g_12_mod * (P2) + f_12 * (P1) + g_23_mod * (P3);
	dP3 = -g_23_mod * (P3) + f_23 * (P2);
	dN1 = kTrop_pn * (P1) - (kTrop_np + g_01_off_mod) * (N1);

    dLTRPNCa = kltrpn_plus * (Cai) * (LTRPNtot - LTRPNCa)- ( kltrpn_minus * LTRPNCa ) * ( 1.0 - twoThirds * FN_Ca );
    dHTRPNCa = khtrpn_plus * (Cai) * (HTRPNtot - (HTRPNCa)) - khtrpn_minus * (HTRPNCa);
	
end

