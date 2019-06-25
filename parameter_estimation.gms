* Parameter estimation problem with stability constraints for Kuramoto model 
* of metabolic and cell cycle oscillators 
* 13th June 2019
* Serdar Ozsezen


option reslim=30;

sets
cid conditions
sid statistics
samp /1*1000/
osc /MET, G1, S, M/
i /1*4/
m /1*18/;
;

alias(i,j,k);
alias(osc,osc1);
alias(osc,osc2);
alias(cid,cid1);

parameters
S_omega_comp(sid,cid)
S_omega_met(sid,cid)
S_dphi_G1(sid,cid)
S_dphi_S(sid,cid)
S_dphi_M(sid,cid)

MS_K(samp,osc1,osc2)
MS_omega(samp,cid,osc)
MS_L2norm(samp)
MS_L2norm_noA(samp)
MS_omega_comp(samp,cid)
MS_dphi(samp,cid,osc)
MS_RSS(samp)
MS_jac(samp,cid,i,j)
MS_a0(samp,cid)
MS_a1(samp,cid)
MS_a2(samp,cid)
;

* Load the experimental data 

$GDXIN input.gdx
$load sid cid
$load S_omega_comp S_omega_met S_dphi_G1 S_dphi_S S_dphi_M
$GDXIN

parameters

S(cid,i,j) For stability calculations
St(cid,i,j)
tst(cid,i,j)
;

S(cid,'1','1') = 0.5;
S(cid,'1','2') = -0.5;
S(cid,'1','3') = -0.5;
S(cid,'1','4') = -0.5;
S(cid,'2','1') = 0.5;
S(cid,'2','2') = 0.833333333333333;
S(cid,'2','3') = -0.166666666666667;
S(cid,'2','4') = -0.166666666666667;
S(cid,'3','1') = 0.5;
S(cid,'3','2') = -0.166666666666667;
S(cid,'3','3') = 0.833333333333333;
S(cid,'3','4') = -0.166666666666667;
S(cid,'4','1') = 0.5;
S(cid,'4','2') = -0.166666666666667;
S(cid,'4','3') = -0.166666666666667;
S(cid,'4','4') = 0.833333333333333;

St(cid,i,j) = S(cid,j,i);
tst(cid,i,j) = sum(k, S(cid,i,k)*St(cid,k,j));

variable
V_L2norm squared loss
V_dphi(cid,osc) the phase shift between metabolism and the cell cycle phase
V_alphaK

V_jac(cid,i,j) Jacobian matrix
V_Sj(cid,i,j)
V_Jt(cid,i,j) For submatrix

V_a2(cid) Coefficients
V_a1(cid)
V_a0(cid)
V_stab(cid)
;

positive variable
V_omega_comp(cid) model angular compromise frequency of the states per hour
V_omega(cid,osc) model angular frequency of the states per hour
V_K(osc1,osc2) Coupling matrix
V_omegasame(osc)
V_zero(cid) For estimating the left out condition
;

scalar
alpha
;



equation
E_L2norm squared error
E_oscillators oscillators for the different states at different conditions
E_sameomega
E_alphaK
E_zero

E_jac1 Jacobian matrix elements
E_jac2
E_jac3
E_jac4
E_jac5
E_jac6
E_jac7
E_jac8
E_jac9
E_jac10
E_jac11
E_jac12
E_jac13
E_jac14
E_jac15
E_jac16

E_Sj Stability calculations
E_Jt
E_a2
E_a1
E_a0

E_stab
;

* Objective function
E_L2norm..
V_L2norm =e= sum(cid,
         sqr((V_omega_comp(cid)-S_omega_comp('mean',cid))/S_omega_comp('sd',cid))
         +sqr((V_omega(cid,'MET')-S_omega_met('mean',cid))/S_omega_met('sd',cid))
         +sqr((V_dphi(cid,'G1')-S_dphi_G1('mean',cid))/S_dphi_G1('sd',cid))
         +sqr((V_dphi(cid,'S')-S_dphi_S('mean',cid))/S_dphi_S('sd',cid))
         +sqr((V_dphi(cid,'M')-S_dphi_M('mean',cid))/S_dphi_M('sd',cid))
) + V_alphaK;

* Penalty term
E_alphaK.. V_alphaK =e= alpha*sum(osc1,sum(osc2,V_K(osc1,osc2)));

* Kuramoto model equation
E_oscillators(cid,osc)..
V_omega_comp(cid) =e= V_omega(cid,osc) - sum(osc1$(not sameas(osc,osc1)),V_K(osc,osc1)*(sin(V_dphi(cid,osc) - V_dphi(cid,osc1))));

* Natural frequencies of the cell cycles should be the same regardless of nutrient conditions
E_sameomega(cid,osc)$(not sameas(osc,'MET')).. V_omegasame(osc) =e= V_omega(cid,osc);

* Jacobian calculations
E_jac1(cid,i,j).. V_jac(cid, '1','1') =e= - V_K('MET','G1')*cos(V_dphi(cid,'MET') - V_dphi(cid,'G1'))
                                          - V_K('MET','S')*cos(V_dphi(cid,'MET') - V_dphi(cid,'S'))
                                          - V_K('MET','M')*cos(V_dphi(cid,'MET') - V_dphi(cid,'M'));

E_jac2(cid,i,j).. V_jac(cid, '2','2') =e= - V_K('G1','MET')*cos(V_dphi(cid,'G1') - V_dphi(cid,'MET'))
                                          - V_K('G1','S')*cos(V_dphi(cid,'G1') - V_dphi(cid,'S'))
                                          - V_K('G1','M')*cos(V_dphi(cid,'G1') - V_dphi(cid,'M'));

E_jac3(cid,i,j).. V_jac(cid, '3','3') =e= - V_K('S','MET')*cos(V_dphi(cid,'S') - V_dphi(cid,'MET'))
                                          - V_K('S','G1')*cos(V_dphi(cid,'S') - V_dphi(cid,'G1'))
                                          - V_K('S','M')*cos(V_dphi(cid,'S') - V_dphi(cid,'M'));

E_jac4(cid,i,j).. V_jac(cid, '4','4') =e= - V_K('M','MET')*cos(V_dphi(cid,'M') - V_dphi(cid,'MET'))
                                          - V_K('M','G1')*cos(V_dphi(cid,'M') - V_dphi(cid,'G1'))
                                          - V_K('M','S')*cos(V_dphi(cid,'M') - V_dphi(cid,'S'));

E_jac5(cid,i,j).. V_jac(cid,'1','2') =e= V_K('MET','G1')*cos(V_dphi(cid,'MET') - V_dphi(cid,'G1'));
E_jac6(cid,i,j).. V_jac(cid,'1','3') =e= V_K('MET','S')*cos(V_dphi(cid,'MET') - V_dphi(cid,'S'));
E_jac7(cid,i,j).. V_jac(cid,'1','4') =e= V_K('MET','M')*cos(V_dphi(cid,'MET') - V_dphi(cid,'M'));

E_jac8(cid,i,j).. V_jac(cid,'2','1') =e= V_K('G1','MET')*cos(V_dphi(cid,'G1') - V_dphi(cid,'MET'));
E_jac9(cid,i,j).. V_jac(cid,'2','3') =e= V_K('G1','S')*cos(V_dphi(cid,'G1') - V_dphi(cid,'S'));
E_jac10(cid,i,j).. V_jac(cid,'2','4') =e= V_K('G1','M')*cos(V_dphi(cid,'G1') - V_dphi(cid,'M'));

E_jac11(cid,i,j).. V_jac(cid,'3','1') =e= V_K('S','MET')*cos(V_dphi(cid,'S') - V_dphi(cid,'MET'));
E_jac12(cid,i,j).. V_jac(cid,'3','2') =e= V_K('S','G1')*cos(V_dphi(cid,'S') - V_dphi(cid,'G1'));
E_jac13(cid,i,j).. V_jac(cid,'3','4') =e= V_K('S','M')*cos(V_dphi(cid,'S') - V_dphi(cid,'M'));

E_jac14(cid,i,j).. V_jac(cid,'4','1') =e= V_K('M','MET')*cos(V_dphi(cid,'M') - V_dphi(cid,'MET'));
E_jac15(cid,i,j).. V_jac(cid,'4','2') =e= V_K('M','G1')*cos(V_dphi(cid,'M') - V_dphi(cid,'G1'));
E_jac16(cid,i,j).. V_jac(cid,'4','3') =e= V_K('M','S')*cos(V_dphi(cid,'M') - V_dphi(cid,'S'));

* Routh–Hurwitz stability criterion

E_Sj(cid,i,j).. V_Sj(cid,i,j) =e= sum(k,St(cid,i,k)*V_jac(cid,k,j));

E_Jt(cid,i,j).. V_Jt(cid,i,j) =e= sum(k,V_Sj(cid,i,k)*S(cid,k,j));

* Coefficients for stability

E_a2(cid).. V_a2(cid) =e= - V_Jt(cid,'2','2') - V_Jt(cid,'3','3') - V_Jt(cid,'4','4');

E_a1(cid).. V_a1(cid) =e= V_Jt(cid,'2','2')*V_Jt(cid,'3','3') - V_Jt(cid,'2','3')*V_Jt(cid,'3','2') + V_Jt(cid,'2','2')*V_Jt(cid,'4','4')
- V_Jt(cid,'2','4')*V_Jt(cid,'4','2') + V_Jt(cid,'3','3')*V_Jt(cid,'4','4') - V_Jt(cid,'3','4')*V_Jt(cid,'4','3');

E_a0(cid).. V_a0(cid) =e= V_Jt(cid,'2','2')*V_Jt(cid,'3','4')*V_Jt(cid,'4','3') - V_Jt(cid,'2','2')*V_Jt(cid,'3','3')*V_Jt(cid,'4','4')
+ V_Jt(cid,'2','3')*V_Jt(cid,'3','2')*V_Jt(cid,'4','4') - V_Jt(cid,'2','3')*V_Jt(cid,'3','4')*V_Jt(cid,'4','2')
- V_Jt(cid,'2','4')*V_Jt(cid,'3','2')*V_Jt(cid,'4','3') + V_Jt(cid,'2','4')*V_Jt(cid,'3','3')*V_Jt(cid,'4','2');

E_stab(cid)$(sameas(cid,'PYR') or sameas(cid,'GLU_H')).. V_a2(cid)*V_a1(cid) - V_a0(cid) =g= 1e-5;

model regression /E_L2norm, E_oscillators, E_sameomega, E_alphaK, E_jac1, E_jac2, E_jac3, E_jac4, E_jac5, E_jac6, E_jac7, E_jac8, E_jac9, E_jac10, E_jac11,
E_jac12, E_jac13, E_jac14, E_jac15, E_jac16, E_Sj, E_Jt, E_a2, E_a1, E_a0, E_stab/;


* We constraint the variables to be within 5 standard deviations of the sample mean
V_omega_comp.lo(cid) = S_omega_comp('mean',cid)-5*S_omega_comp('sd',cid);
V_omega_comp.up(cid) = S_omega_comp('mean',cid)+5*S_omega_comp('sd',cid);
V_omega.lo(cid,'MET') = S_omega_met('mean',cid)-5*S_omega_met('sd',cid);
V_omega.up(cid,'MET') = S_omega_met('mean',cid)+5*S_omega_met('sd',cid);
V_omega.up(cid,osc)$(not sameas(osc,'MET'))  = 1e2;

* We define dphi to be the phase difference relative to the metabolism state
V_dphi.fx(cid,'MET') = 0;
V_dphi.lo(cid,'G1') = S_dphi_G1('mean',cid)-5*S_dphi_G1('sd',cid);
V_dphi.up(cid,'G1') = S_dphi_G1('mean',cid)+5*S_dphi_G1('sd',cid);
V_dphi.lo(cid,'S') = S_dphi_S('mean',cid)-5*S_dphi_S('sd',cid);
V_dphi.up(cid,'S') = S_dphi_S('mean',cid)+5*S_dphi_S('sd',cid);
V_dphi.lo(cid,'M') = S_dphi_M('mean',cid)-5*S_dphi_M('sd',cid);
V_dphi.up(cid,'M') = S_dphi_M('mean',cid)+5*S_dphi_M('sd',cid);
V_dphi.l(cid,'S') = S_dphi_S('mean',cid);
V_dphi.l(cid,'G1') = S_dphi_G1('mean',cid);
V_dphi.l(cid,'M') = S_dphi_M('mean',cid);


V_a2.lo('PYR') = 1e-5;
V_a0.lo('PYR') = 1e-5;
V_a2.lo('GLU_H') = 1e-5;
V_a0.lo('GLU_H') = 1e-5;

* Non-linear programming solver is CONOPT

$onecho > conopt.opt
LFILOG 10
LFILOS 10
$offecho

option nlp=conopt;

regression.solvelink = 5;
regression.holdfixed = 1;
regression.optfile = 1;

execseed = 1 + gmillisec(jnow);

* Regularization factor
alpha = 0.01;

MS_L2norm_noA(samp) = NA;

loop(samp,
  V_omega_comp.l(cid) = S_omega_comp('mean',cid);
  V_omega.l(cid,'MET') = S_omega_met('mean',cid);
  V_omega.l(cid,osc)$(not sameas(osc,'MET'))  =  10**(uniform(-5,2));
  V_dphi.l(cid,'G1') = S_dphi_G1('mean',cid);
  V_dphi.l(cid,'S') = S_dphi_S('mean',cid);
  V_dphi.l(cid,'M') = S_dphi_M('mean',cid);
  V_K.l(osc1,osc2)$(not sameas(osc1,osc2)) =  10**(uniform(-1,7));
  V_jac.l(cid,i,j) = 0;
  V_Sj.l(cid,i,j) = 0;
  V_Jt.l(cid,i,j) = 0;
  V_a0.l(cid) = 5;
  V_a1.l(cid) = 5;
  V_a2.l(cid) = 5;
  solve regression us nlp min V_L2norm;
  if((regression.modelstat eq 1 or regression.modelstat eq 2 or regression.modelstat eq 3 or regression.modelstat eq 7),
    MS_K(samp,osc1,osc2) = V_K.l(osc1,osc2);
    MS_omega(samp,cid,osc) = V_omega.l(cid,osc);
    MS_L2norm(samp) = V_L2norm.l;
    MS_L2norm_noA(samp) = V_L2norm.l - V_alphaK.l;
    MS_omega_comp(samp,cid) = V_omega_comp.l(cid);
    MS_dphi(samp,cid,osc) = V_dphi.l(cid,osc);
    MS_jac(samp,cid,i,j) = V_jac.l(cid,i,j);
    MS_a0(samp,cid) = V_a0.l(cid);
    MS_a1(samp,cid) = V_a1.l(cid);
    MS_a2(samp,cid) = V_a2.l(cid);
    MS_RSS(samp) = sum(cid, sqr(S_dphi_G1('mean',cid) - V_dphi.l(cid, 'G1')) +
                            sqr(S_dphi_S('mean',cid) - V_dphi.l(cid, 'S')) +
                            sqr(S_dphi_M('mean',cid) - V_dphi.l(cid, 'M')) +
                            sqr(S_omega_met('mean',cid) - V_omega.l(cid,'MET')) +
                            sqr(S_omega_comp('mean',cid) - V_omega_comp.l(cid))
                      );
  );
  );

execute_unload './regression_results.gdx';