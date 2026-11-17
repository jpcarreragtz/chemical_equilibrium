$title Si epitaxial por Lagrange con prueba de fase - Sesion 28 (1300 K, 1 bar)
* Ejecutar:  gams s28_si_lagrange.gms        (opcional: --salida=archivo.csv)
* Si la terminal dice "command not found: gams" (no esta en el PATH), usar la ruta completa:
*    /Library/Frameworks/GAMS.framework/Versions/Current/Resources/gams s28_si_lagrange.gms
* Usa IPOPT (corre sin restriccion con la licencia Free Personal de GAMS).
* Por cada clorosilano y R = mol H2 / mol clorosilano:
*   1) gas solo: NLP (minimiza nG/RT) y CNS con las condiciones de Lagrange;
*   2) prueba de fase con los lambda del gas: Delta_Si = gSi/RT + lambda_Si/RT;
*   3) si Delta_Si < 0, CNS con Si(s): su condicion gSi/RT + lambda_Si/RT = 0 y
*      n_Si en el balance de Si; si n_Si sale negativo, se queda el gas solo.
* Variable: ln(n_i) del gas, sin limite inferior en 0; n_Si libre.
* Datos: dGf a 1300 K (NIST-JANAF, como ejercicios/s24_ej3_si_epitaxial.py), J/mol.

Sets
   i   'especies gaseosas' / SiCl2, SiCl4, SiCl3H, SiCl2H2, SiClH3, SiH4, H2, HCl /
   cl(i) 'clorosilanos alimentados' / SiCl4, SiCl3H /
   r   'R = mol H2 / mol clorosilano' / 1, 2, 3, 5, 7, 10, 15, 20, 30, 50, 75, 100, 125, 150 /
   e   'elementos' / Si, Cl, H / ;

Parameters
   dGf(i) 'J por mol a 1300 K'
      / SiCl2 -216012, SiCl4 -492536, SiCl3H -356537, SiCl2H2 -199368,
        SiClH3 -28482, SiH4 151897, H2 0, HCl -102644 /
   ASi(e) 'atomos de Si(s)' / Si 1 / ;

Table A(e,i) 'atomos del elemento e por molecula de i'
        SiCl2  SiCl4  SiCl3H  SiCl2H2  SiClH3  SiH4  H2  HCl
   Si     1      1      1        1       1      1    0    0
   Cl     2      4      3        2       1      0    0    1
   H      0      0      1        2       3      4    2    1 ;

Parameters  n0(i), b(e), gRT(i), lnn_gas(i), lam_gas(e), lnNt_gas,
            nSi(cl,r)    'mol de Si(s) por mol de clorosilano'
            delta(cl,r)  'prueba de fase con el gas solo'
            yR(cl,r,i)   'fraccion mol del gas'
            lamR(cl,r,e) 'lambda entre RT final'
            estado(cl,r,*) 'model status: NLP, CNS gas, CNS con Si(s)' ;
Scalars  Rgas 'J por mol K' / 8.314 /, Tk / 1300 /, P / 1 /, P0 / 1 /, gSi 'dGf de Si(s) / RT' / 0 / ;

Variables  lnn(i), lnNt, z, lam(e), nS ;
Equations  balance(e), balance_sol(e), total, objetivo, potencial(i), solido ;

balance(e)..      sum(i, A(e,i)*exp(lnn(i))) =e= b(e);
balance_sol(e)..  sum(i, A(e,i)*exp(lnn(i))) + ASi(e)*nS =e= b(e);
total..           sum(i, exp(lnn(i))) =e= exp(lnNt);
objetivo..        z =e= sum(i, exp(lnn(i))*(gRT(i) + lnn(i) - lnNt + log(P/P0)));
potencial(i)..    gRT(i) + lnn(i) - lnNt + log(P/P0) + sum(e, A(e,i)*lam(e)) =e= 0;
solido..          gSi + sum(e, ASi(e)*lam(e)) =e= 0;

lnn.lo(i) = -700;

Model gibbs   / balance, total, objetivo / ;
Model kkt_gas / balance, total, potencial / ;
Model kkt_sol / balance_sol, total, potencial, solido / ;
option nlp = ipopt, cns = ipopt, solprint = off;

gRT(i) = dGf(i)/(Rgas*Tk);
lnn_gas(i) = 0;  lnNt_gas = 0;  lam_gas(e) = 0;
loop((cl,r),
   n0(i) = 0;  n0(cl) = 1;  n0('H2') = r.val;
   b(e) = sum(i, A(e,i)*n0(i));
*  1) gas solo
   lnn.l(i) = log(sum(e, b(e))/(3*card(i)));
   lnNt.l   = log(sum(i, exp(lnn.l(i))));
   solve gibbs using nlp minimizing z;
   estado(cl,r,'NLP') = gibbs.modelstat;
*  Con R grande el NLP desde el reparto uniforme a veces no llega (IPOPT: "local
*  infeasibility"); entonces se reintenta desde el gas del R anterior.
   if(gibbs.modelstat > 2 and ord(r) > 1,
      lnn.l(i) = lnn_gas(i);  lnNt.l = lnNt_gas;
      solve gibbs using nlp minimizing z;
      estado(cl,r,'NLP_reintento') = gibbs.modelstat;
   );
   lam.l(e) = -balance.m(e);
   solve kkt_gas using cns;
   estado(cl,r,'CNS_gas') = kkt_gas.modelstat;
   lnn_gas(i) = lnn.l(i);  lnNt_gas = lnNt.l;  lam_gas(e) = lam.l(e);
*  2) prueba de fase
   delta(cl,r) = gSi + sum(e, ASi(e)*lam.l(e));
   nSi(cl,r) = 0;
*  3) Si(s) si baja G
   if(delta(cl,r) < 0,
      nS.l = 0;
      solve kkt_sol using cns;
      estado(cl,r,'CNS_sol') = kkt_sol.modelstat;
      nSi(cl,r) = nS.l;
      if(nS.l < 0,
         nSi(cl,r) = 0;
         lnn.l(i) = lnn_gas(i);  lnNt.l = lnNt_gas;  lam.l(e) = lam_gas(e);
      );
   );
   yR(cl,r,i)   = exp(lnn.l(i) - lnNt.l);
   lamR(cl,r,e) = lam.l(e);
);

option nSi:5:1:1, delta:4:1:1;
display nSi, delta, estado;

$ifThen set salida
file fcsv / '%salida%' /;
fcsv.nr = 2; fcsv.nd = 12; fcsv.pw = 32767; fcsv.nz = 0;
put fcsv 'caso,clave,valor' /;
loop((cl,r),
   put cl.tl:0 '|' r.tl:0 ',nSi,' nSi(cl,r):22:14 / cl.tl:0 '|' r.tl:0 ',delta,' delta(cl,r):22:14 /;
   loop(i, put cl.tl:0 '|' r.tl:0 ',y_' i.tl:0 ',' yR(cl,r,i):22:14 /);
   put cl.tl:0 '|' r.tl:0 ',estado_NLP,' estado(cl,r,'NLP'):22:14 /
       cl.tl:0 '|' r.tl:0 ',estado_NLP_reintento,' estado(cl,r,'NLP_reintento'):22:14 /
       cl.tl:0 '|' r.tl:0 ',estado_CNS_gas,' estado(cl,r,'CNS_gas'):22:14 /
       cl.tl:0 '|' r.tl:0 ',estado_CNS_sol,' estado(cl,r,'CNS_sol'):22:14 /;
);
putclose fcsv;
$endIf
