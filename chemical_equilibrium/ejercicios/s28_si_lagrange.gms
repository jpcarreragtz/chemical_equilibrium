$offListing
$offSymXRef offSymList
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
*   4) R minimo: CNS con R como variable y Delta_Si del gas solo = 0.
* Variable: ln(n_i) del gas, sin limite inferior en 0; n_Si libre.
* Salida: el .lst solo trae los display finales; la tabla de resultados sale en la
* terminal y en s28_si_lagrange_resultados.txt. IPOPT no escribe en la terminal.
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

* IPOPT sin texto en la terminal ni errores de evaluacion en el .lst: las opciones
* van en ipopt.opt, que se escribe aqui mismo en la carpeta de trabajo.
$onEcho > ipopt.opt
print_level 0
sb yes
print_eval_error no
$offEcho

Set sv 'solves de cada punto' / NLP, NLP_reintento, CNS_gas, CNS_sol / ;

Parameters  n0(i), b(e), gRT(i), lnn_gas(i), lam_gas(e), lnNt_gas,
            nSi(cl,r)    'mol de Si(s) por mol de clorosilano'
            delta(cl,r)  'prueba de fase con el gas solo'
            yR(cl,r,i)   'fraccion mol del gas'
            lamR(cl,r,e) 'lambda entre RT final'
            estado(cl,r,*) 'model status: NLP, CNS gas, CNS con Si(s)'
            hecho(cl,r,sv) '1 si ese solve se corrio'
            bien(cl,r,sv)  '1 si salio bien (NLP: model status 1 o 2; CNS: 15 a 17)'
            lnnG(cl,r,i), lnNtG(cl,r), lamG(cl,r,e)  'gas solo de cada punto (arranque de R minimo)'
            bcl(e)       'atomos que aporta 1 mol de clorosilano'
            Rmin(cl)     'R con Delta_Si del gas solo = 0'
            estadoRmin(cl) 'model status del CNS de R minimo'
            bienRmin(cl) ;
Scalars  Rgas 'J por mol K' / 8.314 /, Tk / 1300 /, P / 1 /, P0 / 1 /, gSi 'dGf de Si(s) / RT' / 0 / ;

Variables  lnn(i), lnNt, z, lam(e), nS, Rv 'R como variable (R minimo)' ;
Equations  balance(e), balance_sol(e), balance_R(e), total, objetivo, potencial(i), solido ;

balance(e)..      sum(i, A(e,i)*exp(lnn(i))) =e= b(e);
balance_sol(e)..  sum(i, A(e,i)*exp(lnn(i))) + ASi(e)*nS =e= b(e);
balance_R(e)..    sum(i, A(e,i)*exp(lnn(i))) =e= bcl(e) + A(e,'H2')*Rv;
total..           sum(i, exp(lnn(i))) =e= exp(lnNt);
objetivo..        z =e= sum(i, exp(lnn(i))*(gRT(i) + lnn(i) - lnNt + log(P/P0)));
potencial(i)..    gRT(i) + lnn(i) - lnNt + log(P/P0) + sum(e, A(e,i)*lam(e)) =e= 0;
solido..          gSi + sum(e, ASi(e)*lam(e)) =e= 0;

lnn.lo(i) = -700;

Model gibbs   / balance, total, objetivo / ;
Model kkt_gas / balance, total, potencial / ;
Model kkt_sol / balance_sol, total, potencial, solido / ;
Model kkt_rmin / balance_R, total, potencial, solido / ;
* solprint = silent quita tambien el resumen de cada solve (con off se queda).
option nlp = ipopt, cns = ipopt, limrow = 0, limcol = 0, solprint = silent;
gibbs.optfile = 1;  kkt_gas.optfile = 1;  kkt_sol.optfile = 1;  kkt_rmin.optfile = 1;

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
   hecho(cl,r,'NLP') = 1;
   bien(cl,r,'NLP')  = (gibbs.modelstat <= 2 and gibbs.solvestat = 1);
*  Con R grande el NLP desde el reparto uniforme a veces no llega (IPOPT: "local
*  infeasibility"); entonces se reintenta desde el gas del R anterior.
   if(gibbs.modelstat > 2 and ord(r) > 1,
      lnn.l(i) = lnn_gas(i);  lnNt.l = lnNt_gas;
      solve gibbs using nlp minimizing z;
      estado(cl,r,'NLP_reintento') = gibbs.modelstat;
      hecho(cl,r,'NLP_reintento') = 1;
      bien(cl,r,'NLP_reintento')  = (gibbs.modelstat <= 2 and gibbs.solvestat = 1);
   );
   lam.l(e) = -balance.m(e);
   solve kkt_gas using cns;
   estado(cl,r,'CNS_gas') = kkt_gas.modelstat;
   hecho(cl,r,'CNS_gas') = 1;
   bien(cl,r,'CNS_gas')  = (kkt_gas.modelstat >= 15 and kkt_gas.modelstat <= 17 and kkt_gas.solvestat = 1);
   lnn_gas(i) = lnn.l(i);  lnNt_gas = lnNt.l;  lam_gas(e) = lam.l(e);
   lnnG(cl,r,i) = lnn.l(i);  lnNtG(cl,r) = lnNt.l;  lamG(cl,r,e) = lam.l(e);
*  2) prueba de fase
   delta(cl,r) = gSi + sum(e, ASi(e)*lam.l(e));
   nSi(cl,r) = 0;
*  3) Si(s) si baja G
   if(delta(cl,r) < 0,
      nS.l = 0;
      solve kkt_sol using cns;
      estado(cl,r,'CNS_sol') = kkt_sol.modelstat;
      hecho(cl,r,'CNS_sol') = 1;
      bien(cl,r,'CNS_sol')  = (kkt_sol.modelstat >= 15 and kkt_sol.modelstat <= 17 and kkt_sol.solvestat = 1);
      nSi(cl,r) = nS.l;
      if(nS.l < 0,
         nSi(cl,r) = 0;
         lnn.l(i) = lnn_gas(i);  lnNt.l = lnNt_gas;  lam.l(e) = lam_gas(e);
      );
   );
   yR(cl,r,i)   = exp(lnn.l(i) - lnNt.l);
   lamR(cl,r,e) = lam.l(e);
);

*  4) R minimo: donde Delta_Si del gas solo cambia de signo entre dos R seguidos,
*     CNS con R como variable y lambda_Si/RT = -gSi (Delta_Si = 0), desde el gas de ese R.
Rmin(cl) = 0;  estadoRmin(cl) = 0;  bienRmin(cl) = 0;
loop((cl,r)$(ord(r) < card(r) and delta(cl,r)*delta(cl,r+1) < 0),
   bcl(e) = A(e,cl);
   lnn.l(i) = lnnG(cl,r,i);  lnNt.l = lnNtG(cl,r);  lam.l(e) = lamG(cl,r,e);  Rv.l = r.val;
   solve kkt_rmin using cns;
   Rmin(cl) = Rv.l;
   estadoRmin(cl) = kkt_rmin.modelstat;
   bienRmin(cl) = (kkt_rmin.modelstat >= 15 and kkt_rmin.modelstat <= 17 and kkt_rmin.solvestat = 1);
);

option nSi:5:1:1, delta:4:1:1, Rmin:4;
display nSi, delta, Rmin, estado;

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

* ---------------------------------------------------------------------------------
* Reporte de formato fijo, con los titulos de s28_si_epitaxial_lagrange.py: se guarda
* junto a este archivo y se copia a la terminal. Solo cuando este es el archivo que
* se corre (no cuando otro .gms lo incluye, como la revision automatica de VS Code).
$ifThen "%system.incParent%" == ""
Scalars  nTot 'solves corridos', nBien 'solves que salieron bien' ;
nTot  = sum((cl,r,sv), hecho(cl,r,sv)) + sum(cl$estadoRmin(cl), 1);
nBien = sum((cl,r,sv), bien(cl,r,sv))  + sum(cl, bienRmin(cl));

file frep / '%system.fp%%system.fn%_resultados.txt' /;
frep.pw = 200;
put frep;
loop(cl,
   put / '====================================================================================' /
       cl.tl:0 ' + R H2 a 1300 K y 1 bar: prueba de fase y n_Si(s)   [GAMS, IPOPT]' /
       '====================================================================================' /
       '    R  Δ_Si gas solo              decisión    n_Si(s)   λ_Si/RT   λ_Cl/RT   λ_H/RT' /;
   loop(r,
      put r.tl:>5 delta(cl,r):15:4;
      if(delta(cl,r) >= 0,     put 'solo gas (Δ > 0)':>23;
      elseif nSi(cl,r) > 0,    put 'Si(s) entra':>22;
      else                     put 'solo gas (n_Si < 0)':>22;
      );
      put nSi(cl,r):11:5 lamR(cl,r,'Si'):10:4 lamR(cl,r,'Cl'):10:4 lamR(cl,r,'H'):9:4 /;
   );
);
put / '====================================================================================' /
      'R mínimo para depositar Si(s) (Δ_Si del gas solo = 0)' /
      '====================================================================================' /;
loop(cl,
   if(estadoRmin(cl) > 0,
      put 'Desde ' cl.tl:0 ': R_min = ' Rmin(cl):0:4 /;
   elseif smax(r, delta(cl,r)) < 0,
      put 'Desde ' cl.tl:0 ', Δ_Si < 0 en todo el intervalo: siempre hay depósito.' /;
   else
      put 'Desde ' cl.tl:0 ', Δ_Si > 0 en todo el intervalo: no hay depósito.' /;
   );
);
put / 'Estado: ' nBien:0:0 ' de ' nTot:0:0 ' solves salieron bien.' /;
loop((cl,r,sv)$(hecho(cl,r,sv) and not bien(cl,r,sv)),
   put '  No salió bien: ' sv.tl:0 ', ' cl.tl:0 ', R = ' r.tl:0 ' (model status ' estado(cl,r,sv):0:0 ')';
   if(sameas(sv,'NLP') and bien(cl,r,'NLP_reintento'),
      put '; el reintento desde el gas del R anterior sí salió bien';
   );
   put '.' /;
);
loop(cl$(estadoRmin(cl) > 0 and not bienRmin(cl)),
   put '  No salió bien: CNS de R mínimo, ' cl.tl:0 ' (model status ' estadoRmin(cl):0:0 ').' /;
);
putclose frep;
execute 'cat "%system.fp%%system.fn%_resultados.txt"';
$endIf
