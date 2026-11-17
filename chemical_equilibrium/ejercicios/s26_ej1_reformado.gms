$title Minimizacion de Gibbs - Sesion 26, Ejercicio 1 (reformado 4:1, 1 bar, 600-1100 K)
* Ejecutar:  gams s26_ej1_reformado.gms        (opcional: --salida=archivo.csv)
* Si la terminal dice "command not found: gams" (no esta en el PATH), usar la ruta completa:
*    /Library/Frameworks/GAMS.framework/Versions/Current/Resources/gams s26_ej1_reformado.gms
* Usa IPOPT (corre sin restriccion con la licencia Free Personal de GAMS).
* Variable: ln(n_i), sin limite inferior en 0. Dos etapas en cada T (como s25_ej2):
*   1) NLP: minimiza nG/RT.  2) CNS: condiciones de Lagrange, que fijan las trazas.
* Datos: dGf de JANAF (data/janaf_dgf.json), J/mol.

Sets
   t  'temperaturas' / 600, 700, 800, 900, 1000, 1100 /
   i  'especies'     / CH4, H2O, CO, CO2, H2 /
   e  'elementos'    / C, H, O / ;

Table dGf(t,i) 'J por mol (JANAF C-067, H-064, C-093, C-095; H2 referencia)'
                CH4        H2O         CO        CO2         H2
   600       -22887    -214007    -164486    -395182          0
   700       -12643    -208812    -173518    -395398          0
   800        -2115    -203496    -182497    -395586          0
   900         8616    -198083    -191416    -395748          0
   1000       19492    -192590    -200275    -395886          0
   1100       30472    -187033    -209075    -396001          0 ;

Table A(e,i) 'atomos del elemento e por molecula de i'
        CH4  H2O  CO  CO2  H2
   C     1    0    1   1    0
   H     4    2    0   0    2
   O     0    1    1   2    0 ;

Parameters
   n0(i)  'alimentacion en mol (H2O:CH4 = 4:1)' / CH4 1, H2O 4 /
   b(e)   'moles de atomos en la carga'
   gRT(i) 'dGf entre RT a la T del ciclo'
   yT(t,i)   'fraccion mol'
   lamT(t,e) 'lambda entre RT'
   nT(t)     'moles totales'
   aC(t)     'actividad del carbono, exp(-lambda_C/RT)'
   estado(t,*) 'model status: NLP (2 = Locally Optimal) y CNS (16 = Solved)' ;

Scalars  R 'J por mol K' / 8.314 /,  Tk 'K',  P 'bar' / 1 /,  P0 'bar' / 1 / ;
b(e) = sum(i, A(e,i)*n0(i));

Variables  lnn(i), lnNt, z, lam(e) ;
Equations  balance(e), total, objetivo, potencial(i) ;

balance(e)..  sum(i, A(e,i)*exp(lnn(i))) =e= b(e);
total..       sum(i, exp(lnn(i)))        =e= exp(lnNt);
objetivo..    z =e= sum(i, exp(lnn(i))*(gRT(i) + lnn(i) - lnNt + log(P/P0)));
potencial(i)..  gRT(i) + lnn(i) - lnNt + log(P/P0) + sum(e, A(e,i)*lam(e)) =e= 0;

lnn.lo(i) = -700;

Model gibbs / balance, total, objetivo / ;
Model kkt   / balance, total, potencial / ;
option nlp = ipopt, cns = ipopt, solprint = off;

loop(t,
   Tk = t.val;
   gRT(i) = dGf(t,i)/(R*Tk);
*  Arranque: reparto uniforme en cada T, sin adivinar la respuesta.
   lnn.l(i) = log(sum(e, b(e))/(3*card(i)));
   lnNt.l   = log(sum(i, exp(lnn.l(i))));
   solve gibbs using nlp minimizing z;
   estado(t,'NLP') = gibbs.modelstat;
   lam.l(e) = -balance.m(e);
   solve kkt using cns;
   estado(t,'CNS') = kkt.modelstat;
   yT(t,i)   = exp(lnn.l(i) - lnNt.l);
   lamT(t,e) = lam.l(e);
   nT(t)     = exp(lnNt.l);
   aC(t)     = exp(-lam.l('C'));
);

option yT:6:1:1, lamT:4:1:1;
display yT, lamT, nT, aC, estado;

$ifThen set salida
file fcsv / '%salida%' /;
fcsv.nr = 2; fcsv.nd = 12; fcsv.pw = 32767; fcsv.nz = 0;
put fcsv 'caso,clave,valor' /;
loop((t,i), put t.tl:0 ',y_' i.tl:0 ',' yT(t,i):22:14 /);
loop((t,e), put t.tl:0 ',lam_' e.tl:0 ',' lamT(t,e):22:14 /);
loop(t, put t.tl:0 ',nT,' nT(t):22:14 / t.tl:0 ',estado_NLP,' estado(t,'NLP'):22:14 /
           t.tl:0 ',estado_CNS,' estado(t,'CNS'):22:14 /);
putclose fcsv;
$endIf
