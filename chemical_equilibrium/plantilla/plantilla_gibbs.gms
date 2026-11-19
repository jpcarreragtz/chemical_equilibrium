$offListing
$offSymXRef offSymList
$title PLANTILLA - Minimizacion de Gibbs en GAMS (copiar a ejercicios/sNN_ejN_tema/)
* =================================================================================
* COMO USARLA
*   1) Copia este archivo al folder del ejercicio, ejercicios/sNN_ejN_tema/, junto a su .py.
*   2) Edita SOLO el bloque "DATOS": especies, elementos, dGf, A, n0, T, P y phi.
*      Los datos de la maestra se usan tal cual; si no los da, JANAF (CLAUDE.md).
*   3) Corre:  gams <archivo>.gms           (solo la tabla: agrega lo=0)
*      Ruta completa si gams no esta en el PATH:
*      /Library/Frameworks/GAMS.framework/Versions/Current/Resources/gams <archivo>.gms
*   4) Resultados: en la terminal y en <archivo>_resultados.txt (tabla limpia).
*      El .lst solo trae los display finales.
* METODO (convencion del proyecto, ver CLAUDE.md):
*   variable ln(n_i) (sin limite en 0), NLP que minimiza nG/RT y luego CNS con
*   las condiciones de Lagrange, que fijan las trazas. Solver IPOPT, callado.
* =================================================================================

* ═════════════════════════════ DATOS (edita aqui) ═════════════════════════════════
Sets
   i  'especies'   / CH4, H2O, CO, CO2, H2 /
   e  'elementos'  / C, H, O / ;

Parameters
   dGf(i)  'J por mol a la T del problema (fuente: ...)'
      / CH4 19492, H2O -192590, CO -200275, CO2 -395886, H2 0 /
   n0(i)   'alimentacion en mol' / CH4 1, H2O 4 /
   phi(i)  'coeficiente de fugacidad (1 = gas ideal)' ;
phi(i) = 1;

Table A(e,i) 'atomos del elemento e por molecula de i'
        CH4  H2O  CO  CO2  H2
   C     1    0    1   1    0
   H     4    2    0   0    2
   O     0    1    1   2    0 ;

Scalars  T  'K'   / 1000 /
         P  'bar' / 1 /
         P0 'bar' / 1 /
         R  'J por mol K' / 8.314 / ;
$setGlobal TITULO 'Composicion de equilibrio por minimizacion de Gibbs'
* ══════════════════════════════ FIN DE DATOS ══════════════════════════════════════

Set sv 'solves' / NLP, CNS / ;
Parameters
   b(e)      'moles de atomos en la carga'
   gRT(i)    'dGf entre RT, con phi y P'
   y(i)      'fraccion mol en el equilibrio'
   n(i)      'moles en el equilibrio'
   lamRT(e)  'lambda entre RT'
   estado(sv) 'model status'
   bien(sv)   '1 si el solve salio bien' ;
b(e)   = sum(i, A(e,i)*n0(i));
gRT(i) = dGf(i)/(R*T) + log(phi(i)) + log(P/P0);

* IPOPT sin texto en la terminal ni errores de evaluacion en el .lst
$onEcho > ipopt.opt
print_level 0
sb yes
print_eval_error no
$offEcho

Variables  lnn(i), lnNt, z, lam(e) ;
Equations  balance(e), total, objetivo, potencial(i) ;
balance(e)..    sum(i, A(e,i)*exp(lnn(i))) =e= b(e);
total..         sum(i, exp(lnn(i)))        =e= exp(lnNt);
objetivo..      z =e= sum(i, exp(lnn(i))*(gRT(i) + lnn(i) - lnNt));
potencial(i)..  gRT(i) + lnn(i) - lnNt + sum(e, A(e,i)*lam(e)) =e= 0;

lnn.lo(i) = -700;
lnn.l(i)  = log(sum(e, b(e))/(3*card(i)));
lnNt.l    = log(sum(i, exp(lnn.l(i))));

Model gibbs / balance, total, objetivo / ;
Model kkt   / balance, total, potencial / ;
option nlp = ipopt, cns = ipopt, limrow = 0, limcol = 0, solprint = silent;
gibbs.optfile = 1;  kkt.optfile = 1;

solve gibbs using nlp minimizing z;
estado('NLP') = gibbs.modelstat;
bien('NLP')   = (gibbs.modelstat <= 2 and gibbs.solvestat = 1);
lam.l(e) = -balance.m(e);
solve kkt using cns;
estado('CNS') = kkt.modelstat;
bien('CNS')   = (kkt.modelstat >= 15 and kkt.modelstat <= 17 and kkt.solvestat = 1);

n(i) = exp(lnn.l(i));
y(i) = exp(lnn.l(i) - lnNt.l);
lamRT(e) = lam.l(e);
option y:6:0:1, n:6:0:1;
display y, n, lamRT, estado;

* ------------------------------ REPORTE LIMPIO -----------------------------------
* Se guarda en <archivo>_resultados.txt y se copia a la terminal. Solo si este es
* el archivo que se corre (no en la revision automatica de la extension de VS Code).
$ifThen "%system.incParent%" == ""
file frep / '%system.fp%%system.fn%_resultados.txt' /;
frep.pw = 200;  frep.nz = 0;
put frep;
put / '========================================================================' /
      '%TITULO%' /
      'T = ' T:0:2 ' K, P = ' P:0:4 ' bar   [GAMS, IPOPT]' /
      '========================================================================' /
      '  especie        n_i (mol)              y_i' /;
frep.nr = 2;
loop(i, put '  ' i.tl:<10 n(i):15:5 y(i):17:5 /);
frep.nr = 1;
put '  total     ' sum(i, n(i)):15:5 /;
put / 'λ_e/RT: ';
loop(e, put e.tl:0 ' = ' lamRT(e):0:4 '   ');
frep.nr = 2;
put / / 'Balance de átomos, máx |entra − sale| = '
      smax(e, abs(sum(i, A(e,i)*n(i)) - b(e))):0:2 /;
frep.nr = 1;
put 'Estado: ' sum(sv, bien(sv)):0:0 ' de ' card(sv):0:0 ' solves salieron bien.' /;
loop(sv$(not bien(sv)), put '  No salió bien: ' sv.tl:0 ' (model status ' estado(sv):0:0 ').' /);
putclose frep;
execute 'cat "%system.fp%%system.fn%_resultados.txt"';
$endIf
