$offListing
$offSymXRef offSymList
$title Minimizacion de Gibbs - Sesion 25, Ejercicio 2 (steam cracker, 1000 K, 1 bar)
* Ejecutar:  gams s25_ej2_gibbs.gms
* Usa IPOPT, que corre sin restriccion con la licencia Free Personal de GAMS.
* Variable: ln(n_i), no n_i, para que n_i nunca toque cero.
* Dos etapas:
*   1) NLP: minimiza nG/RT. Las especies mayoritarias salen bien, pero las trazas
*      no: casi no pesan en G y el solver para antes aunque diga "Locally Optimal".
*   2) CNS: resuelve las condiciones de Lagrange, que si fijan las trazas.
* Salida: el .lst solo trae los display finales; la tabla de resultados sale en la
* terminal y en s25_ej2_gibbs_resultados.txt. IPOPT no escribe en la terminal.

Sets
   i  'especies'   / CH4, C2H4, C2H2, CO2, CO, O2, H2, H2O, C2H6 /
   e  'elementos'  / C, H, O / ;

Parameters
   dGf(i)  'J por mol a 1000 K (JANAF; C2H6 no esta en JANAF: SVA + Cp)'
      / CH4 19492, C2H4 119122, C2H2 169607, CO2 -395886, CO -200275
        O2 0, H2 0, H2O -192590, C2H6 110873 /
   n0(i)   'alimentacion en mol' / C2H6 1, H2O 4 /
   gRT(i)  'dGf entre RT'
   y(i)    'fraccion mol en el equilibrio'
   logy(i) 'log10 de y, para leer las trazas' ;

Scalars  R 'J por mol K' / 8.314 /,  T 'K' / 1000 /,  P 'bar' / 1 /,  P0 'bar' / 1 / ;

Table A(e,i) 'atomos del elemento e por molecula de i'
       CH4  C2H4  C2H2  CO2  CO  O2  H2  H2O  C2H6
   C    1    2     2     1   1   0   0   0    2
   H    4    4     2     0   0   0   2   2    6
   O    0    0     0     2   1   2   0   1    0 ;

Parameter b(e) 'moles de atomos en la carga' ;
b(e)   = sum(i, A(e,i)*n0(i));
gRT(i) = dGf(i)/(R*T);

* IPOPT sin texto en la terminal ni errores de evaluacion en el .lst: las opciones
* van en ipopt.opt, que se escribe aqui mismo en la carpeta de trabajo.
$onEcho > ipopt.opt
print_level 0
sb yes
print_eval_error no
$offEcho

Variables  lnn(i) 'ln de n_i',  lnNt 'ln de los moles totales',  z 'nG entre RT',
           lam(e) 'multiplicador de Lagrange entre RT' ;
Equations  balance(e), total, objetivo, potencial(i) ;

balance(e)..  sum(i, A(e,i)*exp(lnn(i))) =e= b(e);
total..       sum(i, exp(lnn(i)))        =e= exp(lnNt);
objetivo..    z =e= sum(i, exp(lnn(i))*(gRT(i) + lnn(i) - lnNt + log(P/P0)));

potencial(i)..  gRT(i) + lnn(i) - lnNt + log(P/P0) + sum(e, A(e,i)*lam(e)) =e= 0;

* Arranque: reparto uniforme. No hay que adivinar la respuesta.
lnn.l(i) = log(sum(e, b(e))/(3*card(i)));
lnNt.l   = log(sum(i, exp(lnn.l(i))));
lnn.lo(i) = -100;

Model gibbs / balance, total, objetivo / ;
Model kkt   / balance, total, potencial / ;
* solprint = silent quita tambien el resumen de cada solve (con off se queda).
option nlp = ipopt, cns = ipopt, limrow = 0, limcol = 0, solprint = silent;
gibbs.optfile = 1;  kkt.optfile = 1;

Scalars  nTot 'solves corridos' / 2 /, bienNLP, bienCNS ;
solve gibbs using nlp minimizing z;
bienNLP = (gibbs.modelstat <= 2 and gibbs.solvestat = 1);
lam.l(e) = -balance.m(e);
solve kkt using cns;
bienCNS = (kkt.modelstat >= 15 and kkt.modelstat <= 17 and kkt.solvestat = 1);

y(i) = exp(lnn.l(i) - lnNt.l);
logy(i) = log10(y(i));
option y:6:0:1, logy:4:0:1;
display y, logy, lam.l;

* ---------------------------------------------------------------------------------
* Reporte de formato fijo, con los titulos de s25_ej2_steam_cracker_gibbs.py: se guarda junto a
* este archivo y se copia a la terminal. Solo cuando este es el archivo que se corre
* (no cuando otro .gms lo incluye, como la revision automatica de VS Code).
$ifThen "%system.incParent%" == ""
Alias (i,j);
Parameters  nEq(i) 'mol en el equilibrio', lugar(i) 'posicion de mayor a menor n_i' ;
nEq(i)   = exp(lnn.l(i));
lugar(i) = 1 + sum(j$(nEq(j) > nEq(i)), 1);

file frep / '%system.fp%%system.fn%_resultados.txt' /;
frep.pw = 200;  frep.nz = 0;
put frep;
put / '  T = 1000 K, P = 1 bar   [GAMS, IPOPT]' /
      '  especie       n_i (mol)           y_i' /;
frep.nr = 2;
loop(j, loop(i$(lugar(i) = ord(j)),
   put '  ' i.tl:<8 nEq(i):15:5 y(i):14:5 /;
));
put '  total   ' exp(lnNt.l):15:5 /;
frep.nr = 1;
put '  λ/RT: C = ' lam.l('C'):0:4 ', H = ' lam.l('H'):0:4 ', O = ' lam.l('O'):0:4 /;
put / 'Estado: ' (bienNLP + bienCNS):0:0 ' de ' nTot:0:0 ' solves salieron bien.' /;
if(not bienNLP, put '  No salió bien: NLP (model status ' gibbs.modelstat:0:0 ').' /;);
if(not bienCNS, put '  No salió bien: CNS (model status ' kkt.modelstat:0:0 ').' /;);
putclose frep;
execute 'cat "%system.fp%%system.fn%_resultados.txt"';
$endIf
