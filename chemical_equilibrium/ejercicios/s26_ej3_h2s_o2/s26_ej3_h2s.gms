$offListing
$offSymXRef offSymList
$title Minimizacion de Gibbs - Sesion 26, Ejercicio 3 (H2S + O2, 6 casos)
* Ejecutar:  gams s26_ej3_h2s.gms        (opcional: --salida=archivo.csv)
* Si la terminal dice "command not found: gams" (no esta en el PATH), usar la ruta completa:
*    /Library/Frameworks/GAMS.framework/Versions/Current/Resources/gams s26_ej3_h2s.gms
* Usa IPOPT (corre sin restriccion con la licencia Free Personal de GAMS).
* Variable: ln(n_i), sin limite inferior en 0 (hay trazas de 1e-68). Dos etapas
* en cada caso: NLP (minimiza nG/RT) y CNS (condiciones de Lagrange).
* Datos: dGf de JANAF (H2S H-080, SO2 O-034, S2 S-012, H2O H-064), J/mol.
* Salida: el .lst solo trae los display finales; la tabla de resultados sale en la
* terminal y en s26_ej3_h2s_resultados.txt. IPOPT no escribe en la terminal.

Sets
   c  'casos (el texto es el nombre del caso en el reporte)'
      / a500   'a 500 K 1 bar 2:1',  a700   'a 700 K 1 bar 2:1',
        b5bar  'b 500 K 5 bar 2:1',  b10bar 'b 500 K 10 bar 2:1',
        c11    'c 500 K 1 bar 1:1',  c12    'c 500 K 1 bar 1:2' /
   i  'especies'  / H2S, O2, H2O, SO2, S2, H2 /
   e  'elementos' / H, O, S /
   sv 'solves de cada caso' / NLP, CNS / ;

Table dGf(c,i) 'J por mol a la T de cada caso (500 K salvo a700)'
                 H2S         O2        H2O        SO2         S2         H2
   a500       -40179          0    -219051    -300871      49019          0
   a700       -44201          0    -208812    -299444      22562          0
   b5bar      -40179          0    -219051    -300871      49019          0
   b10bar     -40179          0    -219051    -300871      49019          0
   c11        -40179          0    -219051    -300871      49019          0
   c12        -40179          0    -219051    -300871      49019          0 ;

Table caso(c,*) 'T (K), P (bar) y alimentacion (mol)'
                   T        P      nH2S      nO2
   a500          500        1         2        1
   a700          700        1         2        1
   b5bar         500        5         2        1
   b10bar        500       10         2        1
   c11           500        1         1        1
   c12           500        1         1        2 ;

Table A(e,i) 'atomos del elemento e por molecula de i'
        H2S  O2  H2O  SO2  S2  H2
   H     2    0   2    0    0   2
   O     0    2   1    2    0   0
   S     1    0   0    1    2   0 ;

Parameters
   n0(i), b(e), gRT(i)
   yC(c,i)    'fraccion mol'
   lamC(c,e)  'lambda entre RT'
   nT(c)      'moles totales'
   estado(c,*) 'model status: NLP (2 = Locally Optimal) y CNS (16 = Solved)'
   bien(c,sv)  '1 si el solve salio bien (NLP: model status 1 o 2; CNS: 15 a 17)' ;

Scalars  R 'J por mol K' / 8.314 /,  Tk, P, P0 / 1 / ;

* IPOPT sin texto en la terminal ni errores de evaluacion en el .lst: las opciones
* van en ipopt.opt, que se escribe aqui mismo en la carpeta de trabajo.
$onEcho > ipopt.opt
print_level 0
sb yes
print_eval_error no
$offEcho

Variables  lnn(i), lnNt, z, lam(e) ;
Equations  balance(e), total, objetivo, potencial(i) ;

balance(e)..  sum(i, A(e,i)*exp(lnn(i))) =e= b(e);
total..       sum(i, exp(lnn(i)))        =e= exp(lnNt);
objetivo..    z =e= sum(i, exp(lnn(i))*(gRT(i) + lnn(i) - lnNt + log(P/P0)));
potencial(i)..  gRT(i) + lnn(i) - lnNt + log(P/P0) + sum(e, A(e,i)*lam(e)) =e= 0;

lnn.lo(i) = -700;

Model gibbs / balance, total, objetivo / ;
Model kkt   / balance, total, potencial / ;
* solprint = silent quita tambien el resumen de cada solve (con off se queda).
option nlp = ipopt, cns = ipopt, limrow = 0, limcol = 0, solprint = silent;
gibbs.optfile = 1;  kkt.optfile = 1;

loop(c,
   Tk = caso(c,'T');  P = caso(c,'P');
   n0(i) = 0;  n0('H2S') = caso(c,'nH2S');  n0('O2') = caso(c,'nO2');
   b(e) = sum(i, A(e,i)*n0(i));
   gRT(i) = dGf(c,i)/(R*Tk);
   lnn.l(i) = log(sum(e, b(e))/(3*card(i)));
   lnNt.l   = log(sum(i, exp(lnn.l(i))));
   solve gibbs using nlp minimizing z;
   estado(c,'NLP') = gibbs.modelstat;
   bien(c,'NLP')   = (gibbs.modelstat <= 2 and gibbs.solvestat = 1);
   lam.l(e) = -balance.m(e);
   solve kkt using cns;
   estado(c,'CNS') = kkt.modelstat;
   bien(c,'CNS')   = (kkt.modelstat >= 15 and kkt.modelstat <= 17 and kkt.solvestat = 1);
   yC(c,i)   = exp(lnn.l(i) - lnNt.l);
   lamC(c,e) = lam.l(e);
   nT(c)     = exp(lnNt.l);
);

option yC:6:1:1, lamC:4:1:1;
display yC, lamC, nT, estado;

$ifThen set salida
file fcsv / '%salida%' /;
fcsv.nr = 2; fcsv.nd = 12; fcsv.pw = 32767; fcsv.nz = 0;
put fcsv 'caso,clave,valor' /;
loop((c,i), put c.tl:0 ',y_' i.tl:0 ',' yC(c,i):22:14 /);
loop((c,e), put c.tl:0 ',lam_' e.tl:0 ',' lamC(c,e):22:14 /);
loop(c, put c.tl:0 ',nT,' nT(c):22:14 / c.tl:0 ',estado_NLP,' estado(c,'NLP'):22:14 /
           c.tl:0 ',estado_CNS,' estado(c,'CNS'):22:14 /);
putclose fcsv;
$endIf

* ---------------------------------------------------------------------------------
* Reporte de formato fijo, con los titulos de s26_ej3_h2s_o2_gibbs.py: se guarda junto a
* este archivo y se copia a la terminal. Solo cuando este es el archivo que se corre
* (no cuando otro .gms lo incluye, como la revision automatica de VS Code).
$ifThen "%system.incParent%" == ""
Scalars  nTot 'solves corridos', nBien 'solves que salieron bien' ;
nTot = card(c)*card(sv);  nBien = sum((c,sv), bien(c,sv));

file frep / '%system.fp%%system.fn%_resultados.txt' /;
frep.pw = 200;  frep.nz = 0;
put frep;
put / '==============================================================================' /
      'Composición de equilibrio (y_i) y moles totales   [GAMS, IPOPT]' /
      '==============================================================================' /
      'caso                        H2S         O2        H2O        SO2         S2         H2     n_T' /;
loop(c,
   put c.te(c):<20;
   frep.nr = 2;  loop(i, put yC(c,i):11:4);
   frep.nr = 1;  put nT(c):8:4 /;
);
put / '==============================================================================' /
      'Azufre: fracción de S alimentado que sale como S2 y como SO2' /
      '==============================================================================' /;
loop(c,
   put c.te(c):<20 ' S2: '               (100*2*yC(c,'S2')*nT(c)/caso(c,'nH2S')):6:2
       '%   SO2: '                       (100*yC(c,'SO2')*nT(c)/caso(c,'nH2S')):6:2
       '%   H2S sin reaccionar: '        (100*yC(c,'H2S')*nT(c)/caso(c,'nH2S')):6:2 '%' /;
);
put / 'Estado: ' nBien:0:0 ' de ' nTot:0:0 ' solves salieron bien.' /;
loop((c,sv)$(not bien(c,sv)),
   put '  No salió bien: ' sv.tl:0 ', caso ' c.te(c):0 ' (model status ' estado(c,sv):0:0 ').' /;
);
putclose frep;
execute 'cat "%system.fp%%system.fn%_resultados.txt"';
$endIf
