$title Fischer-Tropsch por minimizacion de Gibbs - Sesion 27 (500 K, 25 bar, 12 escenarios)
* Ejecutar:  gams s27_fischer_tropsch.gms        (opcional: --salida=archivo.csv)
* Si la terminal dice "command not found: gams" (no esta en el PATH), usar la ruta completa:
*    /Library/Frameworks/GAMS.framework/Versions/Current/Resources/gams s27_fischer_tropsch.gms
* Usa IPOPT (corre sin restriccion con la licencia Free Personal de GAMS).
* Variable: ln(n_i), sin limite inferior en 0 (sin ASF, C15 ~ 1e-37).
* Dos etapas por escenario: NLP (minimiza nG/RT) y CNS (Lagrange).
* E6_ASF (diapositiva): n(Ck) = n(C1)*alpha**(k-1), k = 2..15. En el CNS cada
* restriccion lleva su multiplicador mu(k), que entra en el potencial de Ck y de C1.
* phi entra en el logaritmo: gRT(i) + ln(phi_i) + ln(y_i) + ln(P/P0).
* Datos: Sesion 27, diapositiva Datos (data/ft_s27_maestra.json).
* SUPUESTO: los analisis 1 y 2 no dicen que b usar; base H2/CO = 2 (C 1, H 4, O 1).

Sets
   i  'especies' / CH4, C2H6, C3H8, C4H10, C5H12, C6H14, C7H16, C8H18, C9H20,
                   C10H22, C11H24, C12H26, C13H28, C14H30, C15H32, CO2, H2O, H2, CO /
   k(i) 'parafinas C2..C15 (restricciones ASF)'
         / C2H6, C3H8, C4H10, C5H12, C6H14, C7H16, C8H18, C9H20,
           C10H22, C11H24, C12H26, C13H28, C14H30, C15H32 /
   e  'elementos' / C, H, O /
   s  'escenarios' / base_noASF_nophi, base_noASF_phi, base_ASF_nophi, base_ASF_phi,
                     bio_noASF_nophi,  bio_noASF_phi,  bio_ASF_nophi,  bio_ASF_phi,
                     gn_noASF_nophi,   gn_noASF_phi,   gn_ASF_nophi,   gn_ASF_phi / ;

Parameters
   dG(i)  'J por mol a 500 K (maestra)'
      / CH4 -32740, C2H6 -17480, C3H8 -1180, C4H10 15120, C5H12 31420,
        C6H14 47720, C7H16 64020, C8H18 80320, C9H20 96620, C10H22 112920,
        C11H24 129220, C12H26 145520, C13H28 161820, C14H30 178120, C15H32 194420,
        CO2 -394900, H2O -219050, H2 0, CO -155400 /
   phi(i) 'coeficiente de fugacidad (maestra)'
      / CH4 0.99, C2H6 0.95, C3H8 0.89, C4H10 0.82, C5H12 0.74, C6H14 0.66,
        C7H16 0.58, C8H18 0.50, C9H20 0.43, C10H22 0.36, C11H24 0.30,
        C12H26 0.25, C13H28 0.21, C14H30 0.17, C15H32 0.14,
        CO2 0.97, H2O 0.88, H2 1.01, CO 1.00 /
   nc(i)  'carbonos por molecula de parafina'
      / CH4 1, C2H6 2, C3H8 3, C4H10 4, C5H12 5, C6H14 6, C7H16 7, C8H18 8,
        C9H20 9, C10H22 10, C11H24 11, C12H26 12, C13H28 13, C14H30 14, C15H32 15 / ;

Table A(e,i) 'atomos del elemento e por molecula de i'
       CH4 C2H6 C3H8 C4H10 C5H12 C6H14 C7H16 C8H18 C9H20 C10H22 C11H24 C12H26 C13H28 C14H30 C15H32 CO2 H2O H2 CO
   C    1    2    3    4     5     6     7     8     9    10     11     12     13     14     15     1   0   0  1
   H    4    6    8   10    12    14    16    18    20    22     24     26     28     30     32     0   2   2  0
   O    0    0    0    0     0     0     0     0     0     0      0      0      0      0      0     2   1   0  1 ;

Table esc(s,*) 'b (C, H, O), ASF (1 = si) y phi (1 = si)'
                       bC    bH    bO   asf   usaphi
   base_noASF_nophi     1     4     1     0     0
   base_noASF_phi       1     4     1     0     1
   base_ASF_nophi       1     4     1     1     0
   base_ASF_phi         1     4     1     1     1
   bio_noASF_nophi      1     2     1     0     0
   bio_noASF_phi        1     2     1     0     1
   bio_ASF_nophi        1     2     1     1     0
   bio_ASF_phi          1     2     1     1     1
   gn_noASF_nophi       1     5     1     0     0
   gn_noASF_phi         1     5     1     0     1
   gn_ASF_nophi         1     5     1     1     0
   gn_ASF_phi           1     5     1     1     1 ;

Parameters  b(e), gRT(i), nS(s,i) 'moles', yS(s,i) 'fraccion mol', nTS(s),
            estado(s,*) 'model status: NLP (2 = Locally Optimal) y CNS (16 = Solved)' ;
Scalars  R / 8.314 /, Tk / 500 /, P / 25 /, P0 / 1 /, alpha / 0.80 /,
         asfOn / 0 /, phiOn / 0 / ;

Variables  lnn(i), lnNt, z, lam(e), mu(i) ;
Equations  balance(e), total, objetivo, E6_ASF(i), potencial(i) ;

balance(e)..  sum(i, A(e,i)*exp(lnn(i))) =e= b(e);
total..       sum(i, exp(lnn(i)))        =e= exp(lnNt);
objetivo..    z =e= sum(i, exp(lnn(i))*(gRT(i) + log(phi(i))$phiOn + lnn(i) - lnNt + log(P/P0)));
E6_ASF(k)$asfOn..  exp(lnn(k)) =e= exp(lnn('CH4'))*alpha**(nc(k)-1);
potencial(i)..  gRT(i) + log(phi(i))$phiOn + lnn(i) - lnNt + log(P/P0)
                + sum(e, A(e,i)*lam(e))
                + mu(i)$(asfOn and k(i))
                - sum(k$asfOn, mu(k)*alpha**(nc(k)-1))$sameas(i,'CH4') =e= 0;

lnn.lo(i) = -700;

Model gibbs / balance, total, objetivo, E6_ASF / ;
Model kkt   / balance, total, potencial, E6_ASF / ;
option nlp = ipopt, cns = ipopt, solprint = off;

gRT(i) = dG(i)/(R*Tk);
loop(s,
   b('C') = esc(s,'bC');  b('H') = esc(s,'bH');  b('O') = esc(s,'bO');
   asfOn = esc(s,'asf');  phiOn = esc(s,'usaphi');
   lnn.l(i) = log(sum(e, b(e))/(3*card(i)));
   lnNt.l   = log(sum(i, exp(lnn.l(i))));
   mu.l(i)  = 0;
   solve gibbs using nlp minimizing z;
   estado(s,'NLP') = gibbs.modelstat;
   lam.l(e) = -balance.m(e);
   mu.l(k)  = -E6_ASF.m(k);
   solve kkt using cns;
   estado(s,'CNS') = kkt.modelstat;
   nS(s,i)  = exp(lnn.l(i));
   yS(s,i)  = exp(lnn.l(i) - lnNt.l);
   nTS(s)   = exp(lnNt.l);
);

Parameter resumen(s,*) 'n de CH4, CO2, H2O, suma C5-C15 y H2 libre';
resumen(s,'CH4') = nS(s,'CH4');  resumen(s,'CO2') = nS(s,'CO2');
resumen(s,'H2O') = nS(s,'H2O');  resumen(s,'H2') = nS(s,'H2');
resumen(s,'C5_C15') = sum(i$(nc(i) >= 5), nS(s,i));
option resumen:6:1:1;
display resumen, estado;

$ifThen set salida
file fcsv / '%salida%' /;
fcsv.nr = 2; fcsv.nd = 12; fcsv.pw = 32767; fcsv.nz = 0;
put fcsv 'caso,clave,valor' /;
loop((s,i), put s.tl:0 ',y_' i.tl:0 ',' yS(s,i):22:14 / s.tl:0 ',n_' i.tl:0 ',' nS(s,i):22:14 /);
loop(s, put s.tl:0 ',nT,' nTS(s):22:14 / s.tl:0 ',estado_NLP,' estado(s,'NLP'):22:14 /
           s.tl:0 ',estado_CNS,' estado(s,'CNS'):22:14 /);
putclose fcsv;
$endIf
