$title LendingClub Problem

set
 i  credit grade     / A, B, C, D, E, F, G/
 j  loan category    /L1*L5/
 s  scenarios        /S1*S10000/;


Scalar beta      prob. of not being in the tail               /0.8/;


Table N(i,j)          Number of persons in each category - integer
         L1        L2        L3        L4        L5
A        10        10        10        10        10
B        10        10        10        10        10
C        10        10        10        10        10
D        10        10        10        10        10
E        10        10        10        10        10
F        10        10        10        10        10
G        10        10        10        10        10



;

parameters LA(j) loan_amount
/L1=5018
L2=10473
L3=17546
L4=24990
L5=32970/;

Table intr(i,j)           interest rate each category_did not default
         L1            L2            L3            L4            L5
A        0.0905        0.0913        0.0914        0.0917        0.1015
B        0.1305        0.1320        0.1322        0.1318        0.1323
C        0.1703        0.1712        0.1714        0.1718        0.1756
D        0.2116        0.2138        0.2092        0.2088        0.2117
E        0.2471        0.2453        0.2358        0.2361        0.2322
F        0.2952        0.2814        0.3007        0.2659        0.3086
G        0.3199        0.3701        0.2733        0.3403        0.3663
;


Table PD(i,j)          probability of default
         L1             L2             L3             L4             L5
A        0.05906        0.05519        0.05486        0.04869        0.04073
B        0.11411        0.12463        0.12611        0.11085        0.09326
C        0.18181        0.20421        0.20996        0.19475        0.17072
D        0.23909        0.27811        0.28393        0.27158        0.25844
E        0.30539        0.33930        0.35764        0.32247        0.33952
F        0.38170        0.45714        0.54688        0.42647        0.38462
G        0.40000        0.50588        0.52174        0.50000        0.50000
;

Table MDR(i,j)          mean of return rate for those who default
         L1            L2            L3            L4            L5
A        0.6429        0.6445        0.6427        0.6387        0.6476
B        0.6436        0.6443        0.6432        0.6514        0.6254
C        0.6266        0.6320        0.6238        0.6285        0.6308
D        0.6194        0.6193        0.6138        0.6036        0.5883
E        0.6089        0.6082        0.5790        0.5991        0.5871
F        0.6339        0.6253        0.5594        0.5668        0.5496
G        0.5161        0.5446        0.5505        0.6015        0.6112
;

Table SDDR(i,j)          sd of return rate for those who default

         L1            L2            L3            L4            L5
A        0.2591        0.2556        0.2563        0.2407        0.2492
B        0.2696        0.2624        0.2613        0.2604        0.2463
C        0.2824        0.2723        0.2717        0.2741        0.2837
D        0.2941        0.2936        0.2796        0.2930        0.2821
E        0.3064        0.2958        0.2911        0.3090        0.2786
F        0.3260        0.3280        0.3234        0.3383        0.2837
G        0.3001        0.2936        0.2796        0.5337        0.3817


;


Scalar capital      capital available for borrowing               /1000000/;

parameters  default(i,j,s)      default status
            returnd(i,j,s)      default_percentage returned;
default(i,j,s) = randbinomial(1,(PD(i,j)));
returnd(i,j,s) = normal(MDR(i,j),SDDR(i,j));
returnd(i,j,s) = ifthen(returnd(i,j,s)<0, 0, returnd(i,j,s));

variables
         ppl(i,j) 'number of people for credit grade i, loan amount j'
         cost    'cost minus revenue'
         var    value at Risk
         cvar   conditional value at risk
         z(s)   tail profit;

Integer variables ppl;
positive variable z;
equations   cdef   costdefinition
            cc     capital constraint
            nK(i,j)   number of ppl for each category
            tails(s) calculate the value of the tail
            cvar_eq objective function  ;



cdef..      cost =e= sum((i,j,s), (-default(i,j,s)*returnd(i,j,s)*LA(j)-(1-default(i,j,s))*LA(j)*(1+intr(i,j)))*ppl(i,j))/card(s)+sum((i,j), LA(j)*ppl(i,j));

cc..        sum((i,j), LA(j)*ppl(i,j)) =l= capital;
nK(i,j)..        ppl(i,j) =l= N(i,j);
tails(s)..      z(s)=g=sum((i,j), (-default(i,j,s)*returnd(i,j,s)*LA(j)-(1-default(i,j,s))*LA(j)*(1+intr(i,j)))*ppl(i,j))+sum((i,j), LA(j)*ppl(i,j))-var;
cvar_eq..       cvar =e= var+1/((1-beta)*card(s))*sum(s,z(s));

models chance model        /all/

solve chance minimizing cvar using MIP;

display ppl.l, cost.l,var.l,cvar.l;
