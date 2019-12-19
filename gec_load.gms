$title CSV2GDX Example 2 - Reading CSV Files with CSV2GDX (CSV2GDX2,SEQ=111)

Set
         tgc 'test get cols'
         gecs 'gecs'
         dir 'directions' /x, y/;
*   i 'canning plants'
*   j 'markets';

alias(gecs, g1, g2, g3, g4);

parameter d(gecs, tgc)
         resX  'resolution width'
         resY  'resolution height'
         timestep 'duration of timestep'
         inverseTimestep;
resX = 1280;
resY = 720;
timestep = 1/30;
inverseTimestep = 1/timestep;

variable finalPosition(g1);
variable finalVelocity(g1);


*$call csv2gdx distance.csv id=d index=1 values=2..lastCol useHeader=y trace=0
*$ifE errorLevel<>0 $abort Problems reading distance.csv!
*$gdxIn distance.gdx
*$load i = dim1
*$load j = dim2


$call csv2gdx testGecsData.csv id=d index=1 values=2..lastCol useHeader=y trace=0
$ifE errorLevel<>0 $abort Problems reading testGecsData.csv!
$gdxIn testGecsData.gdx
$load gecs = dim1
$load tgc = dim2

$load d
$gdxIn


* velocity of a gec at the beginning of the time step
parameter velocity(g1, dir);
velocity(g1, 'x') = d(g1, 'velocityx');
velocity(g1, 'y') = d(g1, 'velocityy');

Set gecPair(g1, g2);

* pairs of valid gecPairs (including duplicates)
Set gecPairPair(g1, g2, g3, g4);

gecPair(g1, g2) $(ord(g1)>ord(g2) and ord(g1)<>20 and ord(g2)<>20)=YES;

gecPairPair(g1, g2, g3, g4) $(gecPair(g1, g2)and gecPair(g3, g4))=YES;

* the pairs (g2, g3) containing g1
Set gecPairsWithGec(g1, g2, g3);
gecPairsWithGec(g1, g2, g3) $((ord(g1)=ord(g2) or ord(g1)=ord(g3)) and ord(g2)>ord(g3))=YES;

* inverseMassScalar for each gecPair
parameter pairInverseMassScalar(g1, g2);
pairInverseMassScalar(g1, g2) $gecPair(g1, g2) = d(g1, 'inverseMass') + d(g2, 'inverseMass');


* _massScalar for each gecPair
parameter pair_massScalar(g1, g2);
* max in place to avoid /0 error for trivial calculations, which should not occur
pair_massScalar(g1, g2) $gecPair(g1, g2) = 1 / max(pairInverseMassScalar(g1, g2), 0.00001);

* _radiiSum for each gecPair
parameter pair_radiiSum(g1, g2);
pair_radiiSum(g1, g2) $gecPair(g1, g2) = d(g1, 'radius') + d(g2, 'radius');

* run GecPair.advance()
* actual distances between the centers in x and y components
parameter centerDifference(g1, g2, dir);
centerDifference(g1, g2, 'x')  $gecPair(g1, g2) =  d(g2, 'positionx') - d(g1, 'positionx') ;
centerDifference(g1, g2, 'y')  $gecPair(g1, g2) =  d(g2, 'positiony') - d(g1, 'positiony') ;

* actual distances between the centers
parameter centerDistance(g1, g2);
centerDistance(g1, g2)  $gecPair(g1, g2) =  sqrt(sum(dir, sqr(centerDifference(g1, g2, dir)))) ;

* unit-scaled distance between centers in x and y components
parameter pair_direction(g1, g2, dir) ;
pair_direction(g1, g2, dir)  $gecPair(g1, g2) = centerDifference(g1, g2, dir) / centerDistance(g1, g2);

parameter pair_distanceComponentOfVelocityCushion(g1, g2);
pair_distanceComponentOfVelocityCushion(g1, g2) $gecPair(g1, g2) =  (centerDistance(g1, g2) -  pair_radiiSum(g1, g2)) * inverseTimestep;


* end of GecPair.advance()

* pairVelocityCushion is the amount of velocity that can be applied to a pair without them colliding.
* if negative: how much you need to fix it in velocity units
* if positive: how much extra velocity
parameter pairVelocityCushion(g1, g2);
pairVelocityCushion(g1, g2) = sum(dir, (velocity(g2, dir) -velocity(g1,dir)) *pair_direction(g1, g2, dir));



* q, h, G, P from cone form documentation
parameter q(g1,g2);
q(g1,g2) $gecPair(g1, g2) = pairVelocityCushion(g1, g2);
parameter h(g1,g2);
h(g1,g2) $gecPair(g1, g2) = pairVelocityCushion(g1, g2);

* P and G are 2d matricies of gecPairs
parameter P(g1,g2, g3, g4);
parameter G(g1,g2, g3, g4);

* diagonal entries in P,G
P(g1,g2, g3, g4) $(gecPairPair(g1, g2, g3, g4) and ord(g1)=ord(g3) and ord(g2)=ord(g4)) = 2 * pairInverseMassScalar(g1, g2);
G(g1,g2, g3, g4) $(gecPairPair(g1, g2, g3, g4) and ord(g1)=ord(g3) and ord(g2)=ord(g4)) = -1 * pairInverseMassScalar(g1, g2);


* pairs of pairs that share a gec (aGec)
* calculate how much impule affects otherP's velocity
P(g1,g2, g3, g4) $(ord(g1)=ord(g3) and ord(g2)<>ord(g4)) =  2 * d(g1, 'inverseMass') * sum(dir, pair_direction(g1, g2, dir)*pair_direction(g3, g4, dir));
P(g1,g2, g3, g4) $(ord(g1)=ord(g4))                      = -2 * d(g1, 'inverseMass') * sum(dir, pair_direction(g1, g2, dir)*pair_direction(g3, g4, dir));

G(g1,g2, g3, g4) $(ord(g1)=ord(g3) and ord(g2)<>ord(g4)) = -1 * d(g1, 'inverseMass') * sum(dir, pair_direction(g1, g2, dir)*pair_direction(g3, g4, dir));
G(g1,g2, g3, g4) $(ord(g1)=ord(g4))                      =      d(g1, 'inverseMass') * sum(dir, pair_direction(g1, g2, dir)*pair_direction(g3, g4, dir));


* _appliedImpulse for the pair (the amount of physical impulse )
* this is x in the cone problem
variable pair_appliedImpulse(g1, g2);
pair_appliedImpulse.lo(g1, g2) = 0;

* cone constraint   Gx <= h
equation
         st(g1, g2)
         getXtpx(g1, g2) 'get the (vector) product of x transpose and P and x'
         getObj;
st(g1, g2) $(gecPair(g1, g2))..  sum((g3, g4) $(gecPair(g1, g2) and gecPair(g3, g4) and (ord(g1)<>ord(g3) or ord(g2)<>ord(g4))), pair_appliedImpulse(g3, g4)) =g= h(g1, g2);

variable xtpx(g1, g2);
getXtpx(g1, g2).. Xtpx(g1, g2)$(gecPair(g1, g2)) =e= sum((g3, g4) $(gecPair(g1, g2) and gecPair(g3, g4) and (ord(g1)<>ord(g3) or ord(g2)<>ord(g4))), pair_appliedImpulse(g1, g2)*P(g1, g2, g3, g4))*pair_appliedImpulse(g1, g2);

variable obj;
getObj.. obj =e= sum((g1, g2)$(gecPair(g1, g2)) ,(0.5*xtpx(g1, g2) )+(q(g1, g2)*pair_appliedImpulse(g1, g2)));

display tgc, gecs, pairVelocityCushion;

model cone /all/;
solve cone using qcp minimizing obj;
display pair_appliedImpulse.l, obj.l;

