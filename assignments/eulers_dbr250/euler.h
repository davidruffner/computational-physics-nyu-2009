#ifndef EULER_H
#define EULER_H



struct Conserved
{
  double rho, E, px;
};
struct Primitive
{
  double pre, e, vx;
};
struct Conserved mixed_to_cons(double pre, double rho, double vx);

struct Conserved prim_to_cons(struct Primitive);
struct Primitive cons_to_prim(struct Conserved);


void Eigenvalues(struct Primitive, double*, double*);//(inprim,lamdap,lamdam)
double Alphavalue(double, double );//(lamdal,lamdar)
double MaxEigenVal(struct Conserved *U, int N); 

struct Conserved flux(struct Primitive);
struct Conserved F_HLL(struct Conserved, struct Conserved);//(Ur,Ul)

#endif // EULER_H
