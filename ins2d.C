// =====================================================================
//
// Incompressible Navier Stokes 2D
//
// Finite difference solver
// Dirichlet conditions
// First order in time, second in space
//
//
// Solve with A++ arrays & openMP
//
//
// outputs snowflake positions to a Matlab file to be read
//
// =====================================================================

#include "A++.h"
#include <omp.h>

// define some types
typedef double Real;
typedef doubleSerialArray RealArray;
typedef intSerialArray IntegerArray;

#include <math.h>
#include <float.h>
#include <limits.h>
#define REAL_EPSILON DBL_EPSILON
#define REAL_MIN DBL_MIN

// getCPU() : Return the current wall-clock time in seconds
 #include "getCPU.h"

// include commands tp parse command line arguments
#include "parseCommand.h"












// ------------------------------------------------------------------------
// Return the max-norm residual (for pressure solve)
// ------------------------------------------------------------------------
Real getMaxResidual( doubleArray & pn, doubleArray & f, IntegerArray & gridIndexRange, Real *dx, int numThreads )
{
  
  const Real dxSq = dx[0]*dx[0];
  const Real dySq = dx[1]*dx[1];
  
  const int n1a = gridIndexRange(0,0), n1b = gridIndexRange(1,0);
  const int n2a = gridIndexRange(0,1), n2b = gridIndexRange(1,1);

  int i,i1,i2;
  
  
  Real maxRes=0.;

  #pragma omp parallel for default(shared) private(i1,i2) num_threads(numThreads) reduction(max:maxRes)
  for( i2=n2a+1; i2<=n2b-1; i2++ )
  for( i1=n1a+1; i1<=n1b-1; i1++ )
  {
    Real res = f(i1,i2) + ( (pn(i1+1,i2) -2.*pn(i1,i2) + pn(i1-1,i2))/dxSq
                            + (pn(i1,i2+1) -2.*pn(i1,i2) + pn(i1,i2-1))/dySq );
    maxRes = max( maxRes,abs(res) );
  }


  return maxRes;
  
}















// ======================================================================================
// main function
// ======================================================================================
int main(int argc, char *argv[])
{
  
  // =================== setting and parsing parameters  =================== 
  //defining pi
  const Real pi = 4.*atan2(1.,1.);

  ios::sync_with_stdio(); // Synchronize C++ and C I/O subsystems
  Index::setBoundsCheck(on); // Turn on A++ array bounds checking

  printf("Usage: ins2d -nx=<i> -tFinal=<f> -nf=<i> -filename=<s> -trackSnow=<i> -numThreads=<i>\n");

  const int numberOfDimensions=2;

  Real xa=0., xb=1.; // domain is [xa,xb] X [ya,yb]
  Real ya=0., yb=1.;

  Real nu=.1;
  Real c=1.;
  Real cdv=1.;
  Real cfl=.9;

  int i1,i2;

  int computeErrors=1;
  string filename = "snowPosition.m";

   // solution stuff
  #define polyms 0
  #define forcedbcs 1
  
  #ifndef solution
    #define solution polyms
  //#define solution forcedbcs
  #endif

  // for pressure solve
  Real tol=1.e-12; 
  int  maxIterations=500;
  int intervalToCheckResidual=10;


  // snow options
  int trackSnow=0;
  int nf=110; // number of flakes
  int numPerLine=10;
  Real ydensity=.1;
  Real gravpull = -1; // falling speed
  Real windfactor=1;
  

  // defaults
  int debug= 0;
  Real tFinal=.5;
  int nx=50, ny=nx;
  int numThreads=1;
  


  enum BoundaryConditionsEnum
    {
      periodic=-1,
      dirichlet=1,
      neumann=2
    };
  

  // taking in parameter settings
  string line;
  bool echo = true; // do or do not echo in parseCommand
  for( int i=1; i<argc; i++ )
  {
    line=argv[i];
    
    if(      parseCommand( line,"-nx=",nx,echo)                         ){ny=nx;}
    else if( parseCommand( line,"-debug=",debug,echo)                   ){}
    else if( parseCommand( line,"-tFinal=",tFinal,echo)                 ){}
    else if( parseCommand( line,"-numThreads=",numThreads,echo)         ){}
    else if( parseCommand( line,"-nf=",nf,echo)                         ){}
    else if( parseCommand( line,"-trackSnow=",trackSnow,echo)           ){}
    else if( parseCommand( line,"-filename=",filename,echo)             ){}
  }





















  
  // =================== grid set up and setting bcs ===================
  // indexing
  const int numGhost= 1;
  const int n1a     = 0;
  const int n1b     = n1a + nx;
  const int nd1a    = n1a-numGhost;
  const int nd1b    = n1b+numGhost;
  const int nd1     = nd1b-nd1a+1;

  const int n2a     = 0;
  const int n2b     = n2a + ny;
  const int nd2a    = n2a-numGhost;
  const int nd2b    = n2b+numGhost;
  const int nd2     = nd2b-nd2a+1;

  IntegerArray gridIndexRange(2,numberOfDimensions);
  IntegerArray dimension(2,numberOfDimensions);
  IntegerArray boundaryCondition(2,numberOfDimensions);

  gridIndexRange(0,0)=n1a; gridIndexRange(1,0)=n1b;
  gridIndexRange(0,1)=n2a; gridIndexRange(1,1)=n2b;

  dimension(0,0)=nd1a; dimension(1,0)=nd1b;
  dimension(0,1)=nd2a; dimension(1,1)=nd2b;



  boundaryCondition(0,0)=dirichlet; // left
  boundaryCondition(1,0)=dirichlet; // right
  boundaryCondition(0,1)=dirichlet; // bottom
  boundaryCondition(1,1)=dirichlet; // top


  // Grid points
  Range Rx(nd1a,nd1b), Ry(nd2a,nd2b);
  RealArray x(Rx,Ry,2);

 
   Real dx[2];
   dx[0] = (xb-xa)/nx;
   dx[1] = (yb-ya)/ny;

   #pragma omp parallel for default(shared) private(i1,i2) num_threads(numThreads)
   for( i2=nd2a; i2<=nd2b; i2++ )
     for( i1=nd1a; i1<=nd1b; i1++ )
     {
       x(i1,i2,0) = xa + (i1-n1a)*dx[0];
       x(i1,i2,1) = ya + (i2-n2a)*dx[1];
     }











  





   
 
 

   // =================== some solution stuff  ===================  
   // polynomial manufactured solution
   #if solution == polyms
     string solutionName = "polyManufactured";
     Real ct0=1., ct1=.5;
     Real cu0=1., cux=1., cuy=.5, cuxx=1., cuxy=2., cuyy=1.;
     Real cv0=.5, cvx=.8, cvy=-cux, cvxx=1., cvxy=-2*cuxx, cvyy=-.5*cuxy;
     Real cp0=1., cpx=.5, cpy=.25, cpxx=1., cpxy=-.5, cpyy = 1.5;
     
     // exact solutions
     #define ue(x,y,t) windfactor*(ct0 + ct1*(t))*( cu0 + cux*(x) + cuy*(y) + cuxx*(x)*(x) + cuxy*(x)*(y) + cuyy*(y)*(y) )
     #define ve(x,y,t) windfactor*(ct0 + ct1*(t))*( cv0 + cvx*(x) + cvy*(y) + cvxx*(x)*(x) + cvxy*(x)*(y) + cvyy*(y)*(y) )
     #define pe(x,y,t) windfactor*(ct0 + ct1*(t))*( cp0 + cpx*(x) + cpy*(y) + cpxx*(x)*(x) + cpxy*(x)*(y) + cpyy*(y)*(y) )
     
     // derivatives
     #define uet(x,y,t) windfactor*(ct1)*( cu0 + cux*(x) + cuy*(y) + cuxx*(x)*(x) + cuxy*(x)*(y) + cuyy*(y)*(y) )
     #define uex(x,y,t) windfactor*(ct0 + ct1*(t))*( cux + 2*cuxx*(x) + cuxy*(y) )
     #define uey(x,y,t) windfactor*(ct0 + ct1*(t))*( cuy + 2*cuyy*(y) + cuxy*(x) )
     #define uexx(x,y,t) windfactor*(ct0 + ct1*(t))*( 2*cuxx )
     #define ueyy(x,y,t) windfactor*(ct0 + ct1*(t))*( 2*cuyy )
     
     #define vet(x,y,t) windfactor*(ct1)*( cv0 + cvx*(x) + cvy*(y) + cvxx*(x)*(x) + cvxy*(x)*(y) + cvyy*(y)*(y) )
     #define vex(x,y,t) windfactor*(ct0 + ct1*(t))*( cvx + 2*cvxx*(x) + cvxy*(y) )
     #define vey(x,y,t) windfactor*(ct0 + ct1*(t))*( cvy + 2*cvyy*(y) + cvxy*(x) )
     #define vexx(x,y,t) windfactor*(ct0 + ct1*(t))*( 2*cvxx )
     #define veyy(x,y,t) windfactor*(ct0 + ct1*(t))*( 2*cvyy )
     
     #define pex(x,y,t) windfactor*(ct0 + ct1*(t))*( cpx + 2*cpxx*(x) + cpxy*(y) )
     #define pey(x,y,t) windfactor*(ct0 + ct1*(t))*( cpy + 2*cpyy*(y) + cpxy*(x) )
     #define pexx(x,y,t) windfactor*(ct0 + ct1*(t))*( 2*cpxx )
     #define peyy(x,y,t) windfactor*(ct0 + ct1*(t))*( 2*cpyy )


     // initial conditions
     #define uInit(x,y) ue(x,y,0)
     #define vInit(x,y) ve(x,y,0)
     #define pInit(x,y) pe(x,y,0)

   #elif solution == forcedbcs
     string solutionName = "forcedBoundaries";
     computeErrors=0;
     
     // exact solutions for boundaries
     #define ue(x,y,t) windfactor*sin(4*pi*t)
     #define ve(x,y,t) 0
     #define pe(x,y,t) 0

     // initial conditions
     #define uInit(x,y) 0
     #define vInit(x,y) 0
     #define pInit(x,y) 0
   #endif












    // =================== setting time step ===================
   // Choose dt: assume stability region is an ellipse:
   // (reLambda*dt/aStab)^2 + (imLambda*dt/bStab)^2 =1 
   Real aStab=-2., bStab=1.;  

   // define time-stepping eigenvalue
   Real reLambda, imLambda, maxu0, maxv0, dt;
   int  Nt;
   

   maxu0 = 0;
   maxv0 = 0;  
   for( int i2=nd2a; i2<=nd2b; i2++ )
     for( int i1=nd1a; i1<=nd1b; i1++ )
     {
       maxu0 = max( maxu0,fabs(ue(x(i1,i2,0),x(i1,i2,1),0)) );
       maxv0 = max( maxv0,fabs(ve(x(i1,i2,0),x(i1,i2,1),0)) );
     }
   
   reLambda = -4*nu*( 1/(dx[0]*dx[0]) + 1/(dx[1]*dx[1]) ); 
   imLambda = maxu0/dx[0] + maxv0/dx[1];  

   dt = cfl/sqrt( (reLambda/aStab)*(reLambda/aStab) + (imLambda/bStab)*(imLambda/bStab) ); 
   Nt = round(tFinal/dt);      // number of time-steps 
   dt = tFinal/Nt;         // adjust dt to reach tFinal exactly 






   
   // forcing function has to be defined after dt is set
   #define delta(dt) ( (cdv)/(dt) )

   #if solution == polyms
     #define fu(x,y,t) ( uet(x,y,t)+(ue(x,y,t))*(uex(x,y,t))+(ve(x,y,t))*(uey(x,y,t))+pex(x,y,t)-nu*(uexx(x,y,t)+ueyy(x,y,t)) )
     #define fv(x,y,t) ( vet(x,y,t)+(ue(x,y,t))*(vex(x,y,t))+(ve(x,y,t))*(vey(x,y,t))+pey(x,y,t)-nu*(vexx(x,y,t)+veyy(x,y,t)) )
     #define fp(x,y,t) ( pexx(x,y,t)+peyy(x,y,t)+(uex(x,y,t))*(uex(x,y,t))+2*(uey(x,y,t))*(vex(x,y,t))+(vey(x,y,t))*(vey(x,y,t))-delta(dt)*(uex(x,y,t)+vey(x,y,t)) )
   
   #elif solution == forcedbcs
     #define fu(x,y,t) 0
     #define fv(x,y,t) 0
     #define fp(x,y,t) 0
   #endif

 
 



















  

  
  // =================== initializing arrays  ===================
   // store two time levels
   RealArray ua[2];
   ua[0].redim(Rx,Ry); ua[0]=0.;
   ua[1].redim(Rx,Ry); ua[1]=0.;

   RealArray va[2];
   va[0].redim(Rx,Ry); va[0]=0.;
   va[1].redim(Rx,Ry); va[1]=0.;

   RealArray pa[2];
   pa[0].redim(Rx,Ry); pa[0]=0.;
   pa[1].redim(Rx,Ry); pa[1]=0.;

   RealArray prhs(Rx,Ry);
   prhs=0.;
   
   RealArray flakeLocation(nf,2);
   

   // initial conditions
   RealArray & u0 = ua[0];
   RealArray & v0 = va[0];
   RealArray & p0 = pa[0];
   
   Real t=0.;

   #pragma omp parallel for default(shared) private(i1,i2) num_threads(numThreads)
   for( i2=nd2a; i2<=nd2b; i2++ )
     for( i1=nd1a; i1<=nd1b; i1++ )
     {
       u0(i1,i2)= uInit(x(i1,i2,0),x(i1,i2,1));
       v0(i1,i2)= vInit(x(i1,i2,0),x(i1,i2,1));
       p0(i1,i2)= pInit(x(i1,i2,0),x(i1,i2,1));
     }





   // apply bcs on ics 
   for( int axis=0; axis<numberOfDimensions; axis++ )
       for( int side=0; side<=1; side++ )
       {
         int is = 1-2*side; // is=1 on left, is=-1 on right
         if( boundaryCondition(side,axis)==dirichlet )
         {
           if( axis==0 )
           { // left or right side
             int i1  = gridIndexRange(side,axis);
             int i1g = i1 - is; // index of ghost point
             #pragma omp parallel for default(shared) private(i2) num_threads(numThreads)
             for( int i2=nd2a; i2<=nd2b; i2++ )
             {
               u0(i1,i2)  = ue(x(i1,i2,0),x(i1,i2,1),0);
               u0(i1g,i2) = ue(x(i1g,i2,0),x(i1g,i2,1),0);//3.*un(i1,i2) - 3.*un(i1+is,i2) + un(i1+2*is,i2); // extrap ghost

               v0(i1,i2)  = ve(x(i1,i2,0),x(i1,i2,1),0);
               v0(i1g,i2) = ve(x(i1g,i2,0),x(i1g,i2,1),0);//3.*vn(i1,i2) - 3.*vn(i1+is,i2) + vn(i1+2*is,i2);

               p0(i1,i2)  = pe(x(i1,i2,0),x(i1,i2,1),0);
               p0(i1g,i2) = pe(x(i1g,i2,0),x(i1g,i2,1),0);
             }
           }
           else
           { // bottom or top
             int i2  = gridIndexRange(side,axis);
             int i2g = i2 - is; // index of ghost point
             #pragma omp parallel for default(shared) private(i2) num_threads(numThreads)
             for( int i1=nd1a; i1<=nd1b; i1++ )
             {
               u0(i1,i2)  = ue(x(i1,i2,0),x(i1,i2,1),0);
               u0(i1,i2g) = ue(x(i1,i2g,0),x(i1,i2g,1),0);//3.*un(i1,i2) - 3.*u(i1,i2+is) + un(i1,i2+2*is); // extrap ghost

               v0(i1,i2)  = ve(x(i1,i2,0),x(i1,i2,1),0);
               v0(i1,i2g) = ve(x(i1,i2g,0),x(i1,i2g,1),0);//3.*vn(i1,i2) - 3.*vn(i1,i2+is) + vn(i1,i2+2*is);

               p0(i1,i2)  = pe(x(i1,i2,0),x(i1,i2,1),0);
               p0(i1,i2g) = pe(x(i1,i2g,0),x(i1,i2g,1),0);
             }
           }
         }
         else
         {
           printf("ERROR: unknown boundaryCondition=%d\n",boundaryCondition(side,axis));
           abort();
         }
         
       } 




   
 
   int cur=0;
   int i,n;


   // initial snowflake positions
   if( trackSnow > 0 )
   {
     Real xloc, yloc;
     #pragma omp parallel for default(shared) private(i,xloc,yloc) num_threads(numThreads)
     for( i=0; i<nf; i++)
     {
       xloc = double((i)%numPerLine)*(xb-xa)/double(numPerLine);
       yloc = ydensity*(1 + floor(double(i)/double(numPerLine)));
       
       flakeLocation(i,0) = xloc;
       flakeLocation(i,1) = yloc;
     }
   }




   













   
   // =================== some print statements  ===================
   printf("----- Solve the Incompressible Navier Stokes Equations in two dimensions ------\n");
   printf("      Forward Euler, Dirichelt bcs, solution=%s\n",solutionName.c_str());
   printf("      nx=%d, ny=%d, tFinal=%6.2f",nx,ny,tFinal);
   if( trackSnow > 0 )
     printf(", nf=%d \n",nf);
   else
     printf("\n");


   // writing into files
   FILE *snowfile = NULL;   
   if( trackSnow > 0 )
   {
     snowfile = fopen(filename.c_str(),"w");
     
     fprintf(snowfile,"%% File written by ins2d.C\n");
     fprintf(snowfile,"xa=%g; xb=%g; ya=%g; yb=%g; tf=%g; nf=%d; solution=\'%s\';\n",
             xa,xb,ya,yb,tFinal,nf,solutionName.c_str());


     // writing time vector
     fprintf(snowfile,"t=[0");
     for( int i=1; i<=Nt; i++ )
       fprintf(snowfile,", %g",i*dt);
     fprintf(snowfile,"];\n\n");

     
     // writing locations
     #pragma omp parallel for default(shared) private(i) num_threads(numThreads)
     for( int i=0; i<nf; i++ )
       fprintf(snowfile,"sfLocation(%d,%s,%d)=[%g, %g];\n",i+1,":",1,flakeLocation(i,0),flakeLocation(i,1));
     
   }
















  // =================== time stepping loop ===================
   Real cpu1 = getCPU();
   for( n=0; n<Nt; n++ )
   //for( n=0; n<1; n++ )
   {
     t = n*dt; // cur time

     int next = (cur+1) % 2;

     // making things look nicer
     doubleArray & u  = ua[cur];
     doubleArray & un = ua[next];

     const double *u_p = u.getDataPointer();
     double *un_p = un.getDataPointer();
     #define U(i1,i2) u_p[(i1-nd1a)+nd1*(i2-nd2a)]
     #define UN(i1,i2) un_p[(i1-nd1a)+nd1*(i2-nd2a)]

     doubleArray & v  = va[cur];
     doubleArray & vn = va[next];

     const double *v_p = v.getDataPointer();
     double *vn_p = vn.getDataPointer();
     #define V(i1,i2) v_p[(i1-nd1a)+nd1*(i2-nd2a)]
     #define VN(i1,i2) vn_p[(i1-nd1a)+nd1*(i2-nd2a)]

     doubleArray & p  = pa[cur];
     doubleArray & pn = pa[next];
     

     const double *p_p = p.getDataPointer();
     #define P(i1,i2) p_p[(i1-nd1a)+nd1*(i2-nd2a)]

 


     // --- velocity solve ---
     Real up, ux, uy, uxx, uyy, ku;
     Real vp, vx, vy, vxx, vyy, kv;
     Real px, py;

     #pragma omp parallel for default(shared) private(i1,i2,up,ux,uy,uxx,uyy,ku,vp,vx,vy,vxx,vyy,kv,px,py) num_threads(numThreads)
     for( i2=n2a+1; i2<n2b; i2++ )
       for( i1=n1a+1; i1<n1b; i1++ )
       {
         // derivatives
         up  = U(i1,i2);
         ux  = (U(i1+1,i2) - U(i1-1,i2))/(2*dx[0]);
         uy  = (U(i1,i2+1) - U(i1,i2-1))/(2*dx[1]);
         uxx = (U(i1+1,i2) -2*U(i1,i2) + U(i1-1,i2))/(dx[0]*dx[0]);
         uyy = (U(i1,i2+1) -2*U(i1,i2) + U(i1,i2-1))/(dx[1]*dx[1]);
         
         vp  = V(i1,i2);
         vx  = (V(i1+1,i2) - V(i1-1,i2))/(2*dx[0]);
         vy  = (V(i1,i2+1) - V(i1,i2-1))/(2*dx[1]);
         vxx = (V(i1+1,i2) -2*V(i1,i2) + V(i1-1,i2))/(dx[0]*dx[0]);
         vyy = (V(i1,i2+1) -2*V(i1,i2) + V(i1,i2-1))/(dx[1]*dx[1]);
         
         px  = (P(i1+1,i2) - P(i1-1,i2))/(2*dx[0]);
         py  = (P(i1,i2+1) - P(i1,i2-1))/(2*dx[1]);


         // spatial stuff
         ku = -(up)*(ux) - (vp)*(uy) - px + nu*( uxx + uyy ) + fu(x(i1,i2,0),x(i1,i2,1),t);
         kv = -(up)*(vx) - (vp)*(vy) - py + nu*( vxx + vyy ) + fv(x(i1,i2,0),x(i1,i2,1),t);

         // time update
         UN(i1,i2) = U(i1,i2) + dt*(ku);
         //UN(i1,i2) = ue(x(i1,i2,0),x(i1,i2,1),t+dt);
         VN(i1,i2) = V(i1,i2) + dt*(kv);
         //VN(i1,i2) = ve(x(i1,i2,0),x(i1,i2,1),t+dt);
         
       }

     // --- boundary conditions on u, v---
     for( int axis=0; axis<numberOfDimensions; axis++ )
       for( int side=0; side<=1; side++ )
       {
         int is = 1-2*side; // is=1 on left, is=-1 on right
         if( boundaryCondition(side,axis)==dirichlet )
         {
           if( axis==0 )
           { // left or right side
             i1  = gridIndexRange(side,axis);
             int i1g = i1 - is; // index of ghost point

             #pragma omp parallel for default(shared) private(i2) num_threads(numThreads)
             for( i2=nd2a; i2<=nd2b; i2++ )
             {
               un(i1,i2)  = ue(x(i1,i2,0),x(i1,i2,1),t+dt);
               un(i1g,i2) = ue(x(i1g,i2,0),x(i1g,i2,1),t+dt);//3.*un(i1,i2) - 3.*un(i1+is,i2) + un(i1+2*is,i2); // extrap ghost

               vn(i1,i2)  = ve(x(i1,i2,0),x(i1,i2,1),t+dt);
               vn(i1g,i2) = ve(x(i1g,i2,0),x(i1g,i2,1),t+dt);//3.*vn(i1,i2) - 3.*vn(i1+is,i2) + vn(i1+2*is,i2);
             }
           }
           else
           { // bottom or top
             i2  = gridIndexRange(side,axis);
             int i2g = i2 - is; // index of ghost point

             #pragma omp parallel for default(shared) private(i1) num_threads(numThreads)
             for( i1=nd1a; i1<=nd1b; i1++ )
             {
               un(i1,i2)  = ue(x(i1,i2,0),x(i1,i2,1),t+dt);
               un(i1,i2g) = ue(x(i1,i2g,0),x(i1,i2g,1),t+dt);//3.*un(i1,i2) - 3.*u(i1,i2+is) + un(i1,i2+2*is); // extrap ghost

               vn(i1,i2)  = ve(x(i1,i2,0),x(i1,i2,1),t+dt);
               vn(i1,i2g) = ve(x(i1,i2g,0),x(i1,i2g,1),t+dt);//3.*vn(i1,i2) - 3.*vn(i1,i2+is) + vn(i1,i2+2*is);
             }
           }
         }
         else
         {
           printf("ERROR: unknown boundaryCondition=%d\n",boundaryCondition(side,axis));
           abort();
         }
         
       } // end for axis





     


     // --- pressure solve ---
     // fill in forcing, one point wider than cg stentil (not sure about this)
     Real *prhs_p = prhs.getDataPointer();
     #define PRHS(i1,i2) prhs_p[(i1-nd1a)+nd1*(i2-nd2a)]

     #pragma omp parallel for default(shared) private(i1,i2,ux,uy,vx,vy) num_threads(numThreads)
     for( i2=n2a; i2<=n2b; i2++ )
       for( i1=n1a; i1<=n1b; i1++ )
       {
         // derivatives
         ux  = (UN(i1+1,i2) - UN(i1-1,i2))/(2*dx[0]);
         uy  = (UN(i1,i2+1) - UN(i1,i2-1))/(2*dx[1]);
         
         vx  = (VN(i1+1,i2) - VN(i1-1,i2))/(2*dx[0]);
         vy  = (VN(i1,i2+1) - VN(i1,i2-1))/(2*dx[1]);

         PRHS(i1,i2) = -((ux)*(ux) + 2*(uy)*(vx) + (vy)*(vy)) + delta(dt)*(ux + vy) + fp(x(i1,i2,0),x(i1,i2,1),t+dt);
         PRHS(i1,i2) = -PRHS(i1,i2);
       }




     // --- boundary conditions on pn ---
     pn=0.; //non zero initial guess
     for( int axis=0; axis<numberOfDimensions; axis++ )
       for( int side=0; side<=1; side++ )
       {
         int is = 1-2*side; // is=1 on left, is=-1 on right
         if( boundaryCondition(side,axis)==dirichlet )
         {
           if( axis==0 )
           { // left or right side
             i1  = gridIndexRange(side,axis);
             int i1g = i1 - is; // index of ghost point

             #pragma omp parallel for default(shared) private(i2) num_threads(numThreads)
             for( i2=nd2a; i2<=nd2b; i2++ )
             {
               pn(i1,i2)  = pe(x(i1,i2,0),x(i1,i2,1),t+dt);
               pn(i1g,i2) = pe(x(i1g,i2,0),x(i1g,i2,1),t+dt);
             }
           }
           else
           { // bottom or top
             i2  = gridIndexRange(side,axis);
             int i2g = i2 - is; // index of ghost point

             #pragma omp parallel for default(shared) private(i1) num_threads(numThreads)
             for( i1=nd1a; i1<=nd1b; i1++ )
             {
               pn(i1,i2)  = pe(x(i1,i2,0),x(i1,i2,1),t+dt);
               pn(i1,i2g) = pe(x(i1,i2g,0),x(i1,i2g,1),t+dt);
             }
           }
         }
         else
         {
           printf("ERROR: unknown boundaryCondition=%d\n",boundaryCondition(side,axis));
           abort();
         }
         
       } // end for axis

     
     // conjugate gradient
     // Temporary arrays
     RealArray z(Rx,Ry), r(Rx,Ry), q(Rx,Ry);
 
     Real *z_p = z.getDataPointer();
     Real *r_p = r.getDataPointer();
     Real *q_p = q.getDataPointer();
     #define Z(i1,i2) z_p[(i1-nd1a)+nd1*(i2-nd2a)]
     #define R(i1,i2) r_p[(i1-nd1a)+nd1*(i2-nd2a)]
     #define Q(i1,i2) q_p[(i1-nd1a)+nd1*(i2-nd2a)]

     
     
     // initializing
     Real *pn_p = pn.getDataPointer();
     #define PN(i1,i2) pn_p[(i1-nd1a)+nd1*(i2-nd2a)]
     
     r=0.;

     #pragma omp parallel for default(shared) private(i1,i2) num_threads(numThreads)
     for( i2=n2a+1; i2<n2b; i2++ )
       for( i1=n1a+1; i1<n1b; i1++ )
       {
         R(i1,i2) = PRHS(i1,i2) + ( (pn(i1+1,i2) -2.*pn(i1,i2) + pn(i1-1,i2))/(dx[0]*dx[0])
                                    + (pn(i1,i2+1) -2.*pn(i1,i2) + pn(i1,i2-1))/(dx[0]*dx[0]) );
         Q(i1,i2) = R(i1,i2); // initial search direction
       }

     Real alpha,beta, rNormSquared, rNormSquaredNew, qz;
     rNormSquared=0.;

     #pragma omp parallel for default(shared) private(i1,i2) num_threads(numThreads) reduction(+:rNormSquared)
     for( i2=n2a+1; i2<n2b; i2++ )
       for( i1=n1a+1; i1<n1b; i1++ )
       {
         rNormSquared += R(i1,i2)*R(i1,i2);
       }

     Real res0 = sqrt(rNormSquared);
     Real maxRes=res0, maxResOld=res0;

     // iterate
     for( int m=0; m<maxIterations; m++ )
     {
       qz = 0.;

       #pragma omp parallel for default(shared) private(i1,i2) num_threads(numThreads) reduction(+:qz)
       for( i2=n2a+1; i2<n2b; i2++ )
         for( i1=n1a+1; i1<n1b; i1++ )
         {
           Z(i1,i2) = -( (Q(i1+1,i2) -2.*Q(i1,i2) + Q(i1-1,i2))/(dx[0]*dx[0]) + (Q(i1,i2+1) -2.*Q(i1,i2) + Q(i1,i2-1))/(dx[1]*dx[1]) );
           qz += Q(i1,i2)*Z(i1,i2);
         }
       
       alpha = rNormSquared/qz; // step length
       
       rNormSquaredNew=0.;

       #pragma omp parallel for default(shared) private(i1,i2) num_threads(numThreads) reduction(+:rNormSquaredNew)
       for( i2=n2a+1; i2<n2b; i2++ )
         for( i1=n1a+1; i1<n1b; i1++ )
         {
           PN(i1,i2) += alpha*Q(i1,i2); // new solution
           R(i1,i2) -= alpha*Z(i1,i2); // residual
           rNormSquaredNew += R(i1,i2)*R(i1,i2);
         }
       
       beta = rNormSquaredNew/rNormSquared; // improvement this step
       
       rNormSquared=rNormSquaredNew;
       
       // q = r + beta*p
       #pragma omp parallel for default(shared) private(i1,i2) num_threads(numThreads)
       for( i2=n2a+1; i2<n2b; i2++ )
         for( i1=n1a+1; i1<n1b; i1++ )
         {
           Q(i1,i2) = R(i1,i2) + beta*Q(i1,i2); // new search direction
         }
       
       if( ( m % intervalToCheckResidual) ==0 || m==(maxIterations-1) )
       {
         maxRes = getMaxResidual( pn,prhs,gridIndexRange,dx,numThreads );

         if( debug>0 )
           printf("CG: m=%3d: qz=%9.2e, alpha=%9.2e, beta=%9.2e, maxRes=%8.2e \n",m,qz,alpha,beta,maxRes);
         
         if( maxRes<tol )
         {
           //printf("  -> n=%d, CG pressure solve: numIterations=%d, maxRes=%8.2e \n",n,m,maxRes);
           break;
         }
         
       }
       
       
     }

     // --- update particle positions ---
     if( trackSnow > 0 )
     {
       int ix,iy;
       Real hvel, vvel;

       #pragma omp parallel for default(shared) private(i,ix,iy,hvel,vvel) num_threads(numThreads)
       for( i=0; i<nf; i++)
       {
         // find on grid
         ix = round((flakeLocation(i,0) - (xa - numGhost*dx[0]))/dx[0]) + n1a;
         iy = round((flakeLocation(i,1) - (ya - numGhost*dx[1]))/dx[1]) + n2a;

         // approximate velocities
         if( ix>=nd1a && ix<=nd1b && iy>=nd2a && iy<=nd2b )
         {
           hvel = un(ix,iy); // check this
           vvel = vn(ix,iy);
         }
         else
         {
           hvel = 0;
           vvel = 0;
         }

         flakeLocation(i,0) += dt*hvel;
         flakeLocation(i,1) += dt*(gravpull + vvel);

         //writing locations
         fprintf(snowfile,"sfLocation(%d,%s,%d)=[%g, %g];\n",i+1,":",n+2,flakeLocation(i,0),flakeLocation(i,1));
       }
       fprintf(snowfile,"\n");
     }

     
     cur = next;
   }
   /////// end time stepping loop ////////

   Real cpuTimeStep = getCPU()-cpu1;


  

   







  
   
   // =================== error computation ===================
   if( computeErrors > 0 )
   {
     RealArray & uc = ua[cur];
     RealArray & vc = va[cur];
     RealArray & pc = pa[cur];
     
     RealArray errU(Rx,Ry);
     RealArray errV(Rx,Ry);
     RealArray errP(Rx,Ry);
     
     Real maxErrU=0., maxNormU=0.;
     Real maxErrV=0., maxNormV=0.;
     Real maxErrP=0., maxNormP=0.;
     
     #pragma omp parallel for default(shared) private(i1,i2) num_threads(numThreads) reduction(max:maxErrU) reduction(max:maxErrV) reduction(max:maxErrP) reduction(max:maxNormU) reduction(max:maxNormV) reduction(max:maxNormP)
     for( i2=n2a; i2<=n2b; i2++ )
       for( i1=n1a; i1<=n1b; i1++ )
       {
         errU(i1,i2) = fabs(uc(i1,i2) - ue(x(i1,i2,0),x(i1,i2,1),tFinal));
         errV(i1,i2) = fabs(vc(i1,i2) - ve(x(i1,i2,0),x(i1,i2,1),tFinal));
         errP(i1,i2) = fabs(pc(i1,i2) - pe(x(i1,i2,0),x(i1,i2,1),tFinal));
         
         maxErrU = max(errU(i1,i2),maxErrU);
         maxErrV = max(errV(i1,i2),maxErrV);
         maxErrP = max(errP(i1,i2),maxErrP);
         
         maxNormU = max(abs(uc(i1,i2)),maxNormU);
         maxNormV = max(abs(vc(i1,i2)),maxNormV);
         maxNormP = max(abs(pc(i1,i2)),maxNormP);
       }
     
     
     
     maxErrU /= max(maxNormU,REAL_MIN); // relative error
     maxErrV /= max(maxNormV,REAL_MIN);
     maxErrP /= max(maxNormP,REAL_MIN);
     
     printf("  -> numSteps=%d: maxRelErrU=%8.2e, maxRelErrV=%8.2e, maxRelErrP=%8.2e\n",
            Nt, maxErrU,maxErrV,maxErrP);
   }
   




   printf("Computation finished: cpuTimeStep=%9.2e(s)\n",cpuTimeStep);
   if( trackSnow > 0 )
   {
     fclose(snowfile);
     printf("Wrote file [%s]\n\n",filename.c_str());
   }
   else
     printf("\n\n");






  
  
  return 0;
  
}
