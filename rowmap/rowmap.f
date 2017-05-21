c FILE: rowmap.f
      subroutine rowmap(n,f,ifcn,t,u,tend,hs,rtol,atol,itol,
     1                  jacv,ijacv,fdt,ifdt,solout,iout,work,
     2                  lwork,iwork,liwork,rpar,ipar,idid)
c-----!----------------------------------------------------------------- 
c
c   ROWMAP .. A ROW-code with Krylov techniques for large stiff ODEs.
c
c
c   ROWMAP solves the inital value problem for stiff systems of
c first order ODEs
c 
c      du     
c      -- = f(t,u),    u=(u(1),..,u(n)), f=(f(1),..f(n)).  
c      dt
c 
c
c   ROWMAP is based on the ROW-methods of order 4 of the code ROS4 of Hairer
c and Wanner (see [4]) and uses Krylov techniques for the solution of linear
c systems. By a special multiple Arnoldi process the order of the basic
c method is preserved with small Krylov dimensions.  Step size control is
c done by embedding with a method of order 3.
c
c Version of January 10, 1997.
c
c   Authors: 
c       H. Podhaisky, R. Weiner, 
c       Fachbereich Mathematik und Informatik, Universitaet Halle
c       06099 Halle, Germany
c       email: helmut.podhaisky@mathematik.uni-halle.de
c              ruediger.weiner@mathematik.uni-halle.de,
c     
c       B.A. Schmitt,
c       Fachbereich Mathematik, Universitaet Marburg
c       35032 Marburg, Germany  
c       email: schmitt@mathematik.uni-marburg.de
c
c
c-----!-----------------------------------------------------------------
c References:
c
c [1]  B.A. Schmitt and R. Weiner:
c         Matrix-free W-methods using a multiple Arnoldi Iteration,
c         APNUM 18(1995), 307-320
c
c [2]  R. Weiner, B.A. Schmitt an H. Podhaisky: 
c         ROWMAP - a ROW-code with Krylov techniques for large stiff
c         ODEs. Report 39, FB Mathematik und Informatik, 
c         Universitaet Halle, 1996
c       
c [3]  R.Weiner and B.A.Schmitt:
c         Consistency of Krylov-W-Methods in initial value problems,
c         Tech. Report 14, FB Mathematik und Informatik, 
c         Universitaet Halle, 1995
c
c [4]  E. Hairer and G. Wanner:
c         Solving Ordinary Differential Equations II, Springer-Verlag,
c         Second Edition, 1996
c   
c
c  For retrieving the latest version visit the ROWMAP World Wide Web 
c  homepage at URL :
c   http://www.mathematik.uni-halle.de/institute/numerik/software/
c
c-----!----------------------------------------------------------------- 
c
c   INPUT PARAMETERS
c   ----------------
c   
c   n      = Number of equations of the ODE system. 
c
c   f      = Name (external) of user-supplied subroutine computing the 
c            value f(t,u): 
c                subroutine f(n,t,u,udot,rpar,ipar)
c                real*8 t,u(n),udot(n),rpar(*)
c                integer n,ipar(*) 
c                udot(1)= ...
c            rpar,ipar (user parameters, see below).
c       
c   ifcn   = Switch:
c            ifcn = 0:  f(t,u) independent of t (autonomous)
c            ifcn = 1:  f(t,u) may depend on t (nonautonomous), default
c
c   t      = Initial t-value. Changed on exit.
c
c   u      = Initial values for u. Changed on exit.
c
c   tend   = Final t-value.
c
c   hs     = Initial step size guess. Changed on exit.
c            If hs=0d0 the code puts hs to a default value depending on 
c            atol, rtol and u. 
c           
c   rtol, atol = Relative and absolute error tolerances. They
c            can be both scalars or else both vectors of length n.
c            ROWMAP uses a weighted root-mean-square norm to measure the
c            size of error vectors. It is defined by
c
c                wrms=sqrt(sum[(err(i)/scal(i))**2, i=1,n]/n),
c
c            where
c                 scal(i)=atol+rtol*dabs(u(i))              (itol=0)
c            or
c                 scal(i)=atol(i)+rtol(i)*dabs(u(i))        (itol=1).
c
c   itol   = Switch for rtol,atol: 
c            itol=0: both are scalars (default).
c            itol=1: both are vectors.
c
c   jacv   = Name (external) of user-supplied subroutine computing 
c            the jacobian-vector product. 
c            z := Jac(t,u)*v, Jac=df/du:  
c                subroutine jacv(n,t,u,v,z,rpar,ipar)
c                real*8 t,u(n),v(n),z(n),rpar(*)
c                integer n,ipar(*)
c                z(1)=...   
c            rpar,ipar (see below).
c            This routine is only called if ijacv = 1. 
c            Supply a dummy subroutine in the case ijacv = 0
c            
c   ijacv  = Switch for the computation of the jacobian-vector products:
c            ijacv = 0: The products are approximated by finite differences.
c                       Subroutine jacv is never called (default).
c            ijacv = 1: The products are computed by the user-supplied 
c                       subroutine jacv.
c
c   fdt    = Name (external) of user-supplied subroutine computing the
c            partial derivate of f(t,u) with respect to t.
c            This routine is only called if idft=1 and ifcn=1.
c
c                subroutine fdt(n,t,u,ft,rpar,ipar) 
c                real*8 t,u(n),ft(n),rpar(*)
c                integer ipar(*)
c                ft(1)=...
c
c            rpar,ipar (see below).
c
c   ifdt   = Switch for the computation of df/dt:
c            ifdt = 0: df/dt is computed internally by finite differences,
c                      subroutine fdt is never called (default).
c            ifdt = 1: df/dt is supplied by subroutine fdt.
c
c   solout = Name (external) of user-supplied subroutine exporting the
c            numerical solution during integration. If  iout=1: 
c            solout it is called after every successful step. Supply
c            a dummy subroutine if iout = 0. It must have the form:
c
c               subroutine solout(n,told,tnew,uold,unew,fold,fnew,ucon,
c              1                  intr,rpar,ipar)
c               integer n,intr,ipar(*)
c               real*8 told,tnew,uold(n),unew(n),fold(n),fnew(n),ucon(n)
c               real*8 rpar(*)
c               ...
c               end
c       
c            "intr" serves to interrupt the integration. If solout sets
c             intr .lt. 0, then ROWMAP returns to the calling program.
c
c            solout may produce continuous output by calling the internal
c            subroutine "rowcon" (see below).
c
c   iout   = Gives information on the subroutine solout:
c               iout = 0: subroutine is never called
c               iout = 1: subroutine is used for output (default)
c
c   work   = Array of working space of length "lwork", changed on exit.
c            Serves as working space for all vectors and matrices. 
c            The array is used for optional input. For zero input, the
c            default values are set:
c
c            work(1), work(2) -  Parameters for step size selection. 
c               The new step size is chosen subject to the restriction
c                  work(1) <= hsnew/hsold <= work(2)
c                Default values: work(1)=0.25d0, work(2)=2d0
c   
c            work(3) - The safety factor in step size prediction.
c                Default value: work(3)=0.8d0
c
c            work(4) - UROUND, the machine precision/rounding unit. 
c                Default value: work(4)=1d-16
c
c            work(5) - KTOL, tolerance for the iterative solution of
c                the linear equations. The code keeps the weighted 
c                root-mean-square norm of the residual of the first stage
c                below KTOL/HS, where HS is the current step size
c                (see also ATOL/RTOL).
c                Default value: work(5)=1d-1.
c
c   lwork  = Length of array work in calling program.
c            "lwork" must be at least: 10+n*(mx+11)+mx*(mx+4), 
c            where mx=iwork(3).
c
c   iwork  = Integer working space of length "liwork". The array is used
c            for optional input. For zero input, default values are set:
c            iwork(1)   -  This is the maximal number of allowed steps.
c                Default value is iwork(1)=10000.
c
c            iwork(2)   -  Switch for the choice of integration method:
c                          (see [4], page 110)
c                           1  Method of Shampine
c                           2  Method GRK4T of Kaps-Rentrop 
c                           3  Method of Van Veldhuizen (gamma=1/2)
c                           4  Method of Van Veldhuizen ("D-stable")
c                           5  L-stable Method
c                           6  Method GRK4A of Kaps-Rentrop
c                Default value is iwork(2)=2.
c 
c            iwork(3)   -  The maximum Krylov dimension allowed (=mx). 
c                Default value is iwork(3)=70.
c
c            iwork(4)   -  The logical output unit number "lun" for
c                messages of ROWMAP. Default value is iwork(4)=6.
c                          
c            
c   liwork = Length of array iwork in calling program. 
c            "liwork" must be at least mx+20, where mx is iwork(3). 
c                       
c   rpar, ipar = Real and integer parameters (or parameter arrays)
c            that can be used for communication between your calling
c            program and the subroutines f, fdt, fjac, solout.
c
c   OUTPUT PARAMETERS
c   -----------------
c
c   t     = t-Value where the solution is computed (after successful  
c           return is t=tend)
c
c   u     = Solution at t.
c
c   hs    = Predicted step size from the last accepted step.
c
c   idid  = Reports on success upon return:
c         idid = 1    computation sucessful
c         idid = 2    computation interrupted by solout
c         idid = -1   stepsize too small, i.e. hs < 10*uround*dabs(t)
c         idid = -2   more than iwork(1) steps
c         idid = -3   input is not consistent
c         idid = -4   internal Krylov matrix is repeatedly singular
c
c   iwork(5) = Number of computed (accepted and rejected) steps.
c   iwork(6) = Number of rejected steps.
c   iwork(7) = Number of function evaluations.
c   iwork(8) = Number of jacobian-times-vector products.
c   iwork(9) = Minimum lenght of array work.
c   iwork(10)= Minimum lenght of array iwork.
c
c ROWMAP uses the following basic linear algebra modules (BLAS):
c       Level 1: DAXPY,DCOPY,DNORM2,DROTG,DROT,DDOT,DSCAL
c       Level 2: DGEMV,DTRSV
c
c-----!----------------------------------------------------------------- 
      implicit none
      integer n,lwork,liwork,iwork(liwork)
      integer itol,ijacv,ifdt,iout,ipar(*) ,ifcn,idid
      integer pmq,pmh,pmak1,pmak2,pmak3,pmak4,pmft,pmuu,pmfu0
      integer pmrhs,pml1,pml2,pml3,pml4,pmfm,pmcon,pmscal,pmmax,mx
      real*8 rpar(*) 
      real*8 t,u(n),tend,hs,rtol(*),atol(*),work(lwork),ktol
      external fdt,f,solout,jacv 
c
c Globals.
c
      real*8 uround
      integer cifcn,citol,cijacv,cifdt,ciout,method,lun
      common /rowmap0/uround,cifcn,citol,cijacv,cifdt,ciout,method,lun

c     For statistics.
      integer nfeval,nsteps,nstepsr,nfdt,njacv,nkrydim(4,3)
      common /rowmap1/nfeval,nsteps,nstepsr,nfdt,njacv,nkrydim
c 
c Parameters for step size control.
c
      integer nstpx
      real*8 fac1,fac2,fac3
      common  /rowmap2 /fac1,fac2,fac3,nstpx
c
c Pointers in array "work".
c
      mx=70
      if (iwork(3).gt.19) mx=iwork(3)
      pmq=10
      pmh=pmq+n*mx
      pmak1=pmh+mx*mx
      pmak2=pmak1+n
      pmak3=pmak2+n
      pmak4=pmak3+n
      pmft=pmak4+n
      pmuu=pmft+n
      pmfu0=pmuu+n
      pmrhs=pmfu0+n
      pml1=pmrhs+n
      pml2=pml1+mx
      pml3=pml2+mx
      pml4=pml3+mx
      pmfm=pml4+mx
      pmcon=pmfm+n
      pmscal=pmcon+n
      pmmax=pmscal+n
      iwork(9)=pmmax
c
c Check lwork und liwork.
c
      idid=0
      if (lwork .lt. pmmax) then
          write (lun,9950) lwork, pmmax
 9950     format ('Error in ROWMAP: "lwork"= ',i10,' too small. ', 
     1            'Minimum is',i10,'.')
          idid=-3
      end if
      iwork(10)=mx+20
      if (liwork.lt.iwork(10)) then
          write (lun,9951) liwork,iwork(10)
 9951     format ('Error in ROWMAP: "liwork"=',i6,' too small. ', 
     1           'Minimum is',i6,'.')
          idid=-3
      end if
      if (idid.eq.-3) return
c
c Check optional input.
c
      nstpx=10000
      if (1.lt.iwork(1)) nstpx=iwork(1)
      fac1=work(1)
      fac2=work(2)
      fac3=work(3)
      if (fac1.ge.1.or.fac1.lt.1d-2) fac1=0.25d0
      if (fac2.le.1.or.fac2.gt.1d+2) fac2=2d0
      if (fac3.lt.1d-1.or.fac3.ge.1d0) fac3=0.8d0
c The rounding unit.
      uround=1d-16
      if (work(4).gt.1d-40.and.work(4).lt.1d-5) uround=work(4)
      ktol=1d-1
      if (work(5).gt.0) ktol=work(5)
c     default ROW-method is GRK4T (method=2)
      method=2
      if (1.le.iwork(2).and.6.ge.iwork(2)) method=iwork(2) 
      cifcn=ifcn
      citol=itol
      cijacv=ijacv
      cifdt=ifdt
      ciout=iout
c     logical output unit number for messages
      lun=6
      if (iwork(4).gt.0) lun=iwork(4)
     
c 
c Call to core integrator.
c
      call rowmapc(n,f,t,u,tend,hs,rtol,atol,jacv,fdt,solout,
     1     ktol,work(pmq),work(pmh),
     2     work(pmak1),work(pmak2),
     3     work(pmak3),work(pmak4), work(pmft),
     4     work(pmuu),work(pmfu0),
     5     work(pmrhs),work(pml1),
     6     work(pml2),
     7     work(pml3),work(pml4),
     8     work(pmfm),work(pmcon),work(pmscal),
     9     iwork(20),mx,rpar,ipar,idid) 
c
c Statistics:
c
      iwork(5) = nsteps
      iwork(6) = nstepsr
      iwork(7) = nfeval
      iwork(8) = njacv
      return
      end

c-----!-----------------------------------------------------------------
c --- S U B R O U T I N E   R O W C O N 
c-----!-----------------------------------------------------------------
c This subroutine is for dense output. It can be called by the
c user-supplied subroutine "solout". A Hermite-interpolation 
c is used to compute an approximation (order 3) "ucon" to the solution
c of the ODE system at time s. The value s should lie in the
c interval [told,tnew].
c 
c Input 
c -----
c  n    = Dimension of the ODE system.
c  s    = Point of evaluation.
c  uold = Numerical solution for the ODE at time told, passed by solout.
c  fold = Value of the right hand side at time told, passed by solout.
c  unew = Numerical solution for the ODE at time tnew, passed by solout.
c  fnew = Value of the right hand side at time tnew, passed by solout.
c  
c Output
c ------
c  ucon = Approximation at time s.
c
      subroutine rowcon(n,s,told,tnew,uold,unew,fold,fnew,ucon)
      implicit none
      integer n
      real*8 s,told,tnew,uold(n),unew(n),fold(n),fnew(n),ucon(n)
      real*8 theta,theta2,thetam1,hs,g
      hs=tnew-told
      theta=(s-told)/hs
      theta2=theta**2d0
      thetam1=theta-1d0
      call dcopy(n,0d0,0,ucon,1)
      g=theta2*(3d0-2d0*theta)
      call daxpy(n,1d0-g,uold,1,ucon,1)
      call daxpy(n,g,unew,1,ucon,1)
      g=hs*theta*thetam1**2d0
      call daxpy(n,g,fold,1,ucon,1)
      g=hs*theta2*thetam1
      call daxpy(n,g,fnew,1,ucon,1)
      return
      end

c-----!----------------------------------------------------------------- 
c --- S U B R O U T I N E   R O W M A P C
c-----!----------------------------------------------------------------- 
c Here is the core integrator. ROWMAPC prepares the linear stage equations
c and calls subroutine "STAGE" for solving these. ROWMAPC contains the
c step size control. 
c   
      subroutine rowmapc(n,f,t,u,tend,hs,rtol,atol,jacv,fdt,solout,
     1              ktol,q,h,ak1,ak2,ak3,ak4,ft,uu,fu0,rhs,l1,l2,
     2               l3,l4,fm,ucon,scal,kl,mx,rpar,ipar,idid)
      implicit none
      integer mx,n,kl(mx),asv,mk,ipar(*),idid,i,j,intr,info,nsing
      real*8 t,u(n),tend,hs,rtol(*),atol(*),fm(n),ts,ktol,ktol1
      real*8 q(n,mx),h(mx,mx),ak1(n),ak2(n),ak3(n),ak4(n),
     1       ft(n),uu(n), fu0(n),l1(mx),l2(mx),l3(mx),l4(mx),
     2       rhs(n) , dnrm2,unorm,ucon(n),scal(n),
     3       rpar(*),delt,ehg,err,hnew,told
      external f, solout,fdt,jacv

      real*8 a21,a31,a32,c21,c31,c32,c41,c42,c43,b1,b2,b3,b4,e1,e2,e3,
     1   e4,gamma,c2,c3,d1,d2,d3,d4,g21,g31,g32,g41,g42,g43

      logical reached

c
c Common blocks.
c

C     Globals.   
      real*8 uround
      integer ifcn,itol,ijacv,ifdt,iout,method,lun
      common /rowmap0/uround,ifcn,itol,ijacv,ifdt,iout,method,lun
c     For statistics.
      integer nfeval,nsteps,nstepsr,nfdt,njacv,nkrydim(4,3)
      common /rowmap1/nfeval,nsteps,nstepsr,nfdt,njacv,nkrydim

c     For step size control.
      integer nstpx
      real*8 fac1,fac2,fac3
      common  /rowmap2 /fac1,fac2,fac3,nstpx

c     "unorm" is used in subroutine "ROWDQ".
      common /rownorm/unorm
c
c Initializations.
c 
      call rowcoe(method,a21,a31,a32,c21,c31,c32,c41,c42,c43,
     1            b1,b2,b3,b4,e1,e2,e3,e4,gamma,
     2            c2,c3,d1,d2,d3,d4,g21,g31,g32,g41,g42,g43)
      intr=0
      nfeval=0
      nsteps=0
      nstepsr=0
      njacv=0
      nfdt=0
      nsing=0
      reached=.false.
      do 101 i=1,4
         nkrydim(i,1)=mx
         nkrydim(i,2)=0
         nkrydim(i,3)=0
 101  continue
      call f(n,t,u,fm,rpar,ipar)
      nfeval=nfeval+1
      idid=1
c 
c Begin integration, loop for successful steps:
c 
 1    if (t.ge.tend.or.reached) return
      unorm=dnrm2(n,u,1)/dsqrt(dfloat(n))
      if (itol.ne.1) then
c       atol,rtol are scalars
        do 105 i=1,n
           scal(i)=1d0/(atol(1)+rtol(1)*dabs(u(i)))
 105    continue
      else
c       atol,rtol are vectors
        do 106 i=1,n
            scal(i)=1d0/(atol(i)+rtol(i)*dabs(u(i)))
 106    continue
      end if
      if (nsteps.eq.0.and.hs.le.0d0) then
c          Set initial step size.
           hs=dnrm2(n,scal,1)/dsqrt(dfloat(n))
         hs=(1d0/hs)**0.25*1d-1
      end if
      if ((t+hs*1.005).ge.tend) then
        hs=tend-t
        reached=.true.
      end if
c
c Compute the derivative f_t in the nonautonomous case.
c
      if (ifcn.ne.0) then
          nfdt=nfdt+1
          if (ifdt.ne.1) then
c         finite differences:
          delt=dsqrt(uround*max(1d-5,dabs(t)))
          call f(n,t+delt,u,ft,rpar,ipar)
          nfeval=nfeval+1
          call daxpy(n,-1d0,fm,1,ft,1)
          call dscal(n,1d0/delt,ft,1)
      else
c         user-supplied subroutine:
          call fdt(n,t,u,ft,rpar,ipar)
          end if
      end if
c
c Label for rejected steps:
c
 2    ehg=1.D0/(gamma*hs)
      ktol1=ktol/hs
      ts=t
      asv=1
      mk=0
c     Reset the target indices.
      do 102 j=1,mx
        kl(j)=0
  102 continue
      if (nsteps.ge.nstpx) then 
        write (lun,9990) t
        idid=-2
        write (lun,9992) nstpx 
        return
      end if
      nsteps=nsteps+1
c
c 1. stage.
c
      call dcopy (n,fm,1,rhs,1)
      call dscal (n,ehg,rhs,1)
      if (ifcn.ne.0) call daxpy(n,ehg*hs*d1,ft,1,rhs,1)
      call stage(q,h,rhs,l1,kl,asv,mk,mx,n,fm,fu0,u,t,1,ktol1,
     1           ehg,ak1,scal,f,jacv,fdt,rpar,ipar,info)
      if (info.lt.0) goto 3
c
c 2. stage. 
c
      call dcopy(n,u,1,uu,1)
      call daxpy(n,hs*a21,ak1,1,uu,1)
      ts=t+c2*hs 
      call f(n,ts,uu,rhs,rpar,ipar)
      nfeval=nfeval+1
      call daxpy(n,c21,ak1,1,rhs,1)
      call dscal (n,ehg,rhs,1)
      if (ifcn.ne.0) call daxpy(n,ehg*hs*d2,ft,1,rhs,1)
      call stage(q,h,rhs,l2,kl,asv,mk,mx,n,fm,fu0,u,t,2,ktol1,
     1           ehg,ak2,scal,f,jacv,fdt,rpar,ipar,info)
      if (info.lt.0) goto 3
      call daxpy(n,-c21,ak1,1,ak2,1)
c
c 3. stage.
c
      call dcopy(n,u,1,uu,1)
      call daxpy(n,hs*a31,ak1,1,uu,1)
      call daxpy(n,hs*a32,ak2,1,uu,1)
      ts=t+c3*hs 
      call f(n,ts,uu,rhs,rpar,ipar)
      nfeval=nfeval+1
      call dcopy(n,rhs,1,ak4,1)
      call daxpy(n,c31,ak1,1,rhs,1)
      call daxpy(n,c32,ak2,1,rhs,1)
      call dscal(n,ehg,rhs,1)
      if (ifcn.ne.0) call daxpy(n,ehg*hs*d3,ft,1,rhs,1)
      call stage(q,h,rhs,l3,kl,asv,mk,mx,n,fm,fu0,u,t,3,ktol1,
     1           ehg,ak3,scal,f,jacv,fdt,rpar,ipar,info)    
      if (info.lt.0) goto 3
      call daxpy(n,-c31,ak1,1,ak3,1)
      call daxpy(n,-c32,ak2,1,ak3,1)
c
c 4. stage.
c
      call dcopy(n,ak4,1,rhs,1)
      call daxpy(n,c41,ak1,1,rhs,1)
      call daxpy(n,c42,ak2,1,rhs,1)
      call daxpy(n,c43,ak3,1,rhs,1)
      call dscal(n,ehg,rhs,1)
      if (ifcn.ne.0) call daxpy(n,ehg*hs*d4,ft,1,rhs,1)
      call stage(q,h,rhs,l4,kl,asv,mk,mx,n,fm,fu0,u,t,4,ktol1,
     1           ehg,ak4,scal,f,jacv,fdt,rpar,ipar,info)
      if (info.lt.0) goto 3
      nsing=0
      call daxpy(n,-c41,ak1,1,ak4,1)
      call daxpy(n,-c42,ak2,1,ak4,1)
      call daxpy(n,-c43,ak3,1,ak4,1)
c
c New solution: uu = u + sum (h* b.i * ak.i,i=1..4).
c
      call dcopy (n,u,1,uu,1)
      call daxpy (n,hs*b1,ak1,1,uu,1)
      call daxpy (n,hs*b2,ak2,1,uu,1)
      call daxpy (n,hs*b3,ak3,1,uu,1)
      call daxpy (n,hs*b4,ak4,1,uu,1)
c
c Embedded solution: fu0 = sum (hs* e.i * ak.i, i=1..4).
c
      call dcopy (n,0d0,0,fu0,1)
      call daxpy (n,hs*e1,ak1,1,fu0,1)
      call daxpy (n,hs*e2,ak2,1,fu0,1)
      call daxpy (n,hs*e3,ak3,1,fu0,1)
      call daxpy (n,hs*e4,ak4,1,fu0,1)
c
c Error estimate, step size control.
c
       err=0d0
       do 550 i=1,n
          err=err+(fu0(i)*scal(i))**2 
  550  continue
      err=dsqrt(err/n)
      hnew=hs*dmin1(fac2,dmax1(fac1,(1d0/err)**0.25D0 * fac3))
      if (1d-1*hnew.le.dabs(t)*uround) then
          write (lun,9990) t
          idid=-1
          write (lun,9991) 
         return 
      end if
      if (err .lt. 1d0 ) then
c
c Step is accepted.
c
          told=t
          t=t+hs
          call dcopy(n,fm,1,fu0,1)
          call f(n,t,uu,fm,rpar,ipar)
          nfeval=nfeval+1
          if (iout.ne.0) 
     1        call solout(n,told,t,u,uu,fu0,fm,ucon,intr,rpar,ipar)
          call dcopy (n,uu,1,u,1)
          if (intr.lt.0) then
            idid=2
            write (*,9990) t
            write (*,9993)
            return
          end if
          hs=hnew
          goto 1
      else
c
c Step is rejected.
c
          nstepsr=nstepsr+1
          reached=.false.
          hs=hnew
          goto 2
      end if
c
c Matrix is singular.
c
 3    nsing=nsing+1
      write (lun,9995) 
      if (nsing.ge.3) then
        write (lun,9990) ts
        write (lun,9994)
        idid=-5
        return
      end if
      nstepsr=nstepsr+1
      hs=hs/2d0
      goto 2

 9990 format ('Exit of ROWMAP at t=',d10.3)
 9991 format ('Step size too small.')
 9992 format ('More than ',i7,' steps needed.')
 9993 format ('Computation interrupted by "solout".')
 9994 format ('Matrix is repeatedly singular')
 9995 format ('Warning: Matrix is singular.') 
      end           

c-----!-----------------------------------------------------------------
c --- S U B R O U T I N E   S T A G E 
c-----!-----------------------------------------------------------------
c The subroutine "STAGE" solves the linear equations by using the multiple
c Arnoldi process, see [1,2,3]. The first step incorporates the new right
c hand side into the Krylov subspace. Then the Krylov subspace is extended
c until the residual is small enough or the maximal dimensions are reached.
c Finally the solution "ak" is computed as a linear combination of
c orthogonalized Krylov vectors.
c
c
      subroutine stage(q,h,rhs,l,kl,asv,mk,mx,n,fm,fu0,u,t,st,ktol1,
     1                 ehg,ak,scal,f,jacv,fdt,rpar,ipar,info)
      implicit none
      integer mx,n,mk,i,asv,m,ikrdim,info
      integer kl(mx),st,krydim(4,4),ipar(*)
      logical done
      real*8 q(n,mx),h(mx,mx),rhs(n),l(mx),dnrm2,ktol1,nrmv,nrmrhs
      real*8 fm(n),fu0(n),u(n),t,ehg,dfkt
      real*8 ak(n),scal(n),rpar(*)
      external f,jacv,fdt
c
c     Maximal Krylov dimension increments per stage, see [3]
c     (autonomous/nonautonomous). For the methods shamp, grk4t,
c     veldd, velds and lstab (b3.eq.0)
c
      data (krydim(i,1),i=1,4) /-12,3,1,3/
      data (krydim(i,2),i=1,4) /-16,5,1,5/
c
c    ... and for grk4a (b3.ne.0):
c
      data (krydim(i,3),i=1,4) /-15,3,4,3/
      data (krydim(i,4),i=1,4) /-19,5,4,5/

c     Globals.
      real*8 uround
      integer ifcn,itol,ijacv,ifdt,iout,method,lun
      common /rowmap0/uround,ifcn,itol,ijacv,ifdt,iout,method,lun

c     Statistics.
      integer nfeval,nsteps,nstepsr,nfdt,njacv,nkrydim(4,3)
      common /rowmap1/nfeval,nsteps,nstepsr,nfdt,njacv,nkrydim

      ikrdim=2*(method/6)+ifcn+1
      info=0
      nrmrhs=dnrm2(n,rhs,1)
      if (mk.eq.0) then
          call dcopy(mx,0d0,0,l,1)
      else 
          call orthov(n,mk,q,rhs,l,mx)
      end if
      nrmv=dnrm2(n,rhs,1)
      if (nrmv/(1d-12+nrmrhs) .ge. 1d-8.or.st.lt.4) then
c
c Insert the new (orthogonalized) right hand side "rhs" as 
c new column in Q.
c
         if (mk.gt.0) asv=asv+1
c        Shift last vectors by 1.
         do 100 i=asv-1,1,-1
            call dcopy(n,q(1,mk+i),1,q(1,mk+i+1),1)
 100     continue
         do 101 i=1,mk
            if(kl(i).ge.mk+1) kl(i)=kl(i)+1
 101     continue
         l(mk+1)=nrmv
         call dcopy(n,rhs,1,q(1,mk+1),1)
         call dscal(n,1d0/l(mk+1),q(1,mk+1),1)
      end if
      m=0
      done=.false.
c
c Begin Multiple Arnoldi process.
c
 1000 m=m+1
      if (m.gt.mk) then
          kl(m)=m+asv
      end if
      call kryarn(n,t,u,mk,m,kl,q,h,mx,ehg,l,fm,fu0,f,jacv,rpar,ipar)
      if (m.gt.mk-st+asv.and.m.ge.asv.and.dabs(h(m,m)).lt.uround) then
c        matrix is singular
         info=-1
         return
      end if
      if (m.gt.mk-st+asv.and.m.ge.asv) then
        call krdfkt(n,m,kl,q,h,mx,l,fu0,dfkt,st,asv,ehg,scal)
      else 
        dfkt=1d100
      end if
      if (dfkt.lt.ktol1.and.st.gt.1) done=.true.
      if (dfkt.lt.ktol1.and.m.ge.4) done=.true.
      if(st.gt.1.and.(m-mk).ge.krydim(st,ikrdim)) done=.true.
      if(st.eq.1.and.m.ge.krydim(1,ikrdim)+mx) done=.true.
      if (.not. done) goto 1000
c
c End of multiple Arnoldi process for current stage.
c
      mk=m
c     Used Krylov dimensions.
      nkrydim(st,1)=min0(nkrydim(st,1),mk)
      nkrydim(st,2)=max0(nkrydim(st,2),mk)
      nkrydim(st,3)=nkrydim(st,3)+mk
c     Backsubstitution.
      call dtrsv('u','n','n',m,h,mx,l,1)
c     Compute ak = Q*l.
      call dgemv('n',n,mk,1d0,q,n,l,1,0d0,ak,1)
      return
      end

c-----!-----------------------------------------------------------------
c --- S U B R O U T I N E   O R T H O V 
c-----!-----------------------------------------------------------------
c This subroutine orthogonalizes v with respect to vectors q(.,i),i=1,m,
c by a modified Gram-Schmitt-Algorithm. Dotproducts are stored in "l",
c i.e. l=Q'v.
c    
c
      subroutine orthov(n,m,q,v,l,mx)
      implicit none
      integer n,m,mx,i
      real*8 q(n,*), v(n),l(mx),s,ddot 
      do 100 i=1,m
          s=ddot(n,q(1,i),1,v,1)
          l(i)=s
          call daxpy(n,-s,q(1,i),1,v,1)
 100  continue
      call dcopy(mx-m,0d0,0,l(m+1),1)
      return
      end

c-----!-----------------------------------------------------------------
c --- S U B R O U T I N E   R E O R T N 
c-----!-----------------------------------------------------------------
c The vector q(.,m) is incorporated in the Krylov space, it is scaled to
c norm one with h(m,j), where kl(j)=m. The candidates for further subspace
c extentions have to be re-orthogonalized to q(.,m).
c
      subroutine reortn(n,m,kl,q,h,mx)
      implicit none
      integer n,m,kl(*),mx,j0,j
      real*8 q(n,*),h(mx,*),s,ddot
      j0=m
      do 100 j=m-1,1,-1
 100  if (kl(j).ge.m) j0=j
      do 150 j=j0,m-1
          if (kl(j).eq.m.and.h(m,j).gt.1d-14)
     1        call dscal(n,1d0/h(m,j),q(1,m),1)
 150  continue
      do 200 j=j0,m-1
      if (kl(j).gt.m) then
          s=ddot(n,q(1,m),1,q(1,kl(j)),1)
          h(m,j)=s
          call daxpy(n,-s,q(1,m),1,q(1,kl(j)),1)
      end if
 200  continue
      return
      end
c-----!-----------------------------------------------------------------
c --- S U B R O U T I N E   K R Y A R N 
c-----!-----------------------------------------------------------------
c This subroutine extends the Krylov subspace by one vector and updates
c the QR-decomposition of matrix "h".
c

      subroutine kryarn(n,t,u,mk,m,kl,q,h,mx,ehg,l,fm,fu0,f,
     1          jacv,rpar,ipar)
      implicit none
      integer n,mk,m,kl(*),mx, i,j,j0

      real*8 uround
      integer ifcn,itol,ijacv,ifdt,iout,method,lun
      common /rowmap0/uround,ifcn,itol,ijacv,ifdt,iout,method,lun

      integer nfeval,nsteps,nstepsr,nfdt,njacv,nkrydim(4,3),ipar(*)
      common /rowmap1/nfeval,nsteps,nstepsr,nfdt,njacv,nkrydim
      real*8  u(n),t,q(n,*),h(mx,*),ehg,l(*), c,s,fm(n),fu0(n)
      real*8  rpar(*)
      external f,jacv
      if (m.gt.mk) call reortn(n,m,kl,q,h,mx)
      l(m) = -l(m)
      j0 = m
      do 50 j=m-1,1,-1
  50    if (kl(j).ge.m) j0=j
c     generate new Givens rotations 
      do 100 j=j0,m-1
        if (m .gt. mk) then
         call drotg(h(j,j),h(m,j),c,s)
         call drot(m-j-1,h(j,j+1),mx,h(m,j+1),mx,c,s)
        else
c        reuse old rotations
         call ztocs(h(m,j),c,s)
        endif
c       apply to right hand side "l"
        call drot(1,l(j),1,l(m),1,c,s)
 100  continue

c
c Multiple Arnoldi-Step with QR-Decomposition.
c
      if (m.le.mk) return
c      Krylov-Step: q(.,kl(m)) = (I - QQ') A q(.,m) 
      if (ijacv.ne.1) then
        call rowdq(n,t,u,q(1,m),fu0,fm,q(1,kl(m)),f,rpar,ipar)
        nfeval=nfeval+1
      else
        call jacv(n,t,u,q(1,m),q(1,kl(m)),rpar,ipar)
       end if
      njacv=njacv+1
      call orthov(n,m,q,q(1,kl(m)),h(1,m),mx)
      h(m,m) = h(m,m)-ehg
c     update the QR-decomposition of matrix "h" 
      do 210 i=2,m
       do 200 j=i-1,1,-1
  200   if (kl(j).ge.i) j0=j
       do 210 j=j0,i-1
         call ztocs(h(i,j),c,s)
         call drot(1,h(j,m),1,h(i,m),1,c,s)
  210 continue
      return
      end

c-----!-----------------------------------------------------------------
c --- S U B R O U T I N E   K R D F K T 
c-----!-----------------------------------------------------------------
c This subroutine computes the weighted root-mean-square norm of the residual
c
c               r = (I-gh A) Q l - rhs.
c
      subroutine krdfkt(n,m,kl,q,h,mx,r,fu0,dn,st,asv,egh,scal)
      implicit none
      integer n,m,mx,kl(mx),i,j,li,ls,mi,st,asv,j0
      real*8 q(n,*),h(mx,*),r(m),fu0(*),dn,s,qn,dnrm2,egh,gh
      real*8 scal(n),wrmsq

      dn = 0D0
      gh=1d0/egh
      j0 = m
      do 10 j=m,1,-1
  10    if (kl(j).gt.m) j0=j
      ls = m-j0+1

      do 110 i=0,ls-1
       mi = m-i
       s = r(mi)
       do 100 j=0,i-1
  100   s = s-h(mi,m-j)*fu0(ls-j)
       s = s/h(mi,mi)
       fu0(ls-i) = s
       s = dabs(s)
       li = kl(mi)
       if (li.gt.m) then
         qn = dnrm2(n,q(1,li),1)
         h(m+1,mi) = qn
         wrmsq=0d0
         do 101 j=1,n
            wrmsq=wrmsq+(q(j,li)*scal(j))**2 
 101     continue
         wrmsq=dsqrt(wrmsq/n)
         dn = dn + s*wrmsq 
       endif
  110 continue
       dn=dn*gh
      return
      end

c-----!-----------------------------------------------------------------
c --- S U B R O U T I N E   Z T O C S 
c-----!-----------------------------------------------------------------
c Reconstructs old Givens rotations computed by DROTG.
c
      subroutine ztocs(z,c,s)
      implicit none
      real*8 z,c,s
      if (dabs(z).lt.1d0) then
          s=z
          c=dsqrt(1d0-s*s)
      else if (dabs(z).gt.1d0) then
          c=1d0/z
          s=dsqrt(1d0-c*c)
      else 
          c=0d0
          s=1d0
      end if
      return
      end

c-----!----------------------------------------------------------------- 
c --- S U B R O U T I N E   R O W D Q
c-----!----------------------------------------------------------------- 
c This subroutine approximates the Jacobian-vector product
c by finite differences
c
c     v = (f(u+delta y)-f(u))/delta = A y + O(||u||delta^2)     
c

      subroutine rowdq(n,t,u,y,fu0,fm,v,f,rpar,ipar)
      implicit none
      integer n,ipar(*)
      real*8  t,u(n),y(n),fu0(n),fm(n),v(n),delta
      real*8  rpar(*),eddelta,unorm
      external f
      common /rownorm/ unorm
      delta=1.d-7*dmax1(1d-5,unorm)
      eddelta=1d0/delta
      call dcopy(n,u,1,fu0,1)
      call daxpy(n,delta,y,1,fu0,1)
      call f(n,t,fu0,v,rpar,ipar)
      call daxpy(n,-1d0,fm,1,v,1)
      call dscal(n,eddelta,v,1)
      return
      end

c-----!----------------------------------------------------------------- 
c --- S U B R O U T I N E   R O W C O E
c-----!----------------------------------------------------------------- 
c This subroutine loads the coefficients of the chosen method.
c
      subroutine rowcoe(method,a21,a31,a32,c21,c31,c32,c41,c42,c43,
     1            b1,b2,b3,b4,e1,e2,e3,e4,gamma,
     2            c2,c3,d1,d2,d3,d4,g21,g31,g32,g41,g42,g43)
      implicit none
      integer method
      real*8 a21,a31,a32,c21,c31,c32,c41,c42,c43,b1,b2,b3,b4,e1,e2,e3,
     1   e4,gamma,c2,c3,d1,d2,d3,d4,g21,g31,g32,g41,g42,g43
      goto (1,2,3,4,5,6),method
 10   continue
      c21=g21/gamma
      c31=g31/gamma
      c32=g32/gamma
      c41=g41/gamma
      c42=g42/gamma
      c43=g43/gamma
      c2=a21
      c3=a31+a32
      d1=gamma
      d2=gamma+g21
      d3=gamma+g31+g32
      d4=gamma+g41+g42+g43
      e1=e1-b1
      e2=e2-b2
      e3=e3-b3
      e4=e4-b4
      return

c -- SHAMP
 1    a21= 1.00000000000000  D0
      a31=0.48000000000000   D0 
      a32=0.12000000000000   D0
      g21=-2.00000000000000  D0
      g31=1.32000000000000   D0 
      g32=0.60000000000000   D0
      g41=-0.05600000000000  D0
      g42= -0.22800000000000 D0 
      g43= -0.10000000000000 D0
      b1=0.29629629629630    D0
      b2=0.12500000000000    D0
      b3=0                   D0
      b4=0.57870370370370    D0
      e1=0                   D0
      e2= -0.04166666666667  D0 
      e3= -0.11574074074074  D0
      e4= 1.15740740740741   D0
      GAMMA= 0.50000000000000D0
      goto 10
c --  GRK4T
 2    g21=-0.27062966775244  D0
      g31=0.31125448329409   D0
      g32=0.00852445628482   D0
      g41=0.28281683204353   D0 
      g42=-0.45795948328073  D0
      g43=-0.11120833333333  D0
      b1=0.21748737165273    D0
      b2=0.48622903799012    D0
      b3=0.00000000000000    D0
      b4=0.29628359035715    D0
      a21=0.46200000000000   D0
      a31=-0.08156681683272  D0
      a32=0.96177515016606   D0
      e1=-0.71708850449933   D0
      e2=1.77617912176104    D0
      e3=-0.05909061726171   D0
      e4=0.00000000000000    D0
      GAMMA=0.23100000000000D0
      goto 10
c -- VELDS
 3    g21=-2.00000000000000  D0
      g31=-1.00000000000000  D0
      g32=-0.25000000000000  D0
      g41=-0.37500000000000  D0
      g42=-0.37500000000000  D0
      g43= 0.50000000000000  D0
      b1=0.16666666666667    D0
      b2=0.16666666666667    D0
      b3=0                   D0
      b4=0.66666666666667    D0
      a21=1.00000000000000   D0
      a31=0.37500000000000   D0
      a32=0.12500000000000   D0
      e1=1.16666666666667    D0
      e2=0.50000000000000    D0
      e3=-0.66666666666667   D0
      e4=0                   D0
      GAMMA=0.50000000000000 D0
      goto 10
c -- VELDD
 4    g21=-0.27170214984937   D0
      g31=0.20011014796684    D0
      g32=0.09194078770500    D0
      g41=0.35990464608231    D0
      g42=-0.52236799086101   D0
      g43=-0.10130100942441   D0
      b1=0.20961757675658     D0
      b2=0.48433148684810     D0
      b3=0.0                  D0
      b4=0.30605093639532     D0
      a21=0.45141622964514    D0
      a31=-0.15773202438639   D0
      a32=1.03332491898823    D0
      e1=-0.74638173030838    D0
      e2=1.78642253324799     D0
      e3=-0.04004080293962    D0
      e4=0.0                  D0
      GAMMA=0.22570811482257  D0
      goto 10
c -- LSTAB
 5    g21=-2.34201389131923   D0
      g31=-0.02735980356646   D0
      g32=0.21380314735851    D0
      g41=-0.25909062216449   D0
      g42=-0.19059462272997   D0
      g43=-0.22803686381559   D0
      b1=0.32453574762832     D0
      b2=0.04908429214667     D0
      b3=0.00000000000000     D0
      b4=0.62637996022502     D0
      a21=1.14564000000000    D0
      a31=0.52092209544722    D0
      a32=0.13429476836837    D0
      e1=0.61994881642181     D0
      e2=0.19268272217757     D0
      e3=0.18736846140061     D0
      e4=0                    D0
      GAMMA=0.57282000000000  D0
      goto 10
c -- GRK4A
 6    a21=0.43800000000000   D0
      a31=0.79692045793846   D0
      a32=0.07307954206154   D0
      g21=-0.76767239548409  D0
      g31=-0.85167532374233  D0
      g32=0.52296728918805   D0
      g41=0.28846310954547   D0
      g42=0.08802142733812   D0
      g43=-0.33738984062673  D0
      b1=0.19929327570063    D0
      b2=0.48264523567374    D0
      b3=0.06806148862563    D0
      b4=0.25000000000000    D0
      e1=0.34632583375795    D0
      e2=0.28569317571228    D0
      e3=0.36798099052978    D0
      e4= 0.0                D0
      GAMMA=0.39500000000000 D0
      goto 10
      end
