MODULE mod_gfsolver
  IMPLICIT NONE
  ! Contains:
  ! 1. subroutine sub_newton for basic Newton solver. It calls hybrid 
  !    Newton-bisection if max number of iterations is exceeded
  ! 2. subroutine sub_hybrid_nb for hybrid Newton-bisection solver

  ! Set kind of real arrays
  INTEGER, PARAMETER :: ki = SELECTED_REAL_KIND(p = 15, r = 30)

CONTAINS

  SUBROUTINE sub_newton(fcn,dfcn,xguess,x,minx_in,errrel_in,errabs_in,max_iter_in,fnorm)
    ! Find zeros of function fcn, given derivative dfcn and initial condition 
    ! xguess using the basic Newton method.
    ! Switches to hybrid Newton-bisection if max number of iterations is exceeded

    ! If the USER specificies the optional argument minx_in, only the region to
    ! the right of minx_in is searched.
    !
    ! Compulsory arguments:
    !	Inputs: 
    !		fcn: the name of the subroutine containing the function
    !           dfcn: the name of the subrouting containing the
    !		      derivative of the function
    !		xguess: initial condition for the function argument 
    !   Output:
    !		x: solution
    !
    ! Optional arguments:
    ! 	minx_in: 	 lower bound on solution x - Default -1.d15 
    ! 	errrerl_in:  convergence criterion on the argument - Default 1.d-8 
    ! 	errabs_in:	 convergence criterion on the function - Default 1.d-8
    !   max_iter_in: max # of iterations
    !   fnorm:       function value at the solution

    EXTERNAL fcn, dfcn

    REAL(ki),INTENT(in)                   :: xguess(1)
    REAL(ki), INTENT(out)                 :: x(1)
    REAL(ki),INTENT(in), OPTIONAL         :: errrel_in, errabs_in, minx_in
    INTEGER, intent(in), OPTIONAL         :: max_iter_in

    REAL(ki), INTENT(out), OPTIONAL       :: fnorm
    REAL(ki) :: dx, x_prime(1), errrel, errabs, minx, random_n, zero(1), df(1), scale
    REAL(ki), DIMENSION(2) :: f, xx 

    INTEGER  ::  iter, iter2, max_iter, max_iter2

    IF (PRESENT(minx_in)) THEN
      minx = minx_in
    ELSE
      minx = -1.d15
    END IF

    IF (PRESENT(errrel_in)) THEN
      errrel = errrel_in
    ELSE
      errrel = 1.d-8
    END IF

    IF (PRESENT(errabs_in)) THEN
      errabs = errabs_in
    ELSE
      errabs = 1.d-8 
    END IF

    IF (PRESENT(max_iter_in)) THEN
      max_iter = max_iter_in
    ELSE
      max_iter = 30
    END IF
    
    max_iter2 = 100
    x = xguess
    DO iter = 1, max_iter
      CALL fcn(1,x,zero)
      CALL dfcn(1,x,df)

      dx = zero(1)/df(1)
      x_prime = x - dx

      IF (ABS(df(1)) .LE. 1.d-20) THEN 
        PRINT *, "Zero derivative"
        STOP
      END IF

      IF (ABS(dx) .LE. errrel*(1.d0+ABS(x_prime(1)))) THEN
        IF (ABS(zero(1))<errabs) THEN
          x = x_prime
          IF (PRESENT(fnorm)) fnorm = zero(1)
          RETURN
        ENDIF
      END IF
      
      IF (x_prime(1).LE.minx) THEN
        ! Ensure the algorithm stays right of lower bound on x
        CALL RANDOM_NUMBER(random_n)
        x_prime = minx + random_n
      END IF

      x = x_prime

    END DO

    IF (iter>max_iter) THEN

      ! Locate bracketing interval before switching to 
      ! Newton-bisection hybrid     
      xx(1) = x(1)
      CALL fcn(1,xx(1),f(1))
      CALL dfcn(1,xx(1),df)
      dx = f(1)/df(1)
      scale = 3.d-2
      xx(2) = xx(1)
      DO iter2 = 1,max_iter2
        xx(2) = xx(2) - SIGN(scale,dx)
        CALL fcn(1,xx(2),f(2))
        IF (f(1)*f(2).LE.0.d0) EXIT
      ENDDO

      IF (iter2.GE.max_iter2 .AND. f(1)*f(2)>0) THEN
        PRINT *, "Exceeded max number of iterations in Newton"
        PRINT *, dx, xx, f
        STOP
      END IF
      
      IF (f(1)>f(2)) THEN
        xx = CSHIFT(xx,1)
      END IF
      !PRINT *, "Switching to Newton-bisection hybrid"
      CALL sub_hybrid_nb(fcn,dfcn,xx,x)
      IF (PRESENT(fnorm)) THEN
        CALL fcn(1,x(1),zero(1))
        fnorm = zero(1)
      END IF
        
    END IF
    
      
  END SUBROUTINE sub_newton


SUBROUTINE sub_hybrid_nb(fcn,dfcn,xguess_in,x,errrel_in,errabs_in,max_iter_in,fnorm)
    ! Find zeros of function fcn, given two points xguess_in(1), xguess_in(2) 
    ! using a hybrid Newton-bisection method.
    !
    ! Compulsory arguments:
    !	Inputs: 
    !		fcn: subroutine containing the function
    !   dfcn: subroutine containing the derivative of the function
    !		xguess: array of initial condition for the function argument 
    !   Output:
    !		x: solution
    !
    ! Optional arguments:
    ! 	errrerl_in:     convergence criterion on the argument - Default 1.d-8
    ! 	errabs_in:	convergence criterion on the function - Default 1.d-8
    !   max_iter_in:	max # of iterations
    !   fnorm:          function value at the solution

    EXTERNAL fcn,dfcn
    
    REAL(ki),INTENT(in)              :: xguess_in(2)
    REAL(ki), INTENT(out)            :: x(1)
    REAL(ki),INTENT(in), OPTIONAL    :: errrel_in, errabs_in
    INTEGER, INTENT(in), OPTIONAL    :: max_iter_in
    REAL(ki), INTENT(out), OPTIONAL  :: fnorm

    REAL(ki), DIMENSION(1) 	     :: df, f_x, zero
    REAL(ki), DIMENSION(2)           :: xguess, f
    REAL(ki)			     :: errabs, errrel, h, dx

    INTEGER  			     :: iter, max_iter

    
    IF (PRESENT(errrel_in)) THEN
      errrel = errrel_in
    ELSE
      errrel = 1.d-8
    END IF

    IF (PRESENT(errabs_in)) THEN
      errabs = errabs_in
    ELSE
      errabs = 1.d-8
    END IF

    IF (PRESENT(max_iter_in)) THEN
      max_iter = max_iter_in
    ELSE
      max_iter = 2000
    END IF

    xguess = xguess_in
    x = 0.5d0*SUM(xguess)
    CALL fcn(1,x(1),f_x(1))
        

    DO iter = 1,max_iter
      CALL fcn(1,x,zero)
      CALL dfcn(1,x,df)
      dx = zero(1)/df(1)
      x = x - dx
      
      IF ( (xguess(2)-x(1))*(x(1)-xguess(1))<0.d0 ) THEN
        ! If Newton step is outside interval bisect
        x = 0.5d0*SUM(xguess)
        h = xguess(2) - xguess(1)
      ELSE
        h = dx
      ENDIF

      CALL fcn(1,x(1),f_x(1))
      
      IF ( ABS(h) .LE. errrel*(1.d0+ABS(x(1))) ) THEN
        IF (ABS(f_x(1))<errabs) THEN
          IF (PRESENT(fnorm)) fnorm = f_x(1)
          RETURN
        ENDIF
      END IF
       
      IF (f_x(1)<0.d0) THEN 
        xguess(1) = x(1)
      ELSE
        xguess(2) = x(1)
      END IF
    END DO
    
    PRINT *, "Exceeded maximum number of iterations in hybrid method"
    STOP

  END SUBROUTINE sub_hybrid_nb

END MODULE mod_gfsolver
