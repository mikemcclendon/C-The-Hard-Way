#include "fintrf.h"

 ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function logsumexp(v, m)
      integer m
      double precision, dimension(m) :: v !need to declare dim to work
      double precision :: logsumexp, c
      c = maxval(v)
      logsumexp = log(sum(exp(v-c))) + c
      end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine EMFA(a, R, mu, sr, irs, jcs, P, Rtype,  n,
     + aout, muout, Rout, Pout, sumlogp, l_sr, l_irs, l_jcs, k,f,m)

      implicit none
 
      integer k, f, m, l_sr, l_irs, l_jcs 
      double precision, intent(in) :: a(k,f,m), R(k,m), mu(k,m)
      double precision, intent(in) :: sr(l_sr) , P(m)
      integer, intent(in) :: irs(l_irs)        
      integer, intent(in) :: jcs(l_jcs)  
      double precision, intent(out) :: aout(k,f, m), Rout(k,m), muout(k,m)
      double precision, intent(out) :: Pout(m), sumlogp

      integer, dimension(k) :: I, J 
      logical, dimension(k) :: msk
      double precision, dimension (k, f + 1, m) :: YX
      double precision, dimension ((f+1)*(f+2)/2, m) :: XXtmp
      double precision, dimension ((f+1)*(f+2)/2, k, m) :: XX
      double precision, dimension (k, f, m) :: aiR, YXtmp, aout
      double precision, dimension (k, m) :: iR, S, N_mix, muout
      double precision, dimension (k, m) :: Rout
      double precision, dimension (k, k) :: Ryt
      double precision, dimension (k) :: y, ym, ymiR
      double precision, dimension (k, f) :: Qa_t, RxyRyti
      double precision, dimension (f, f) :: iQ, Q
      double precision, dimension (f) :: q_lh
      double precision, dimension (m) :: lh, xi, sumP
      double precision, dimension(f + 1) :: x
      double precision, dimension(f + 1, f + 1) :: xx_lower
      double precision logdet, one, zero, logp, sumlogp
      double precision ddot, final_lh, logsumexp
      
      integer, dimension(k)::  packarr
      integer INFO, oneI, kt, t, u, findloc, l, v, counter
     
      character uplo, trans, side
      character(len=*) Rtype
      logical zeromean
	  
Cf2py intent(in) k
Cf2py intent(in) n
Cf2py intent(in) f
Cf2py intent(in) m
Cf2py intent(in) l_sr
Cf2py intent(in) l_irs
Cf2py intent(in) l_jcs
Cf2py intent(in) a
Cf2py depend(k,f,m) a
Cf2py intent(in) R
Cf2py depend(k,m) R
Cf2py intent(in) mu
Cf2py depend(k,m) mu
Cf2py intent(in) sr
Cf2py depend(l_sr) sr
Cf2py intent(in) irs
Cf2py depend(l_irs) irs
Cf2py intent(in) jcs
Cf2py depend(l_jcs) jcs
Cf2py intent(in) P
Cf2py depend(m) P
Cf2py intent(out) Rout
Cf2py depend(k,m) Rout
Cf2py intent(out) Pout
Cf2py depend(m) Pout
Cf2py intent(out) muout
Cf2py depend(k,m) muout
Cf2oy intent(out) sumlogp

      uplo = 'U'
      trans = 'N'
      side = 'L'
      oneI=1

      one = 1.0d+0
      zero = 0.0d+0
      iQ = zero
      q_lh = zero
      YXtmp = zero
      XXtmp = zero
      YX = zero
      XX = zero
      S = zero
      N_mix = zero
      xi = zero
	  muout= zero

c     check if mean = 0 and create flag
      zeromean = .FALSE.
      if (sum(abs(mu)) == zero) zeromean = .TRUE. 

c     get iR and aiR
      iR = one/R
      do t=1, m
        do l=1,f
          aiR(:,l,t) = iR(:,t) * a(:,l,t)
        enddo
      enddo

c     begin loop through records
      do t=1,n
        kt = jcs(t+1)-jcs(t)
       
        I(1:kt)= irs(jcs(t)+1: jcs(t+1))+1
        y(1:kt) = sr( jcs(t)+1: jcs(t+1))  
        msk = 0
        msk(I(1:kt))= 1
        J(1:(k-kt)) = pack(packarr,msk==0)

        do l=1, m
          ym(1:kt) = y(1:kt) - mu(I(1:kt), l)
          logdet = zero
          x = zero
          Q = zero
          counter = 1
          xx_lower = zero
          YXtmp(:,:,l) = zero
          iQ = zero
          q_lh = zero
          Qa_t = zero
          IF (kt < f) THEN
            Ryt = zero
            RxyRyti = zero
            do u = 1,kt !top ktxkt of Ryt is full
              Ryt(u,u) = R(I(u), l)
            enddo
            call dsyrk(UPLO, 'N', kt, f, one, a(I(1:kt),:,l), kt, one, 
     +      Ryt, k)
            call dpotrf(UPLO, kt, Ryt, k, INFO) !Ryt now C
            call dtrsv(UPLO, 'T', 'N', kt, Ryt, k, ym, 1) !ym now icym
            do u = 1,kt
              logdet = logdet + log(Ryt(u,u))
            enddo
            lh(l) = log(P(l)) - kt*0.918938533204673d+0 - logdet
     +      - 0.5*sum(ym(1:kt)**2)
            call dtrsv(UPLO, 'N', 'N', kt, Ryt, k, ym, 1) !ym is iCym\C
            call dgemv('T', kt, f, one, a(I(1:kt),:,l) , kt, ym, 1, 
     +      zero, x, 1)  
            x(f+1) = one   
            RxyRyti(1:kt, :) = a(I(1:kt),:,l)
            call dtrsm('L', UPLO, 'T', 'N', kt, f, one, Ryt, k, 
     +      RxyRyti, k)
            do u = 1,f
              Q(u,u) = one
            enddo
            call dsyrk(UPLO, 'T', f, kt, -one, RxyRyti, k, one, Q, f)
            call dsyr(UPLO, f+1,one, x, 1, xx_lower, f+1)
            xx_lower(1:f,1:f) = xx_lower(1:f,1:f) + Q !upper triangular
          ELSE       !kt > f
            call dgemm('T', 'N', f, f, kt, one, a(I(1:kt),:,l), kt, 
     +      aiR(I(1:kt),:,l), kt ,zero, iQ, f)
            do u=1,f
              iQ(u,u) = iQ(u,u) + 1.0
            enddo
            call  dpotrf(UPLO, f, iQ, f, INFO )  
            ymiR(1:kt) = ym(1:kt) * iR(I(1:kt),l)
            call dgemv('T', kt, f, one, a(I(1:kt),:,l), kt, ymiR, 1,  
     +      zero, q_lh, 1)
            call dtrsv(UPLO, 'T', 'N', f, iQ, f, q_lh, 1)
            do u = 1,f
              logdet = logdet + log(iQ(u,u))
            enddo
            lh(l) = log(P(l)) -kt*0.918938533204673d+0 + .5*(-2*logdet 
     +      - sum(log(R(I(1:kt),l))) - ddot(kt, ym, 1, ymiR, 1) +
     +      ddot(f, q_lh, 1, q_lh, 1))
            call dpotri(uplo, f, iQ, f, info) !iQ now cholesky
            call dsymm('R', uplo, kt, f, one, iQ, f, a(I(1:kt),:,l),
     +      kt, zero, Qa_t, k)
            call dgemv('T',kt, f, one, Qa_t, k, ymiR, 1, zero, x, 1) 
            x(f+1) = one
            call dsyr(UPLO,f+1,one, x, 1, xx_lower, f+1)
            xx_lower(1:f,1:f) = xx_lower(1:f,1:f) + iQ !upper triangular
          END IF

          call dger(kt, f, one, y, 1, x, 1,
     +    YXtmp(:,:,l), k)

          do u=1,f+1 !u are rows
            do v=1,f+1 !v are cols
              if (u <= v) then
                XXtmp(counter,l) = xx_lower(u,v)
                counter = counter +1
              end if
            enddo
          enddo
          
        enddo !finish first loop 

        logp = logsumexp(lh, m) 
        sumlogp = sumlogp + logp
        xi = exp(lh - logp)
        sumP = sumP + xi

        do l=1,m 
          YX(I(1:kt),1:f,l) = YX(I(1:kt),1:f,l) + xi(l)*YXtmp(1:kt,:,l)
          YX(I(1:kt),f+1,l) = YX(I(1:kt),f+1,l) + xi(l)*y(1:kt)
          S(I(1:kt),l) = S(I(1:kt),l) + xi(l)*(y(1:kt)**2)
          N_mix(I(1:kt),l)= N_mix(I(1:kt),l) +  xi(l)
 
          do u = 1,kt   !update xx
            XX(:,I(u),l) = XX(:,I(u),l) + xi(l)*XXtmp(:,l)
          enddo
        enddo
      enddo
      call est_lambda(YX, XX, m, f, k, muout, aout, zeromean)
      call est_R(aout, muout, S, N_mix, YX, Rtype,
     +  k, f, m, Rout)

      Pout = sumP/n

      return
      end
      
