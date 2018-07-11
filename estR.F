ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine est_R(aout, muout, S, N, YX, Rtype, k, f, m, 
     + Rout)
      implicit none
      mwSize k, f, m
      double precision S(k,m), N(k,m), Rout(k,m), aout(k,f,m)
      double precision YX(k, f + 1, m), Sout(k,m), muout(k,m)
      double precision, dimension(k, f+1, m) :: atemp
      double precision sigma2, sigma
      character(len=*) Rtype
      integer t
      atemp = 0.0d+0
      atemp(:,1:f,:) = aout(:,:,:)
      
      do t = 1,m
        atemp(:, f+1, t) = muout(:,t)
        S(:,t) = S(:,t) - sum((YX(:,:,t) * atemp(:,:,t)), DIM = 2)
      enddo

      sigma2 = sum(S)/sum(N)
      sigma = sqrt(sigma2)
      if (Rtype == 'iid') then
        Rout = sigma2
      else if (Rtype == 'diagonal') then
        Rout = S/N
        where (N == 0) Rout = 0 !avoid division by zero error when a line has no obs in a mix
      else if (Rtype == 'common') then
        do t=1,m
          Rout(:,t) = sum(S,2)/sum(N,2)  
        enddo
      end if
      return
      end