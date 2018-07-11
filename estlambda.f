
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine est_lambda(YX, XX, m, f, k, muout, aout, zeromean)
      implicit none
      integer k, f, m
      double precision muout(k,m), aout(k,f,m), YX(k, f + 1, m)
      double precision XX((f+1)*(f+2)/2, k, m)
      double precision, dimension(f+1, f+1) :: xx_lower
      double precision, dimension(k, f+1, m) :: lambda_temp
      integer t, u, v, w, info, counter
      logical zeromean
      lambda_temp = YX !set equal to yx to overwrite in solve
      do t = 1,m !change to 1,m
        do w = 1,k! change to 1,f
          xx_lower = 0.0
          counter = 1
              do u=1,f+1 !u are rows
                do v=1,f+1 !v are cols
                  if (u <= v) then
                    xx_lower(u,v) = XX(counter, w, t)
                    counter = counter +1
                  end if
                enddo
              enddo 
          if (zeromean) then
            call dposv('U', f,1,xx_lower,f+1, lambda_temp(w,:,t),
     +      f+1, info)
          else
            call dposv('U',f + 1, 1, xx_lower, f+1,  
     +      lambda_temp(w,:, t), f + 1, info)
          end if
        enddo
      enddo
      aout(:,:,:)  = lambda_temp(:, 1:f, :)
	  if (.not. zeromean) then
	    muout = lambda_temp(:,f+1,:)
	  endif
      return
      end