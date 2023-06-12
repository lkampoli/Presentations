



  
      subroutine find_throat(x_bl_throat,y_bl_throat,cp_bl_throat)
      implicit real*8 (a-h,o-z)
      
      cp_bl_throat = 100
      x_bl_throat  = 0
      
      
      open(1,file='blade.dat',status='old',action='read')
      read(1,*)
      do while (.true.)
      read(1,*,end=100)x_bl,y_bl,x_norm,ps_bl,cp_bl
 
      if(cp_bl<0.and.x_bl>0.3d0.and.x_bl<0.75d0)then
          !print*,'on the suction side'
          if(cp_bl<cp_bl_throat)then
              cp_bl_throat = cp_bl
              x_bl_throat  = x_bl
              y_bl_throat  = y_bl
          endif
      endif


      enddo
      100 continue
      close(1)
      end

      subroutine find_te(x_bl_throat,y_bl_throat,cp_bl_throat)
      implicit real*8 (a-h,o-z)
      
      cp_bl_throat = 100
      x_bl_throat  = 0
      
      
      open(1,file='blade.dat',status='old',action='read')
      read(1,*)
      do while (.true.)
      read(1,*,end=100)x_bl,y_bl,x_norm,ps_bl,cp_bl
 
      if(cp_bl<0.and.x_bl>0.8d0.and.x_bl<0.85d0)then
          !print*,'on the suction side'
          if(cp_bl<cp_bl_throat)then
              cp_bl_throat = cp_bl
              x_bl_throat  = x_bl
              y_bl_throat  = y_bl
          endif
      endif


      enddo
      100 continue
      close(1)
      end




      subroutine get_thickness(xth,iofile,delta_1_th,delta_2_th)

      implicit real*8 (a-h,o-z)
      real*8 :: xth !target location 
      character*25 :: iofile
      real*8,dimension(:),allocatable:: x_bl,delta_1,delta_2

      !write(*,*)iofile
      !write(*,*)trim(iofile)
      open(1,file=trim(iofile))
      read(1,*)
      
      i=0
      do while (.true.)
      read(1,*,end=100)
      i=i+1
      enddo
      100 continue
      nlength=i
      !write(*,*)i
      allocate(x_bl(i),delta_1(i),delta_2(i))
      close(1)

      open(1,file=trim(iofile))
      read(1,*)
      do i=1,nlength
      read(1,*)dummy,x_bl(i),dummy,delta_1(i),delta_2(i)
      enddo
      close(1)

      do i=1,nlength-1
      if(xth>x_bl(i).and.xth<x_bl(i+1))then
        
!        |___________________________________|
!        |   x_bl(i)    |        delta_1(i)  |
!        |___________________________________|
!        |   xth        |           ?        |
!        |___________________________________|
!        |   x_bl(i+1)  |       delta_1(i+1) |
!        |___________________________________|

        delta_1_th = delta_1(i)+  (delta_1(i+1)-delta_1(i))/    &
                                  (x_bl(i+1)   -x_bl(i)   )*    &
                                  (xth         -x_bl(i)   )

        delta_2_th = delta_2(i)+  (delta_2(i+1)-delta_2(i))/    &
                                  (x_bl(i+1)   -x_bl(i)   )*    &
                                  (xth         -x_bl(i)   )
!       write(*,*)i,x_bl(i),x_bl(i+1)
      endif
      enddo
      close(1) 
      

      end subroutine
      
      
      
      
      
      program main

      implicit real*8 (a-h,o-z)
      real*8:: d1_ss,d2_ss,xth
      real*8:: d1_ps,d2_ps,xte
      character*25 :: iofile


      write(*,*)'On the suction side:'
      call find_throat(xth,yth,cpth)
      write(*,*)'(',xth,yth,')',cpth
      iofile='bound_ss_p1.dat'
      call get_thickness(xth,iofile,d1_ss,d2_ss)
      write(*,*)d1_ss,d2_ss
      
      write(*,*)'On the pressure side:'
      !call find_te(xte,yth,cpth)
      !write(*,*)'(',xte,yth,')',cpth 
      xte=0.842
      iofile='bound_ps_p1.dat'
      call get_thickness(xte,iofile,d1_ps,d2_ps)
      write(*,*)d1_ps,d2_ps

      write(*,*)'==========================================='
      write(*,*)'Final combined boundary layer thickness:'
      write(*,*)'delta:',d1_ss+d1_ps
      write(*,*)'theta:',d2_ss+d2_ps
      write(*,*)'==========================================='


      end
