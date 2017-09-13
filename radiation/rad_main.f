

      program rad_main

      implicit none

      include 'rad_coms.in'

      double precision alphas,dNg_over_dxdydt,splittingP,tau_f,try
      external dNg_over_dxdydt,try
      double precision delta_Ng_over_dt,error_Ng
      double precision kT,x,y,argu(2)
      integer ndim,itmx,ncall,nprn,init
      double precision region(4),sd,chi2a,tgral
      integer i,j,k
      character*77 file23

      time_init=0.d0
      HQmass=cMass

      ndim=2
      itmx=5
      ncall=100000
      nprn=0
      init=0

      region(1)=0.d0
      region(2)=0.d0

      region(3)=1.d0
      region(4)=1.d0
      
      delta_HQener=HQener_max/HQener_gn
      delta_t=(t_max-t_min)/t_gn
      delta_temp=(temp_max-temp_min)/temp_gn

      file23='      '
      call getenv('ftn23',file23)
      if (file23(1:4).ne.'    ') then
         OPEN(23,FILE=file23,STATUS='unknown',FORM='FORMATTED')
      else
         OPEN(23,FILE='dNg_over_dt.dat',STATUS='UNKNOWN')
      endif

      k=1
      do while(k.lt.t_gn+1)
         time=t_min+delta_t*k
         i=1
         write(22,1111) "timestep: ",k,"time: ",time
         write(23,1111) "timestep: ",k,"time: ",time
         do while(i.lt.temp_gn+1)
            temp=temp_min+delta_temp*i
            j=1
            do while(j.lt.HQener_gn+1)
               
               HQenergy=delta_HQener*j
c               region(1)=PI*temp/HQenergy
c               region(1)=sqrt(6d0*PI*fixAlphas)*temp/HQenergy
               max_dNg_over_dxdydt=0d0 

               call vegas8(region,ndim,dNg_over_dxdydt,init,ncall,
     &              itmx,nprn,tgral,sd,chi2a)
               dNg_over_dt(i,j)=tgral
               write(22,1112) time,temp,HQenergy,dNg_over_dt(i,j),
     &                        max_dNg_over_dxdydt
               write(23,1112) dNg_over_dt(i,j),max_dNg_over_dxdydt
               j=j+1            
               
            enddo
            i=i+1
         enddo
         write(6,*) "timestep: ",k
         k=k+1
      enddo

 1111 format(A10,2X,I3,2X,A6,2X,F8.3)
 1112 format(e10.4,4(2X,e10.4))
      
      write(6,*) "PROGRAM ENDS SUCCESSFULLY :)"

      end
         
