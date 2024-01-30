!=============================================================================!
!=============================================================================!
!=============================================================================!
module spectra
  implicit none
  save
  integer, parameter :: n1=42,n2=140,n3=19,n4=22,n5=21,n6=4
  integer :: sym=1
  double precision, parameter :: m_p = 938.272d6        ! proton mass/eV 
  double precision, parameter :: E0=1.d10,d=1.5848,E_min=1d5,dn=0.1d0
  double precision xa(n1),ya(n2),fqgs(n1,n2,n6),elqgs(n1,n2,n6),nuqgs(n1,n2,n6,n6)
  double precision d_f,x0
  double precision, parameter :: Etrans=40.d9,width=0.01d0
end module spectra
!=============================================================================!
!=============================================================================!
!program secondary_spectra
!  use spectra, only : sym,m_p 
!  implicit none
!  integer i,k
!  double precision E_p,E_s,T,s1,s2,s3
!  double precision spec_gam,spec_gam_QGS,spec_el_QGS,spec_nu_tot_QGS
!  
!  call init
!
! input:  E_p/eV of primary, energy E_s/eV of secondary
!         reac: 1 pp, 2 pHe, 3 Hep, 4 HeHe (nuclei not implemented)
!         sym = 1, Kamae symmetrised, sym=0 not
! output: E_s*dsigma(E_p,E_s)/dE_s in mbarn 
! spec_gam is a mix of Kamae (low-energies) and QGSJET (high-energies)
!
!  sym=1
!  k=1
!  
!  E_p = 2.d9
!  open(21,file='test1')
!  open(22,file='test2')                            
!  do i=1,20000,100
!     E_s = 0.9979**i * E_p
!     s1 = spec_nu_tot_QGS(E_p,E_s,k)
!     s2 = spec_gam_QGS(E_p,E_s,k)
!     s3 = spec_el_QGS(E_p,E_s,k)
!     write(21,*) real(E_s),real(s1),real(s2),real(s3)
!     s2 = spec_gam(E_p,E_s,k)
!     write(22,*) real(E_s),real(s2)
!  end do
!  close(21)
!  close(22)
!
!end program secondary_spectra
!=============================================================================!
!=============================================================================!
double precision function spec_gam(E_p,E_g,reac)
  use spectra
  implicit none
  integer reac
  double precision :: E_p,E_g,E,x,s,q,k
  double precision, parameter :: Esym=67.285d6
  double precision kamae_nd,spec_gam_QGS

  if (E_g<Esym.and.sym==1) then
     E = Esym**2/E_g
  else
     E = E_g
  end if
  
  q = spec_gam_QGS(E_p,E_g,reac) /E_g

  select case (reac)
  case (1)                          ! pp
     x = log10(E_p)
     s = (tanh((x-x0)/(width*x0))+1.d0)/2.d0
     k = kamae_nd(E_p/1.d9,E/1.d9) /E 
  case (2)                          ! p+He
     x = log10(E_p)
     s = (tanh((x-x0)/(width*x0))+1.d0)/2.d0
     k = 3.03d0 * kamae_nd(E_p/1.d9,E/1.d9)/E 
  case (3)                          ! He+p
     x = log10(E_p/4.d0)
     s = (tanh((x-x0)/(width*x0))+1.d0)/2.d0
     k = 4.d0 * kamae_nd(E_p/4.d9,E/1.d9)/E  
  case (4)                          ! He+He
     x = log10(E_p/4.d0)
     s = (tanh((x-x0)/(width*x0))+1.d0)/2.d0
     k = 12.1d0 * kamae_nd(E_p/4.d9,E/1.d9)/E 
  end select

  spec_gam = k*(1.d0-s) +q*s
  spec_gam = spec_gam * E_g

end function spec_gam
!=============================================================================!
!=============================================================================!
double precision function spec_gam_QGS(E_p,E_g,k)
  use spectra
  implicit none
  integer i,j,k
  double precision E_p,E_g
  double precision x,y,y1,y2,y3,y4,t,u,r,iex,jex

  if (E_p<1.5d10.or.E_g<1.d6) then
     spec_gam_QGS = 0.d0
     return
  end if

  x=log10(E_p)
  y=log10(E_g)

  iex=log(E_p/E0)/log(d)
  jex=(y-d_f)/dn
  i=int(iex)
  j=int(jex)
  jex=jex-j
  iex=iex-i

  y1 = fqgs(i,j,k)
  y2 = fqgs(i+1,j,k)
  y4 = fqgs(i,j+1,k)
  y3 = fqgs(i+1,j+1,k)

!corrected (may intruduce rounding errors)
!  t = (x-xa(i))/(xa(i+1)-xa(i))
!  u = (y-ya(j))/(ya(j+1)-ya(j))
!  r = (1.d0-t)*(1.d0-u)*y1 + t*(1.d0-u)*y2 + t*u*y3  + (1.d0-t)*u*y4
  r = (1.d0-iex)*(1.d0-jex)*y1 + iex*(1.d0-jex)*y2 +(1.d0-iex)*jex*y4 +iex*jex*y3;
  spec_gam_QGS = 10.d0**r 

end function spec_gam_QGS
!=============================================================================!
!=============================================================================!
double precision function spec_el_QGS(E_p,E_el,k)
  use spectra
  implicit none
  integer i,j,k
  double precision E_p,E_el
  double precision x,y,y1,y2,y3,y4,t,u,r,iex,jex

  if (E_p<1.5d10.or.E_el<1.d6) then
     spec_el_QGS = 0.d0
     return
  end if

  x=log10(E_p)
  y=log10(E_el)

  iex=log(E_p/E0)/log(d)
  jex=(y-d_f)/dn
  i=int(iex)
  j=int(jex)
  jex=jex-j
  iex=iex-i

  y1 = elqgs(i,j,k)
  y2 = elqgs(i+1,j,k)
  y4 = elqgs(i,j+1,k)
  y3 = elqgs(i+1,j+1,k)

!corrected (may intruduce rounding errors)
!  t = (x-xa(i))/(xa(i+1)-xa(i))
!  u = (y-ya(j))/(ya(j+1)-ya(j))
!  r = (1.d0-t)*(1.d0-u)*y1 + t*(1.d0-u)*y2 + t*u*y3  + (1.d0-t)*u*y4
  r = (1.d0-iex)*(1.d0-jex)*y1 + iex*(1.d0-jex)*y2 +(1.d0-iex)*jex*y4 +iex*jex*y3;

  spec_el_QGS = 10.d0**r 

end function spec_el_QGS
!=============================================================================!
!=============================================================================!
double precision function spec_nu_tot_QGS(E_p,E_nu,k)
  use spectra
  implicit none
  integer i,j,fl,k
  double precision E_p,E_nu
  double precision x,y,y1,y2,y3,y4,t,u,r,s,iex,jex

  if (E_p<1.5d10.or.E_nu<1.d8) then
     spec_nu_tot_QGS = 0.d0
     return
  end if

  x=log10(E_p)
  y=log10(E_nu)

  iex=log(E_p/E0)/log(d)
  jex=(y-d_f)/dn
  i=int(iex)
  j=int(jex)
  jex=jex-j
  iex=iex-i

  s=0.d0

  do fl=1,4

     y1 = nuqgs(i,j,fl,k)
     y2 = nuqgs(i+1,j,fl,k)
     y4 = nuqgs(i,j+1,fl,k)
     y3 = nuqgs(i+1,j+1,fl,k)

     !t = (x-xa(i))/(xa(i+1)-xa(i))!may intruduce rounding errors
     !u = (y-ya(j))/(ya(j+1)-ya(j))!may intruduce rounding errors
     !r = (1.d0-t)*(1.d0-u)*y1 + t*(1.d0-u)*y2 + t*u*y3  + (1.d0-t)*u*y4
     r = (1.d0-iex)*(1.d0-jex)*y1 + iex*(1.d0-jex)*y2 +(1.d0-iex)*jex*y4 +iex*jex*y3;

     s = s + 10.d0**r 

  end do

  spec_nu_tot_QGS = s

end function spec_nu_tot_QGS
!=============================================================================!
!=============================================================================!
!=============================================================================!
!=============================================================================!
! fl = 1 electron neutrino
! fl = 2 electron antineutrino
! fl = 3 muon neutrino
! fl = 4 muon antineutrino
double precision function spec_nu_QGS(E_p,E_nu,k,fl)
  use spectra
  implicit none
  integer i,j,fl,k
  double precision E_p,E_nu
  double precision x,y,y1,y2,y3,y4,r,iex,jex

  if (E_p<1.5d10.or.E_nu<1.d8) then
     spec_nu_QGS = 0.d0
     return
  end if

  x=log10(E_p)
  y=log10(E_nu)

  iex=log(E_p/E0)/log(d)
  jex=(y-d_f)/dn
  i=int(iex)
  j=int(jex)
  jex=jex-j
  iex=iex-i

  y1 = nuqgs(i,j,fl,k)
  y2 = nuqgs(i+1,j,fl,k)
  y4 = nuqgs(i,j+1,fl,k)
  y3 = nuqgs(i+1,j+1,fl,k)

  r = (1.d0-iex)*(1.d0-jex)*y1 + iex*(1.d0-jex)*y2 +(1.d0-iex)*jex*y4 +iex*jex*y3;

  spec_nu_QGS =  10.d0**r

end function spec_nu_QGS
!=============================================================================!
!=============================================================================!
subroutine pp_init
  use spectra
  implicit none
  integer i,j 
  double precision x,y,f,f1,f2,f3,f4
  character*256 c_line
  
  call banner
   
  d_f = log10(E_min)-0.1d0         
  x0 = log10(Etrans)

  open(21,file='tables/pp/cross_sec/gam_pp')
  do i=1,n1
     do j=1,n2
        read(21,*) x,y,f
        xa(i)=log10(x)
        ya(j)=log10(y)
        if (f>0.d0) then
           fqgs(i,j,1)=log10(f)
        else
           fqgs(i,j,1)=-20.d0
        end if
     end do
  end do
!  write(*,*) x,y,f
  close(21)


  
  open(21,file='tables/pp/cross_sec/el_pp')
  do i=1,n1
     do j=1,n2
        read(21,*) x,y,f
        xa(i)=log10(x)
        ya(j)=log10(y)
        if (f>0.d0) then
           elqgs(i,j,1)=log10(f)
        else
           elqgs(i,j,1)=-20.d0
        end if
     end do
     read(21,*) 
  end do
!  write(*,*) x,y,f
  close(21)

  open(21,file='tables/pp/cross_sec/nu_pp')
  do i=1,n1
     do j=1,n2
        read(21,*) x,y,f1,f2,f3,f4
        xa(i)=log10(x)
        ya(j)=log10(y)
        if (f1>0.d0) then
           nuqgs(i,j,1,1)=log10(f1)
        else
           nuqgs(i,j,1,1)=-20.d0
        end if
        if (f1>0.d0) then
           nuqgs(i,j,2,1)=log10(f2)
        else
           nuqgs(i,j,2,1)=-20.d0
        end if
        if (f1>0.d0) then
           nuqgs(i,j,3,1)=log10(f3)
        else
           nuqgs(i,j,3,1)=-20.d0
        end if
        if (f1>0.d0) then
           nuqgs(i,j,4,1)=log10(f4)
        else
           nuqgs(i,j,4,1)=-20.d0
        end if
     end do
     read(21,*) 
  end do
!  write(*,*) x,y,f1,f2,f3,f4
  close(21)

  write(*,*) 
  write(*,*) ' tables (for pp) read'
  write(*,*) 
  write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,*) 

end subroutine pp_init
!=============================================================================!
!=============================================================================!
subroutine banner
  write(*,*) 
  write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,*) '!  ppfrag v3    --  for more info see                    !'
  write(*,*) '!  http://sourceforge.net/projects/ppfrag                !'
  write(*,*) '!  Authors: M.Kachelriess, S. Ostapchenko                !'
  write(*,*) '!                                                        !'
  write(*,*) '!  if you use these data, please cite                    !'
  write(*,*) '!  Refs.: * MK, SO, Phys.Rev. D83 (2011) 014018          !'     
  write(*,*) '!           arxiv[1010.1869]                             !'
  write(*,*) '!         * SO, Phys.Rev. D86 (2012) 043004              !'
  write(*,*) '!           arxiv[1206.4705]                             !'
  write(*,*) '!                                                        !'
  write(*,*) '!  if you use photon spectra with Kamae at E<Etrans,     !'
  write(*,*) '!  cite please additionally                              !'
  write(*,*) '!  T. Kamae et al., ApJ 647 (2006) 692.                  !'     
  write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
end subroutine banner
!=============================================================================!
!=============================================================================!
