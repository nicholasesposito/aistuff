function expm(x) result(e)!
real(dp),intent(in ):: x
real(dp)            :: e
real(dp):: p
integer :: i
if(abs(x)>o2)then
   e=exp(x)-u1
else
   p=x; e=p
   do i=2,19; p=p*x/i; e=e+p; if(abs(p)<=abs(e*eps))return; enddo
endif
end function expm


function expmm(x) result(e)!                                           [expmm]
real(dp),intent(in ):: x
real(dp)            :: e

real(dp):: p
integer :: i

if(abs(x)>o2)then
   e=exp(x)-u1-x
else
   p=x*x*o2; e=p
   do i=3,25; p=p*x/i; e=e+p; if(abs(p)<=abs(e*eps))return; enddo
endif
end function expmm

function coshm(x) result(c)
real(dp),intent(in ):: x
real(dp)            :: c

c=2*sinh(x*o2)**2
end function coshm


function sinhm(x) result(s)
real(dp),intent(in ):: x
real(dp)            :: s

real(dp):: p,xx
integer :: i

if(abs(x)>o2)then
   s=sinh(x)-x
else
   p=x**3/6;  s=p;  xx=x*x
   do i=5,19,2; p=p*xx/(i*(i-1)); s=s+p; if(abs(p)<=abs(s*eps))return; enddo
endif
end function sinhm

function coshmm(x) result(c)
real(dp),intent(in ):: x
real(dp)            :: c

real(dp)            :: xh

xh=x*o2
c=sinhm(xh)*(2*sinh(xh)+x)
end function coshmm


function xcms(x) result(e)

real(dp),intent(in ):: x
real(dp)            :: e

real(dp):: p,xx
integer :: i,i2

if(abs(x)>o2)then
   e=x*coshm(x)-sinhm(x)
else
   p=x**3/3;  e=p;  xx=x*x
   do i=2,15
      i2=i*2; p=p*xx/(i2*(i2+1)); e=e+i*p; if(abs(p)<=abs(e*eps))return
   enddo
endif
end function xcms


function enbase_t(tspan,hspan)result(r)
real(dp),intent(in ):: tspan,hspan
real(dp)            :: r

if(tspan<u0)stop 'In enbase_t; thspan must be positive'
if(hspan==u0)then; r=u1; return; endif
r=hspan**2/expmm(-tspan)*o2
end function enbase_t


subroutine tbnewton(nh,m,bigT,halfgate,hgts,hs,hgtp,p,q, te,dhdt,FF)
integer,               intent(in ):: nh,m
real(dp),              intent(in ):: bigT,halfgate
integer ,dimension(nh),intent(in ):: hgts
real(dp),dimension(nh),intent(in ):: hs
integer, dimension(m), intent(in ):: hgtp
real(dp),dimension(m) ,intent(in ):: p,q
real(dp),dimension(nh),intent(out):: dhdt, te
logical,               intent(out):: FF

integer,parameter    :: nit=12
real(dp),dimension(m):: tr
real(dp)             :: gate,tee,he,hac,dhadt,dh,dt
integer              :: i,it

gate=2*halfgate/bigT
tr=hgtp*halfgate/bigT
do i=1,nh
   tee=hgts(i)*halfgate/bigT
   he=hs(i)
!  Use Newton iterations to estimate the rescaled time, tee, at which the
!  height is he
   it = 1
   do while (it <= nit)
      call eval_tspline(m,tr,p,q, tee,hac,dhadt)
      if(it==1)dhdt(i)=dhadt/bigT
      if(dhadt==u0)exit
      dh=hac-he
      dt=-dh/dhadt
      if(abs(dt)>gate)then
         write(41,*) 'WARNING! In tbnewton; i,it,dt/gate = ',i,it,dt/gate
         exit
      endif
      if(abs(dh)<heps)then
         dhdt(i)=dhadt/bigT
         exit
      endif
      tee=tee+dt
      it = it + 1
   enddo
   FF=(it>nit)
   if(FF)then
      write(41,'("In tbnewton; Newton iterations seem not to be")')
      write(41,'("converging at i=",i3)'),i
      write(41,'("tee,he,hac,heps,dhadt:",5(1x,e11.4))'),tee,he,hac,heps,dhadt
   endif
   te(i) = tee
enddo
end subroutine tbnewton


subroutine ubnewton(nh,m,halfgate,hgts,hs,hgtp,p,q, te,dhdt,FF)
integer,               intent(in ):: nh,m
real(dp),              intent(in ):: halfgate
integer, dimension(nh),intent(in ):: hgts
real(dp),dimension(nh),intent(in ):: hs
integer, dimension(m), intent(in ):: hgtp
real(dp),dimension(m) ,intent(in ):: p,q
real(dp),dimension(nh),intent(out):: dhdt, te
logical,               intent(out):: FF

integer,parameter    :: nit=12
real(dp),dimension(m):: tr
real(dp)             :: gate,tee,he,hac,dhadt,dh,dt
integer              :: i,it

gate=2*halfgate
tr=hgtp*halfgate
do i=1,nh
   tee=hgts(i)*halfgate
   he=hs(i)
!  Use Newton iterations to estimate the rescaled time, tee, at which the
!  height is he
   it = 1
   do while (it <= nit)
      call eval_uspline(m,tr,p,q, tee,hac,dhadt)
      if(it==1)dhdt(i)=dhadt
      if(dhadt==u0)exit
      dh=hac-he
      dt=-dh/dhadt
      if(abs(dt)>gate)then
         write(41,*) 'WARNING! In ubnewton; i,it,dt/gate = ',i,it,dt/gate
         exit
      endif
      if(abs(dh)<heps)then
         dhdt(i)=dhadt
         exit
      endif
      tee=tee+dt
      it = it + 1
   enddo
   FF=(it>nit)
   if(FF)then
      write(41,'("In ubnewton; Newton iterations seem not to be")')
      write(41,'("converging at i=",i3)'),i
      write(41,'("tee,he,hac,heps,dhadt:",5(1x,e11.4))'),tee,he,hac,heps,dhadt
   endif
   te(i) = tee
enddo
end subroutine ubnewton


subroutine fit_gtspline(n,xs,ys,on,q,j,yac,en,FF)
integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,ys
logical, dimension(n),intent(in ):: on
real(dp),dimension(n),intent(out):: q,j,yac
real(dp),             intent(out):: en
logical,              intent(out):: FF

real(dp),dimension(n):: xa,ya,qa,ja
integer              :: i,k,m

m=0
do i=1,n
   if(on(i))then; m=m+1; xa(m)=xs(i); ya(m)=ys(i); endif
enddo
call fit_tspline(m,xa(1:m),ya(1:m),qa(1:m),ja(1:m),en,FF)
if(FF)then
   write(41,*) 'In fit_gtspline; failure flag raised at call to fit_tspline'
   return
endif
k=0
do i=1,n
   if(on(i))then
      k=k+1
      q(i)=qa(k)
      j(i)=ja(k)
      yac(i)=ys(i)
   else
      call eval_tsplined(m,xa(1:m),ya(1:m),qa(1:m),xs(i), yac(i),q(i))
      j(i)=0
   endif
enddo
end subroutine fit_gtspline


subroutine fit_tspline(n,xs,p,q,j,en,FF)
use pmat2, only: ldltb, ltdlbv
integer,                intent(in ):: n
real(dp),dimension(  n),intent(in ):: xs,p
real(dp),dimension(  n),intent(out):: q,j
real(dp),               intent(out):: en
logical,                intent(out):: FF
!----------------------------------------------------------------------------
integer                   :: i,ip
real(dp)                  :: x,ch,sh,sa,sb,sap,ccc,xcmsx2,egg,ehh
real(dp),dimension(n-1)   :: difp,sumq,cpp,cqp
real(dp),dimension(n,-1:0):: qq ! <- Tridiagonal, stored as rows of nonupper
!=============================================================================
FF=F
if(n<1)stop 'In fit_tspline; size of data array must be positive'
if(n==1)then; q=0; j=0; en=0; return; endif
! apply a strict monotonicity check on the xs:
do i=2,n
   if(xs(i-1)>=xs(i)) then
      FF=T
      write(41,*) 'In fit_tspline; xs data must increase strictly monotonically'
      return
   end if
enddo
qq=0
do i=1,n-1
   ip=i+1
   difp(i)=p(ip)-p(i)
   x=(xs(ip)-xs(i))*o2 
   ch=cosh(x);  sh=sinh(x)
   xcmsx2=xcms(x)*2
   egg=x*sh/xcmsx2; ehh=ch/(2*sh)
   ccc=egg+ehh
   cpp(i)=ch/xcmsx2 
   cqp(i)=-difp(i)*sh/xcmsx2 
   qq(i,0)=qq(i,0)+ccc; qq(ip,-1)=qq(ip,-1)+egg-ehh; qq(ip,0)=qq(ip,0)+ccc
enddo
qq(1,0)=qq(1,0)+1
qq(n,0)=qq(n,0)+1

q(1:n-1)=-cqp; q(n)=0
q(2:n)=q(2:n)-cqp

call ldltb(n,1,qq) 
call ltdlbv(n,1,qq,q)
sumq=q(1:n-1)+q(2:n) 

en=o2*(dot_product(difp**2,cpp)+dot_product(sumq,cqp))
sb=q(1)
do i=1,n-1
   ip=i+1
   x=o2*(xs(ip)-xs(i))
   xcmsx2=xcms(x)*2
   ch=cosh(x);  sh=sinh(x)
   sap=(sh*sumq(i)-ch*difp(i))/xcmsx2
   sa=sap+q(i)
   j(i)=sa-sb
   sb  =sap+q(ip)
enddo
j(n)=q(n)-sb 
end subroutine fit_tspline


subroutine int_tspline(n,xs,p,q, m)
integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q
real(dp),dimension(n),intent(out):: m
!-----------------------------------------------------------------------------
real(dp):: a,b,c,d,e,t2,x,pa,pd,qa,qd,shx,chmx,shmx,chmmx,xcmsx
integer :: i,ip


e=u0
do i=1,n-1
   ip=i+1
   x=(xs(ip)-xs(i))*o2 !<- interval half-width
   t2=x*x*o2
   shx  =sinh  (x)
   chmx =coshm (x)
   shmx =sinhm (x)
   chmmx=coshmm(x)
   xcmsx=xcms  (x)
   pa=(p(ip)+p(i))*o2
   pd=(p(ip)-p(i))*o2/x
   qa=(q(ip)+q(i))*o2
   qd=(q(ip)-q(i))*o2/shx
   c=qd
   a=pa-c*chmx
   d=(qa-pd)*x/xcmsx
   b=qa-d*chmx
   m(i)=e+a*x -b*t2 +c*shmx -d*chmmx
   e=e+2*(a*x+c*shmx)
enddo
m(n)=e
end subroutine int_tspline


subroutine fit_guspline(n,xs,ys,on,q,j,yac,en,FF)
integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,ys
logical, dimension(n),intent(in ):: on
real(dp),dimension(n),intent(out):: q,j,yac
real(dp),             intent(out):: en
logical,              intent(out):: FF

real(dp),dimension(n):: xa,ya,qa,ja
integer              :: i,k,m

m=0
do i=1,n
   if(on(i))then; m=m+1; xa(m)=xs(i); ya(m)=ys(i); endif
enddo
call fit_uspline(m,xa(1:m),ya(1:m),qa(1:m),ja(1:m),en,FF)
if(FF)then
   write(41,*) 'In fit_guspline; failure flag raised at call to fit_uspline'
   return
endif
k=0
do i=1,n
   if(on(i))then
      k=k+1
      q(i)=qa(k)
      j(i)=ja(k)
      yac(i)=ys(i)
   else
      call eval_usplined(m,xa(1:m),ya(1:m),qa(1:m),xs(i), yac(i),q(i))
      j(i)=0
   endif
enddo
end subroutine fit_guspline


subroutine fit_uspline(n,xs,p,q,j,en,FF)

