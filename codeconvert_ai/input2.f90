subroutine fit_uspline(n,xs,p,q,j,en,FF)
use pietc, only: o3
use pmat2, only: ldltb, ltdlbv
integer,                intent(in ):: n
real(dp),dimension(  n),intent(in ):: xs,p
real(dp),dimension(  n),intent(out):: q,j
real(dp),               intent(out):: en
logical,                intent(out):: FF

integer                   :: i,ip
real(dp)                  :: x,x2,sa,sb,ccc,xcmsx2
real(dp),dimension(n-1)   :: difp,sumq,cpp,cqp
real(dp),dimension(n,-1:0):: qq 

FF=F
if(n<1)stop 'In fit_uspline; size of data array must be positive'
if(n==1)then; q=0; j=0; en=0; return; endif

do i=2,n
   if(xs(i-1)>=xs(i)) then
      FF=T
      write(41,*) 'In fit_uspline; xs data must increase strictly monotonically'
      return
   end if
enddo

qq=0 

do i=1,n-1
   ip=i+1
   difp(i)=p(ip)-p(i)
   x2=xs(ip)-xs(i); x=o2*x2
   xcmsx2=o3*x**3*2


   ccc=2/x
   cpp(i)=u1/xcmsx2 
   cqp(i)=-difp(i)*x/xcmsx2 
   qq(i,0)=qq(i,0)+ccc; qq(ip,-1)=qq(ip,-1)+1/x; qq(ip,0)=qq(ip,0)+ccc
enddo
q(1:n-1)=-cqp; q(n)=0
q(2:n)=q(2:n)-cqp


call ldltb(n,1,qq) 
call ltdlbv(n,1,qq,q)
sumq=q(1:n-1)+q(2:n) 

en=o2*(dot_product(difp**2,cpp)+dot_product(sumq,cqp))

sb=0
do i=1,n-1
   ip=i+1
   x=o2*(xs(ip)-xs(i))
   xcmsx2=o3*x**3*2
   sa=(x*sumq(i)-difp(i))/xcmsx2
   j(i)=sa-sb
   sb  =sa
enddo
j(n)=-sb 
end subroutine fit_uspline
  
subroutine int_uspline(n,xs,p,q, m)
integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q
real(dp),dimension(n),intent(out):: m

real(dp),parameter:: u3o2=3*o2
real(dp):: a,b,c,d,e,t2,t3,t4,x,pa,pd,qa,qd
integer :: i,ip
e=u0
do i=1,n-1
   ip=i+1
   x=(xs(ip)-xs(i))*o2 
   t2=x*x/2
   t3=t2*x/3
   t4=t3*x/4
   pa=(p(ip)+p(i))*o2
   pd=(p(ip)-p(i))*o2/x
   qa=(q(ip)+q(i))*o2
   qd=(q(ip)-q(i))*o2/x
   c=qd
   a=pa-c*t2
   d=(qa-pd)*u3o2/t2
   b=qa-d*t2
   m(i)=e+a*x-b*t2+c*t3-d*t4
   e=e+2*(a*x+c*t3)
enddo
m(n)=e
end subroutine int_uspline

subroutine eval_tspline(n,xs,p,q, x,y)
integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q
real(dp),             intent(in ):: x
real(dp),             intent(out):: y

integer :: ia,ib
real(dp):: xr,xh,pm,qm,qah,qbh,qxh,qdh,shh,chh,sh,ch,xcmsh,shm,chm,shhm,chhm

if(x<=xs(1))then; xr=x-xs(1); y=p(1)+q(1)*expm( xr); return; endif
if(x>=xs(n))then; xr=x-xs(n); y=p(n)-q(n)*expm(-xr); return; endif
do ib=2,n
   if(xs(ib)<=xs(ib-1))cycle
   if(xs(ib)>=x)exit 
enddo
ia=ib-1
xh=(xs(ib)-xs(ia))*o2 
xr=x-xs(ia)-xh 
pm=(p(ia)+p(ib))*o2 
qm=(p(ib)-p(ia))/(xh*2) 
qah=q(ia)*o2; qbh=q(ib)*o2
qxh=qah+qbh-qm 
qdh=qbh-qah   
shh=sinh(xh);   chh=cosh(xh)
sh =sinh(xr);    ch=cosh(xr)
shm=sinhm(xr);  chm=coshm(xr)
shhm=sinhm(xh); chhm=coshm(xh)
xcmsh=xcms(xh)
qdh=qdh/shh; qxh=qxh/xcmsh 
y=pm+xr*qm +qdh*(chm-chhm) +  qxh*(xh*shm-xr*shhm)
end subroutine eval_tspline
   

subroutine eval_tsplined(n,xs,p,q, x,y,dydx)



integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q
real(dp),             intent(in ):: x
real(dp),             intent(out):: y,dydx

integer :: ia,ib
real(dp):: xr,xh,pm,qm,qah,qbh,qxh,qdh,shh,chh,sh,ch,xcmsh,shm,chm,shhm,chhm,&
           qemxr

if(x<=xs(1))then
   xr=x-xs(1); qemxr=q(1)*expm( xr); y=p(1)+qemxr; dydx=qemxr+q(1); return
endif
if(x>=xs(n))then
   xr=x-xs(n); qemxr=q(n)*expm(-xr); y=p(n)-qemxr; dydx=qemxr+q(n); return
endif
do ib=2,n
   if(xs(ib)<=xs(ib-1))cycle 
   if(xs(ib)>=x)exit    
enddo
ia=ib-1
xh=(xs(ib)-xs(ia))*o2 
xr=x-xs(ia)-xh          
pm=(p(ia)+p(ib))*o2   
qm=(p(ib)-p(ia))/(xh*2) 
qah=q(ia)*o2;  qbh=q(ib)*o2
qxh=qah+qbh-qm 
qdh=qbh-qah    
shh=sinh(xh);   chh=cosh(xh)
sh =sinh(xr);   ch=cosh(xr)
shm=sinhm(xr);  chm=coshm(xr)
shhm=sinhm(xh); chhm=coshm(xh)
xcmsh=xcms(xh)
qdh=qdh/shh; qxh=qxh/xcmsh 
y=pm+xr*qm +qdh*(chm-chhm)  +  qxh*(xh*shm-xr*shhm)
dydx=qm+qdh*sh +qxh*(xh*chm-shhm)
end subroutine eval_tsplined


subroutine eval_tsplinedd(n,xs,p,q, x,y,dydx,ddydxx)



integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q
real(dp),             intent(in ):: x
real(dp),             intent(out):: y,dydx,ddydxx

integer :: ia,ib
real(dp):: xr,xh,pm,qm,qah,qbh,qxh,qdh,shh,chh,sh,ch,xcmsh,shm,chm,shhm,chhm,&
           qemxr

if(x<=xs(1))then
   xr=x-xs(1); qemxr=q(1)*expm( xr); y=p(1)+qemxr; dydx=qemxr+q(1)
   ddydxx=dydx; return
endif
if(x>=xs(n))then
   xr=x-xs(n); qemxr=q(n)*expm(-xr); y=p(n)-qemxr; dydx=qemxr+q(n)
   ddydxx=-dydx; return
endif
do ib=2,n
   if(xs(ib)<=xs(ib-1))cycle 
   if(xs(ib)>=x)exit    
enddo
ia=ib-1
xh=(xs(ib)-xs(ia))*o2 
xr=x-xs(ia)-xh          
pm=(p(ia)+p(ib))*o2   
qm=(p(ib)-p(ia))/(xh*2) 
qah=q(ia)*o2;  qbh=q(ib)*o2
qxh=qah+qbh-qm 
qdh=qbh-qah    
shh=sinh(xh);   chh=cosh(xh)
sh =sinh(xr);   ch=cosh(xr)
shm=sinhm(xr);  chm=coshm(xr)
shhm=sinhm(xh); chhm=coshm(xh)
xcmsh=xcms(xh)
qdh=qdh/shh; qxh=qxh/xcmsh 
y=pm+xr*qm +qdh*(chm-chhm)  +  qxh*(xh*shm-xr*shhm)
dydx=qm+qdh*sh +qxh*(xh*chm-shhm)
ddydxx=qdh*ch +qxh*xh*sh
end subroutine eval_tsplinedd


subroutine eval_tsplineddd(n,xs,p,q, x,y,dydx,ddydxx,dddydxxx)



integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q
real(dp),             intent(in ):: x
real(dp),             intent(out):: y,dydx,ddydxx,dddydxxx

integer :: ia,ib
real(dp):: xr,xh,pm,qm,qah,qbh,qxh,qdh,shh,chh,sh,ch,xcmsh,shm,chm,shhm,chhm,&
           qemxr

if(x<=xs(1))then
   xr=x-xs(1); qemxr=q(1)*expm( xr); y=p(1)+qemxr; dydx=qemxr+q(1)
   ddydxx=dydx; dddydxxx=dydx; return
endif
if(x>=xs(n))then
   xr=x-xs(n); qemxr=q(n)*expm(-xr); y=p(n)-qemxr; dydx=qemxr+q(n)
   ddydxx=-dydx; dddydxxx=dydx; return
endif
do ib=2,n
   if(xs(ib)<=xs(ib-1))cycle 
   if(xs(ib)>=x)exit    
enddo
ia=ib-1
xh=(xs(ib)-xs(ia))*o2 
xr=x-xs(ia)-xh          
pm=(p(ia)+p(ib))*o2   
qm=(p(ib)-p(ia))/(xh*2) 
qah=q(ia)*o2;  qbh=q(ib)*o2
qxh=qah+qbh-qm 
qdh=qbh-qah    
shh=sinh(xh);   chh=cosh(xh)
sh =sinh(xr);   ch=cosh(xr)
shm=sinhm(xr);  chm=coshm(xr)
shhm=sinhm(xh); chhm=coshm(xh)
xcmsh=xcms(xh)
qdh=qdh/shh; qxh=qxh/xcmsh 
y       =pm+xr*qm +qdh*(chm-chhm)  +  qxh*(xh*shm-xr*shhm)
dydx    =      qm +qdh*sh          +  qxh*(xh*chm-   shhm)
ddydxx  =          qdh*ch          +  qxh* xh*sh
dddydxxx=          qdh*sh          +  qxh* xh*ch
end subroutine eval_tsplineddd


subroutine eval_itspline(n,xs, p,q,m,  x,y)



integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q,m
real(dp),             intent(in ):: x
real(dp),             intent(out):: y

real(dp):: a,b,c,d,t2,xh,shx,chmx,shmx,chmmx,xcmsx,xr,pa,pd,qa,qd
integer :: ia,ib

if(x<=xs(1))then; xr=x-xs(1); y=     p(1)*xr+q(1)*expmm( xr); return; endif
if(x>=xs(n))then; xr=x-xs(n); y=m(n)+p(n)*xr+q(n)*expmm(-xr); return; endif
do ib=2,n
   if(xs(ib)<=xs(ib-1))cycle 
   if(xs(ib)>=x)exit 
enddo
ia=ib-1
xh=(xs(ib)-xs(ia))*o2   
shx =sinh (xh)
chmx=coshm(xh)
xcmsx=xcms(xh)
xr=x-xs(ia)-xh          
pa=(p(ib)+p(ia))*o2
pd=(p(ib)-p(ia))*o2/xh
qa=(q(ib)+q(ia))*o2
qd=(q(ib)-q(ia))*o2/shx

c=qd
a=pa-c*chmx
d=(qa-pd)*xh/xcmsx
b=qa-d*chmx

t2=xr**2/2
shmx =sinhm (xr)
chmmx=coshmm(xr)
y=m(ia)+a*xr+b*t2+c*shmx+d*chmmx
end subroutine eval_itspline


subroutine eval_uspline(n,xs,p,q, x,y)









integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q
real(dp),             intent(in ):: x
real(dp),             intent(out):: y

integer :: ia,ib
real(dp):: xr,xh,pm,qm,qah,qbh,qxh,qdh,shh,chh,sh,ch,xcmsh,shm,chm,shhm,chhm

if(x<=xs(1))then; xr=x-xs(1); y=p(1)+q(1)*xr; return; endif
if(x>=xs(n))then; xr=x-xs(n); y=p(n)+q(n)*xr; return; endif
do ib=2,n
   if(xs(ib)<=xs(ib-1))cycle 
   if(xs(ib)>=x)exit 
enddo
ia=ib-1
xh=(xs(ib)-xs(ia))*o2   
xr=x-xs(ia)-xh          
pm=(p(ia)+p(ib))*o2     
qm=(p(ib)-p(ia))/(xh*2) 
qah=q(ia)*o2;  qbh=q(ib)*o2
qxh=qah+qbh-qm 
qdh=qbh-qah    
shh=xh;        chh=u1
sh =xr;        ch =u1
shm =xr**3/6;  chm =xr**2*o2
shhm=xh**3/6;  chhm=xh**2*o2
xcmsh=xh**3/3
qdh=qdh/shh; qxh=qxh/xcmsh 
y=pm+xr*qm +qdh*(chm-chhm) +  qxh*(xh*shm-xr*shhm)

end subroutine eval_uspline
   

