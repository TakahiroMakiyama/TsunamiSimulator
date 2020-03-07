!…[‚ªˆê’èê‡‚Ì’Ã”g‚Ì‰^“®‚ğ·•ª‹ß—®‚Å‰ğ‚­
implicit none
!•Ï”‚ÌŒ^éŒ¾ --------------------------------------
integer, parameter :: iq=500
real*8 ::   uf(iq),  up(iq),  ub(iq)
real*8 ::   zf(iq+1),zp(iq+1),zb(iq+1)
real*8 ::   depth,width,grav,time_end,time_out,dx,dt,eps
integer ::  i,n,nend,nout,index
! ƒpƒ‰ƒ[ƒ^‚Ìİ’è -------------------------
depth=100.           !…[ (m)
width=100.*1000.     !•(m)
grav=9.8             !d—Í‰Á‘¬“x(m/s^2)
!------------------------------
time_end=2.*60.*60.  !ŒvZŠÔ(sec)
time_out=60.         !ƒf[ƒ^o—ÍŠÔŠÔŠu(sec)
!----------------------------
dx=width/iq      !‹óŠÔƒOƒŠƒbƒhƒTƒCƒY
dt=1             !ŠÔƒXƒeƒbƒvƒTƒCƒY
do i=1,iq 
 up(i)=0.
end do
do i=1,iq+1
 zp(i)=exp(-dble(i-iq/2)**2/dble(iq/30)**2)
end do
!----o—Íƒf[ƒ^ƒtƒ@ƒCƒ‹‚ğw’è--------
open(10,file='z.data')
open(20,file='u.data')
!------------------------------------
nend=time_end/dt  !ŠÔƒXƒeƒbƒv”
nout=time_out/dt  !o—Íƒf[ƒ^ƒXƒeƒbƒvŠÔŠu
index=0
! ŠÔƒ‹[ƒv(n=0,nend‚Ü‚ÅŒJ‚è•Ô‚·) 
do n=0,nend
!**************************************************************
! Å‰‚ÌƒXƒeƒbƒv‚Ì‚İŒ»İƒXƒeƒbƒv’l‚©‚ç–¢—ˆƒXƒeƒbƒv’l‚ğŒvZi‘O•û·•ªj
    if(n==0) then
       do i=1,iq
         uf(i)=up(i)-grav*(dt/dx)*(zp(i+1)-zp(i))
       end do
       do i=2,iq
        zf(i)=zp(i)-depth*(dt/dx)*(up(i)-up(i-1))
       end do
    end if
! ‰ß‹ƒXƒeƒbƒv’l‚ÆŒ»İƒXƒeƒbƒv’l‚©‚ç–¢—ˆƒXƒeƒbƒv’l‚ğŒvZi’†‰›·•ª)
    if(n>=1) then
      do i=1,iq
        uf(i)=ub(i)-2*grav*(dt/dx)*(zp(i+1)-zp(i))
      end do
      do i=2,iq
        zf(i)=zb(i)-2*depth*(dt/dx)*(up(i)-up(i-1))
      end do
    end if
!---‹«ŠEğŒ----------
    uf(1 )=0.
    uf(iq)=0.
! ŒvZ‚ÌˆÀ’è‰»‚Ì‚½‚ß‚Ì‚¨‚Ü‚¶‚È‚¢iAsselin filterj-------
    if(n>=1) then
      eps=0.01
      do i=1,iq
        up(i)=up(i)+eps*(uf(i)-2*up(i)+ub(i))
      end do
      do i=1,iq+1
        zp(i)=zp(i)+eps*(zf(i)-2*zp(i)+zb(i))
      end do
    end if
 !---noutã‚¹ãƒ†ãƒƒãƒ—ã”ã¨ã«ãƒ‡ãƒ¼ã‚¿ã‚’ãƒ•ã‚¡ã‚¤ãƒ«ã«å‡ºåŠ› --------------
    if(mod(n,nout).eq.0) then
       do i=2,iq
         write(10,*) dx*(i-2)/1000.,zp(i)
         write(20,*) dx*(i-2)/1000.,up(i)
       end do
       write(10,*)
       write(10,*)
       write(20,*)
       write(20,*)
       write(*,*) 'time (sec)=',dt*n,index
       index=index+1
    end if
! ƒf[ƒ^‚ğƒtƒ@ƒCƒ‹o—Í(noutƒXƒeƒbƒv‚²‚Æ‚É) ------------------
    do i=1,iq
      ub(i)=up(i)
    end do
    do i=1,iq+1
      zb(i)=zp(i)
    end do
    do i=1,iq
      up(i)=uf(i)
    end do
    do i=1,iq+1
      zp(i)=zf(i)
    end do
!*****************************************************************
end do
stop
end
   

