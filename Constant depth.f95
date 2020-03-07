!水深が一定場合の津波の運動を差分近似式で解く
implicit none
!変数の型宣言 --------------------------------------
integer, parameter :: iq=500
real*8 ::   uf(iq),  up(iq),  ub(iq)
real*8 ::   zf(iq+1),zp(iq+1),zb(iq+1)
real*8 ::   depth,width,grav,time_end,time_out,dx,dt,eps
integer ::  i,n,nend,nout,index
! パラメータの設定 -------------------------
depth=100.           !水深 (m)
width=100.*1000.     !幅(m)
grav=9.8             !重力加速度(m/s^2)
!------------------------------
time_end=2.*60.*60.  !計算時間(sec)
time_out=60.         !データ出力時間間隔(sec)
!----------------------------
dx=width/iq      !空間グリッドサイズ
dt=1             !時間ステップサイズ
do i=1,iq 
 up(i)=0.
end do
do i=1,iq+1
 zp(i)=exp(-dble(i-iq/2)**2/dble(iq/30)**2)
end do
!----出力データファイルを指定--------
open(10,file='z.data')
open(20,file='u.data')
!------------------------------------
nend=time_end/dt  !時間ステップ数
nout=time_out/dt  !出力データステップ間隔
index=0
! 時間ループ(n=0,nendまで繰り返す) 
do n=0,nend
!**************************************************************
! 最初のステップのみ現在ステップ値から未来ステップ値を計算（前方差分）
    if(n==0) then
       do i=1,iq
         uf(i)=up(i)-grav*(dt/dx)*(zp(i+1)-zp(i))
       end do
       do i=2,iq
        zf(i)=zp(i)-depth*(dt/dx)*(up(i)-up(i-1))
       end do
    end if
! 過去ステップ値と現在ステップ値から未来ステップ値を計算（中央差分)
    if(n>=1) then
      do i=1,iq
        uf(i)=ub(i)-2*grav*(dt/dx)*(zp(i+1)-zp(i))
      end do
      do i=2,iq
        zf(i)=zb(i)-2*depth*(dt/dx)*(up(i)-up(i-1))
      end do
    end if
!---境界条件----------
    uf(1 )=0.
    uf(iq)=0.
! 計算の安定化のためのおまじない（Asselin filter）-------
    if(n>=1) then
      eps=0.01
      do i=1,iq
        up(i)=up(i)+eps*(uf(i)-2*up(i)+ub(i))
      end do
      do i=1,iq+1
        zp(i)=zp(i)+eps*(zf(i)-2*zp(i)+zb(i))
      end do
    end if
! データをファイル出力(noutステップごとに) ------------------
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
! データをファイル出力(noutステップごとに) ------------------
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
