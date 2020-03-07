!----------------------------------------------------------
!水深が変化する場合の津波の運動を差分近似式で解く
!----------------------------------------------------------
implicit none
!変数の型宣言 --------------------------------------
integer,parameter::iq=988
real*8:: uf(iq), up(iq), ub(iq)         !流速。f,p,bはそれぞれ未来，現在，過去に対応。
real*8:: zf(iq+1), zp(iq+1), zb(iq+1)   !海表面変位
real*8:: x(iq+1), depth(iq+1)
real*8:: grav, time_end, time_out, dx,&
         dt, eps, depth0, depth1
integer:: i, n, nend, nout, index
!
! パラメータの設定 -------------------------
grav=9.8              !重力加速度[m/s^2]
time_end=3.*60.*60.   !計算時間[s]
time_out=60.          !データ出力時間間隔[s]
! 水深データの読み込み ---------------------
open(10,file="DEPTH.data")
!
do i=1,iq+1
 read(10,*) x(i), depth(i)  !x[km]は岸からの距離。水深の符号はマイナス。
     depth(i)=-depth(i)
         x(i)=x(i)*1000.    !単位をキロからメートルへ
 enddo
! サイズ設定 --------------------------------
dx=(x(iq+1)-x(1))/iq  !空間グリッドサイズ
dt=0.5                ! 時間ステップサイズ
! 初期値の設定 ------------------------------      
do i=1,iq
up(i)=0.
enddo
 do i=1,iq+1
 zp(i)=-2.*exp(-dble(x(i)-115.*1000.)**2/dble(40.*1000.)**2)&   ! 引き波の初期値（x=115kmで-2mの沈降分） 
       +5.*exp(-dble(x(i)-215.*1000.)**2/dble(40.*1000.)**2)    ! 押し波の初期値（x=215kmで+5mの隆起分） 
! 津波の初期値は、国土地理院が推定した海底地形の隆起・沈降分布から定規を使って読み取って与える。
! http://www.gsi.go.jp/common/000060406.pdf 
enddo
!
! 出力データファイルを指定 ---------------------------------------------------
open(15,file='z.data')      ! 海表面の変位
open(25,file='u.data')      ! 流速
open(35,file='topo.data')   ! 地面または海底
!
! スッテプ指定 --------------------------------------------------------
nend=time_end/dt !時間ステップ数
nout=time_out/dt !出力データステップ間隔
index=0
!
! 時間ループ(n=0,nendまで繰り返す) ----------------------------------------
do n=0,nend
! 最初のステップのみ現在ステップ値から未来ステップ値を計算（前方差分
if(n==0) then                                        !n=0は初回の計算なので，過去の値がない
do i=1,iq
uf(i)=up(i)-grav*(dt/dx)*(zp(i+1)-zp(i))             !∂U/∂t=-g*∂Z/∂x
enddo
do i=2,iq
 depth0=(depth(i)+depth(i-1))*0.5
 depth1=(depth(i)+depth(i+1))*0.5
 zf(i)=zp(i)-(dt/dx)*(depth1*up(i)-depth0*up(i-1))   !∂Z/∂t=-∂(HU)/∂x
enddo
endif
! 過去ステップ値と現在ステップ値から未来ステップ値を計算（中央差分)
if(n>=1) then
do i=1,iq
uf(i)=ub(i)-2.*grav*(dt/dx)*(zp(i+1)-zp(i))           !∂U/∂t=-g*∂Z/∂x
enddo
do i=2,iq
 depth0=(depth(i)+depth(i-1))*0.5
 depth1=(depth(i)+depth(i+1))*0.5
 zf(i)=zb(i)-2.*(dt/dx)*(depth1*up(i)-depth0*up(i-1)) !∂Z/∂t=-∂(HU)/∂x
enddo
endif
! 陸上で(depth<0)で流速を0にする ------------------
do i=1,iq
 if((depth(i)+depth(i+1))*0.5<=0.) uf(i)=0.
enddo
!
! 沖側の津波を強制的に減衰させる ------------------
! iがiq-50以上のときは沖に行くほど波が小さくなるように
 do i=1,iq
  if(i>=iq-50) zf(i)=zf(i)*(iq-i)/50.
  if(i>=iq-50) uf(i)=uf(i)*(iq-i)/50.
 enddo
! 境界条件 -------
zf(1)=0.
zf(iq+1)=0.
!
! 計算の安定化のためのおまじない（Asselin filter）-------
if(n>=1) then
 eps=0.01
 do i=1,iq
  up(i)=up(i)+eps*(uf(i)-2.*up(i)+ub(i))
 enddo
 do i=1,iq+1
  zp(i)=zp(i)+eps*(zf(i)-2.*zp(i)+zb(i))
 enddo
 endif
!
! データをファイル出力(noutステップごとに) ------------------
if(mod(n,nout).eq.0) then
 do i=2,iq
  write(15,*) x(i)/1000., zp(i)
  write(25,*) x(i)/1000., up(i)  
  write(35,*) x(i)/1000., -depth(i)
 enddo
 write(15,*)
  write(15,*)
   write(25,*)
    write(25,*)
     write(35,*)
      write(35,*) 
       write(*,*) 'time(sec)=', dt*n, index
       index=index+1
       endif
!
!  ステップを進める
 do i=1,iq
  ub(i)=up(i)
 enddo
 do i=1,iq+1
  zb(i)=zp(i)
 enddo
 do i=1,iq
  up(i)=uf(i)
 enddo
 do i=1,iq+1
  zp(i)=zf(i)
 enddo
!
enddo
stop
end









