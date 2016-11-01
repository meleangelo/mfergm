
subroutine entropy(n,mu,ent)

!************************************
!* subroutine: entropy
!* version: 2.0
!* date: 05/15/2015
!* last modified: 
!* author: Angelo Mele
!* description: computes entropy of matrix mu.
!*              
!************************************
! entropy: function to compute entropy of network
!
! input variables
! n  = size of network
! mu = matrix of mu_ij 

! output variables
! ent = entropy
!

implicit none

!input variables
integer, intent(in):: n
real(8), intent(in), dimension(n,n)::mu

! output variables
real(8), intent(out):: ent 

! work variables
integer i,j,k  ! counters
real(8) logar

ent = 0.0
do i=1,n
    do j=i+1,n
        !ent = ent - mu(i,j)*log(mu(i,j) ) - (1-mu(i,j))*log(1-mu(i,j) )
        ent = ent - mu(i,j)*logar(mu(i,j) ) - (1-mu(i,j))*logar(1-mu(i,j) )
    enddo ! i
enddo ! j

end subroutine entropy




subroutine Eedges(n,mu,Emu)

!************************************
!* subroutine: Eedges
!* version: 2.0
!* date: 05/15/2015
!* last modified: 
!* author: Angelo Mele
!* description: computes expected number of edges of matrix mu.
!*              
!************************************
! Eedges: function to compute expected number of edges of matrix mu
!
! input variables
! n  = size of network
! mu = matrix of mu_ij 

! output variables
! Emu = entropy
!

implicit none

!input variables
integer, intent(in):: n
real(8), intent(in), dimension(n,n)::mu

! output variables
real(8), intent(out):: Emu 

! work variables
integer i,j,k  ! counters
Emu = 0.0
Emu = sum(mu)/2.0

end subroutine Eedges





subroutine E2star(n,mu,Emu2)

!************************************
!* subroutine: E2star
!* version: 2.0
!* date: 05/15/2015
!* last modified: 
!* author: Angelo Mele
!* description: computes expected number of 2-stars of matrix mu.
!*              
!************************************
! E2star: function to compute expected number of 2starss of matrix mu
!
! input variables
! n  = size of network
! mu = matrix of mu_ij 

! output variables
! Emu2 = exprected number of 2stars
!

implicit none

!input variables
integer, intent(in):: n
real(8), intent(in), dimension(n,n)::mu

! output variables
real(8), intent(out):: Emu2

! work variables
integer i,j,k  ! counters
Emu2 = 0.0
do i=1,n
    do j=i,n
        do k=1,n
            !if (k /= j .AND. k /= i) then
                Emu2 = Emu2 + mu(i,j)*mu(j,k) 
            !endif 
        enddo !k
    enddo !j
enddo !i
Emu2 = Emu2/n

end subroutine E2star



!************************************
!* function: dtfu
!* version: 2.0
!* date: 05/22/2009
!* last modified: 05/22/2009
!* author: Angelo Mele
!* description: contains the formula for the 
!*              change statistics of u_ij
!*              implemented with sparse network matrix
!************************************

function dtfu(nu,qu,x,xcolu,iu,ju,tu)

! dtfu: contains the formula for all the change statistics of u_ij
!
! INPUT
! nu  = number of agents in the network (INTEGER) 
! iu  = player i (INTEGER)
! ju  = player j (INTEGER)
! xcolu = column of x matrix to use
! tu = type of statistics (INTEGER)
!       1 = single link
!       2 = difference: xu(i)-xu(j)
!       3 = absolute difference: |xu(i)-xu(j)|
!       4 = sum: xu(i)+xu(j)
!       5 = indicator: I(x(i)=0,x(j)=0)
!       6 = indicator: I(x(i)=0,x(j)=1)
!       7 = indicator: I(x(i)=1,x(j)=0)
!       8 = indicator: I(x(i)=1,x(j)=1)
!       9 = indicator: I(x(i)=x(j))
!      10 = indicator: I(x(i) .ne. x(j))
!
! NOTE: if you want to add some statistics you can do it here
!       by adding additional terms to the IF conditions

implicit none

!**** input vbls ****
integer, intent(in)::nu,qu
real(8), intent(in), dimension(nu,qu)::x
integer, intent(in):: iu,ju,tu, xcolu
!**** output vbls ****
real(8) dtfu
!**** working vbls ****
!real(8) indic


dtfu=0.0
if (tu==1) then
	dtfu=1.0         !**** edges ****
else if (tu==2) then
	dtfu = x(iu, xcolu) - x(ju, xcolu)       !**** difference ****
else if (tu==3) then
	dtfu = abs(x(iu, xcolu) - x(ju, xcolu))         !**** absolute difference ****
else if (tu==4) then
	dtfu = x(iu, xcolu) + x(ju, xcolu)       !**** sum ****
else if (tu==5) then
	if (x(iu, xcolu)==0 .and. x(ju, xcolu)==0) dtfu = 1.0    !**** indicator both zeros ****
else if (tu==6) then
	if (x(iu, xcolu)==0 .and. x(ju, xcolu)==1) dtfu = 1.0      !**** indicator zero/one ****
else if (tu==7) then
	if (x(iu, xcolu)==1 .and. x(ju, xcolu)==0) dtfu = 1.0      !**** indicator one/zero ****
else if (tu==8) then
	if (x(iu, xcolu)==1 .and. x(ju, xcolu)==1) dtfu = 1.0      !**** indicator both ones ****
else if (tu==9) then
	if (x(iu, xcolu)==x(ju, xcolu))  dtfu = 1.0    !**** same ****
else if (tu==10) then 
	if (x(iu, xcolu) .ne. x(ju, xcolu))  dtfu = 1.0     !**** different ****
else if (tu==11) then
	dtfu = x(iu, xcolu)*x(ju, xcolu)        !**** product ****
else if (tu==12) then
	dtfu = x(iu, xcolu)             !**** x_i ****
else if (tu==13) then
	dtfu = x(ju, xcolu)            !**** x_j ****

! **** DIFFERENTIAL HOMOPHILY ***
else if (tu==21) then
	if (x(iu, xcolu)==1 .and. x(ju, xcolu)==1)  dtfu = 1.0
else if (tu==22) then
	if (x(iu, xcolu)==1 .and. x(ju, xcolu)==2)  dtfu = 1.0
else if (tu==23) then
	if (x(iu, xcolu)==1 .and. x(ju, xcolu)==3) dtfu = 1.0
else if (tu==24) then
	if (x(iu, xcolu)==1 .and. x(ju, xcolu)==4) dtfu = 1.0
else if (tu==25) then
	if (x(iu, xcolu)==1 .and. x(ju, xcolu)==5) dtfu = 1.0
else if (tu==31) then
	if (x(iu, xcolu)==2 .and. x(ju, xcolu)==1) dtfu = 1.0
else if (tu==32) then
	if (x(iu, xcolu)==2 .and. x(ju, xcolu)==2) dtfu = 1.0
else if (tu==33) then
	if (x(iu, xcolu)==2 .and. x(ju, xcolu)==3) dtfu = 1.0
else if (tu==34) then
	if (x(iu, xcolu)==2 .and. x(ju, xcolu)==4) dtfu = 1.0
else if (tu==35) then
	if (x(iu, xcolu)==2 .and. x(ju, xcolu)==5) dtfu = 1.0
else if (tu==41) then
	if (x(iu, xcolu)==3 .and. x(ju, xcolu)==1) dtfu = 1.0
else if (tu==42) then
	if (x(iu, xcolu)==3 .and. x(ju, xcolu)==2) dtfu = 1.0
else if (tu==43) then
	if (x(iu, xcolu)==3 .and. x(ju, xcolu)==3) dtfu = 1.0
else if (tu==44) then
	if (x(iu, xcolu)==3 .and. x(ju, xcolu)==4) dtfu = 1.0
else if (tu==45) then
	if (x(iu, xcolu)==3 .and. x(ju, xcolu)==5) dtfu = 1.0
else if (tu==51) then
	if (x(iu, xcolu)==4 .and. x(ju, xcolu)==1) dtfu = 1.0
else if (tu==52) then
	if (x(iu, xcolu)==4 .and. x(ju, xcolu)==2) dtfu = 1.0
else if (tu==53) then
	if (x(iu, xcolu)==4 .and. x(ju, xcolu)==3) dtfu = 1.0
else if (tu==54) then
	if (x(iu, xcolu)==4 .and. x(ju, xcolu)==4) dtfu = 1.0
else if (tu==55) then
	if (x(iu, xcolu)==4 .and. x(ju, xcolu)==5) dtfu = 1.0
else if (tu==61) then
	if (x(iu, xcolu)==5 .and. x(ju, xcolu)==1) dtfu = 1.0
else if (tu==62) then
	if (x(iu, xcolu)==5 .and. x(ju, xcolu)==2) dtfu = 1.0
else if (tu==63) then
	if (x(iu, xcolu)==5 .and. x(ju, xcolu)==3) dtfu = 1.0
else if (tu==64) then
	if (x(iu, xcolu)==5 .and. x(ju, xcolu)==4) dtfu = 1.0
else if (tu==65) then
	if (x(iu, xcolu)==5 .and. x(ju, xcolu)==5) dtfu = 1.0

!! **** DIFF HOMOPHILY interacted with SHARES ****
!else if (tu==101) then
!	if (x(iu, xcolu)==1 .and. x(ju, xcolu)==1)  dtfu = 1.0*share(1)
!else if (tu==102) then
!	if (x(iu, xcolu)==2 .and. x(ju, xcolu)==2) dtfu = 1.0*share(2)
!else if (tu==103) then
!	if (x(iu, xcolu)==3 .and. x(ju, xcolu)==3) dtfu = 1.0*share(3)
!else if (tu==104) then
!	if (x(iu, xcolu)==4 .and. x(ju, xcolu)==4) dtfu = 1.0*share(4)
!else if (tu==105) then
!	if (x(iu, xcolu)==5 .and. x(ju, xcolu)==5) dtfu = 1.0*share(5)
!else if (tu==111) then
!	if (x(iu, xcolu)==1 .and. x(ju, xcolu)==1) dtfu = 1.0*(share(1)**2)
!else if (tu==112) then
!	if (x(iu, xcolu)==2 .and. x(ju, xcolu)==2) dtfu = 1.0*(share(2)**2)
!else if (tu==113) then
!	if (x(iu, xcolu)==3 .and. x(ju, xcolu)==3) dtfu = 1.0*(share(3)**2)
!else if (tu==114) then
!	if (x(iu, xcolu)==4 .and. x(ju, xcolu)==4) dtfu = 1.0*(share(4)**2)
!else if (tu==115) then
!	if (x(iu, xcolu)==5 .and. x(ju, xcolu)==5) dtfu = 1.0*(share(5)**2)
!
!! **** SHARES 
!else if (tu==201) then
!	dtfu = share(1)
!else if (tu==202) then
!	dtfu = share(2)
!else if (tu==203) then
!	dtfu = share(3)
!else if (tu==204) then
!	dtfu = share(4)
!else if (tu==205) then
!	dtfu = share(5)
!
!! **** SHARES SQUARED
!else if (tu==211) then
!	dtfu = share(1)**2
!else if (tu==212) then
!	dtfu = share(2)**2
!else if (tu==213) then
!	dtfu = share(3)**2
!else if (tu==214) then
!	dtfu = share(4)**2
!else if (tu==215) then
!	dtfu = share(5)**2

endif

!dtfu = dtfu/(nu*(nu-1))
!dtfu = dtfu/nu
!dtfu = dtfu/10000.0
return
end  function dtfu



!************************************
!* function: dtfv
!* version: 2.0
!* date: 05/22/2009
!* last modified: 6/25/2014
!* author: Angelo Mele
!* description: contains the formula for the 
!*              change statistics for v_ij  
!*              
!************************************

function dtfv(nv,qv,muv,x,xcolv,iv,jv,tvq)

! dtfvq: contains the formula for all the change statistics
!
! INPUT
! nv   = number of agents in the network (INTEGER) 
! qv   = number of columns in x (integer)
! mu   = mean-field matrix (real(8), dimension(nv,nv))
! x    = matrix with exogenous characteristics (real, dimension(nv,qv))
! xcolv= column of x matrix to use
! iv   = player i (INTEGER)
! jv   = player j (INTEGER)
! tv   = type of statistics (INTEGER)
!       1 = single link
!       2 = difference: x(i)-x(j)
!       3 = absolute difference: |x(i)-x(j)|
!       4 = sum: x(i)+x(j)
!       5 = indicator: I(x(i)=0,x(j)=0)
!       6 = indicator: I(x(i)=0,x(j)=1)
!       7 = indicator: I(x(i)=1,x(j)=0)
!       8 = indicator: I(x(i)=1,x(j)=1)
!       9 = indicator: I(x(i)=x(j))
!      10 = indicator: I(x(i) .ne. x(j))
!
! NOTE: if you want to add some statistics you can do it here
!       by adding additional terms to the IF conditions


implicit none

! input vbls
integer, intent(in)::nv,qv
real(8), intent(in), dimension(nv,nv)::muv
real(8), intent(in), dimension(nv,qv)::x
integer, intent(in)::xcolv
integer, intent(in):: iv,jv,tvq
! output vbls
real(8) dtfv
! working vbls
integer i0,k0,ch
!real(8) indic

dtfv=0.0


!*********************
!*    single link    *
!*********************
if (tvq==1) then
    do k0=1,nv
        !if (k0 .ne. iv) then
            dtfv = dtfv + muv(jv,k0) + muv(k0,iv)
            dtfv = dtfv/nv
        !endif !(k0 .ne. iv)
    enddo !k0
	
!!********************
!!*    difference    *
!!********************
!else if (tvq==2) then  
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv) - x(k0, xcolv) )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!
!        if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv) - x(jv, xcolv) )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!
!    enddo !k0
!
!
!!***********************
!!* absolute difference *
!!***********************
!else if (tvq==3) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + abs( x(iv, xcolv) - x(k0, xcolv) )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!        if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + abs( x(k0, xcolv) - x(jv, xcolv) )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!
! 
!    enddo !k0
!
!!*************
!!*    sum    *
!!*************
!else if (tvq==4) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv) + x(k0, xcolv) )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv) + x(jv, xcolv) )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!
!         
!    enddo !k0
!	
!!*****************************
!!*    Indicator zero/zero    *
!!*****************************
!else if (tvq==5) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                if ( x(iv, xcolv)==0 .and. x(k0, xcolv) ==0 ) then
!                    dtfv = dtfv + 1.0
!                endif !( x(iv, xcolv)==0 .and. x(k0, xcolv) ==0 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                if ( x(k0, xcolv)==0 .and. x(jv, xcolv) ==0 ) then
!                    dtfv = dtfv + 1.0
!                endif !( x(k0, xcolv)==0 .and. x(jv, xcolv) ==0 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!
!         
!    enddo !k0
!	
!!*****************************
!!*    Indicator zero/one     *
!!*****************************
!else if (tvq==6) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==0 )*( x(k0, xcolv) ==1 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==0)*(x(jv, xcolv) ==1 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!	
!!*****************************
!!*    Indicator one/zero     *
!!*****************************
!else if (tvq==7) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==1 )*( x(k0, xcolv) ==0 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==1)*(x(jv, xcolv) ==0 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!	
!!*****************************
!!*    Indicator one/one      *
!!*****************************
!else if (tvq==8) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                if ( x(iv, xcolv)==1 .and. x(k0, xcolv) ==1 ) then
!                    dtfv = dtfv + 1.0
!                endif !( x(iv, xcolv)==1 .and. x(k0, xcolv) ==1 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                if ( x(k0, xcolv)==1 .and. x(jv, xcolv) ==1 ) then
!                    dtfv = dtfv + 1.0
!                endif !( x(k0, xcolv)==1)*(x(jv, xcolv) ==1 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!*****************************
!!*    Indicator x(i)=x(j)    *
!!*****************************
!else if (tvq==9) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                if ( x(iv, xcolv) == x(k0, xcolv) ) then
!                    dtfv = dtfv + 1.0
!                endif !( x(iv, xcolv) == x(k0, xcolv) )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                if ( x(k0, xcolv)==x(jv, xcolv) ) then
!                    dtfv = dtfv + 1.0
!                endif !( x(k0, xcolv)==x(jv, xcolv) )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!*****************************
!!*    Indicator x(i)/=x(j)   *
!!*****************************
!else if (tvq==10) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                if ( x(iv, xcolv) .ne. x(k0, xcolv) ) then
!                    dtfv = dtfv + 1.0
!                endif ! ( x(iv, xcolv) .ne. x(k0, xcolv) )               
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                if ( x(k0, xcolv) .ne. x(jv, xcolv) ) then
!                    dtfv = dtfv + 1.0
!                endif !( x(k0, xcolv) .ne. x(jv, xcolv) )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!*****************************
!!*    Product x(i)*x(k)      *
!!*****************************
!else if (tvq==11) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)*x(k0, xcolv) )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)*x(jv, xcolv) )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!*****************************
!!*     x(j)      *
!!*****************************
!else if (tvq==12) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + x(jv, xcolv)
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!        if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + x(k0, xcolv)
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!   enddo !k0
!
!!******************************
!!* Differential Homophily 1/1 *
!!******************************
!else if (tvq==21) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==1 )*( x(k0, xcolv) ==1)
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==1)*(x(jv, xcolv) ==1 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 1/2 *
!!******************************
!else if (tvq==22) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==1 )*( x(k0, xcolv) ==2 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==1)*(x(jv, xcolv) ==2 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 1/3 *
!!******************************
!else if (tvq==23) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==1 )*( x(k0, xcolv) ==3 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==1)*(x(jv, xcolv) ==3 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 1/4 *
!!******************************
!else if (tvq==24) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==1 )*( x(k0, xcolv) ==4 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==1)*(x(jv, xcolv) ==4 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 1/5 *
!!******************************
!else if (tvq==25) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==1 )*( x(k0, xcolv) ==5 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==1)*(x(jv, xcolv) ==5 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 2/1 *
!!******************************
!else if (tvq==31) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==2 )*( x(k0, xcolv) ==1 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==2)*(x(jv, xcolv) ==1 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 2/2 *
!!******************************
!else if (tvq==32) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==2 )*( x(k0, xcolv) ==2 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==2)*(x(jv, xcolv) ==2 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 2/3 *
!!******************************
!else if (tvq==33) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==2 )*( x(k0, xcolv) ==3 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==2)*(x(jv, xcolv) ==3 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 2/4 *
!!******************************
!else if (tvq==34) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==2 )*( x(k0, xcolv) ==4 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==2)*(x(jv, xcolv) ==4 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 2/5 *
!!******************************
!else if (tvq==35) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==2 )*( x(k0, xcolv) ==5 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==2)*(x(jv, xcolv) ==5 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 3/1 *
!!******************************
!else if (tvq==41) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv= dtfv + ( x(iv, xcolv)==3 )*( x(k0, xcolv) ==1 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==3)*(x(jv, xcolv) ==1 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 3/2 *
!!******************************
!else if (tvq==42) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==3 )*( x(k0, xcolv) ==2 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==3)*(x(jv, xcolv) ==2 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 3/3 *
!!******************************
!else if (tvq==43) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==3 )*( x(k0, xcolv) ==3 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==3)*(x(jv, xcolv) ==3 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 3/4 *
!!******************************
!else if (tvq==44) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==3 )*( x(k0, xcolv) ==4 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==3)*(x(jv, xcolv) ==4 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 3/5 *
!!******************************
!else if (tvq==45) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==3 )*( x(k0, xcolv) ==5 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==3)*(x(jv, xcolv) ==5 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 4/1 *
!!******************************
!else if (tvq==51) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==4 )*( x(k0, xcolv) ==1 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==4)*(x(jv, xcolv) ==1 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 4/2 *
!!******************************
!else if (tvq==52) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==4 )*( x(k0, xcolv) ==2 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==4)*(x(jv, xcolv) ==2 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 4/3 *
!!******************************
!else if (tvq==53) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==4 )*( x(k0, xcolv) ==3 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==4)*(x(jv, xcolv) ==3 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 4/4 *
!!******************************
!else if (tvq==54) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==4 )*( x(k0, xcolv) ==4 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==4)*(x(jv, xcolv) ==4 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 4/5 *
!!******************************
!else if (tvq==55) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==4 )*( x(k0, xcolv) ==5 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==4)*(x(jv, xcolv) ==5 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 5/1 *
!!******************************
!else if (tvq==61) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==5 )*( x(k0, xcolv) ==1 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==5)*(x(jv, xcolv) ==1 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 5/2 *
!!******************************
!else if (tvq==62) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==5 )*( x(k0, xcolv) ==2 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==5)*(x(jv, xcolv) ==2 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 5/3 *
!!******************************
!else if (tvq==63) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==5 )*( x(k0, xcolv) ==3 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==5)*(x(jv, xcolv) ==3 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 5/4 *
!!******************************
!else if (tvq==64) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==5 )*( x(k0, xcolv) ==4 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==5)*(x(jv, xcolv) ==4 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 5/5 *
!!******************************
!else if (tvq==65) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==5 )*( x(k0, xcolv) ==5 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==5)*(x(jv, xcolv) ==5 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!else
!	dtfv=0.0
!endif
!
!
!!*********************
!!*    cyclic triangles   *
!!*********************
!elseif (tvq==100) then
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            dtfv = dtfv + 1.0*mu(jv,k0)*mu(k0,iv)
!        endif !(k0 .ne. iv)
! 
!    enddo !k0
!	
else
	dtfv=0.0
endif

!dtfv = dtfv/(nv*(nv-1)*(nv-2))
!dtfv = dtfv/(nv*(nv-1))
!dtfv = dtfv/10000.0

return
end function dtfv


!************************************
!* function: dtfvq
!* version: 2.0
!* date: 05/22/2009
!* last modified: 6/25/2014
!* author: Angelo Mele
!* description: contains the formula for the 
!*              change statistics for v_ij  
!*              
!************************************

function dtfvq(nv,qv,muv,x,xcolv,iv,jv,tvq)

! dtfvq: contains the formula for all the change statistics
!
! INPUT
! nv   = number of agents in the network (INTEGER) 
! qv   = number of columns in x (integer)
! mu   = mean-field matrix (real(8), dimension(nv,nv))
! x    = matrix with exogenous characteristics (real, dimension(nv,qv))
! xcolv= column of x matrix to use
! iv   = player i (INTEGER)
! jv   = player j (INTEGER)
! tv   = type of statistics (INTEGER)
!       1 = single link
!       2 = difference: x(i)-x(j)
!       3 = absolute difference: |x(i)-x(j)|
!       4 = sum: x(i)+x(j)
!       5 = indicator: I(x(i)=0,x(j)=0)
!       6 = indicator: I(x(i)=0,x(j)=1)
!       7 = indicator: I(x(i)=1,x(j)=0)
!       8 = indicator: I(x(i)=1,x(j)=1)
!       9 = indicator: I(x(i)=x(j))
!      10 = indicator: I(x(i) .ne. x(j))
!
! NOTE: if you want to add some statistics you can do it here
!       by adding additional terms to the IF conditions


implicit none

! input vbls
integer, intent(in)::nv,qv
real(8), intent(in), dimension(nv,nv)::muv
real(8), intent(in), dimension(nv,qv)::x
integer, intent(in)::xcolv
integer, intent(in):: iv,jv,tvq
! output vbls
real(8) dtfvq
! working vbls
integer i0,k0,ch
!real(8) indic

dtfvq=0.0


!*********************
!*    single link    *
!*********************
if (tvq==1) then
    do k0=1,nv
        !if (k0 .ne. iv) then
            dtfvq = dtfvq + muv(jv,k0)/nv
        !endif !(k0 .ne. iv)
    enddo !k0
	
!!********************
!!*    difference    *
!!********************
!else if (tvq==2) then  
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv) - x(k0, xcolv) )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!
!        if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv) - x(jv, xcolv) )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!
!    enddo !k0
!
!
!!***********************
!!* absolute difference *
!!***********************
!else if (tvq==3) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + abs( x(iv, xcolv) - x(k0, xcolv) )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!        if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + abs( x(k0, xcolv) - x(jv, xcolv) )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!
! 
!    enddo !k0
!
!!*************
!!*    sum    *
!!*************
!else if (tvq==4) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv) + x(k0, xcolv) )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv) + x(jv, xcolv) )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!
!         
!    enddo !k0
!	
!!*****************************
!!*    Indicator zero/zero    *
!!*****************************
!else if (tvq==5) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                if ( x(iv, xcolv)==0 .and. x(k0, xcolv) ==0 ) then
!                    dtfv = dtfv + 1.0
!                endif !( x(iv, xcolv)==0 .and. x(k0, xcolv) ==0 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                if ( x(k0, xcolv)==0 .and. x(jv, xcolv) ==0 ) then
!                    dtfv = dtfv + 1.0
!                endif !( x(k0, xcolv)==0 .and. x(jv, xcolv) ==0 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!
!         
!    enddo !k0
!	
!!*****************************
!!*    Indicator zero/one     *
!!*****************************
!else if (tvq==6) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==0 )*( x(k0, xcolv) ==1 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==0)*(x(jv, xcolv) ==1 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!	
!!*****************************
!!*    Indicator one/zero     *
!!*****************************
!else if (tvq==7) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==1 )*( x(k0, xcolv) ==0 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==1)*(x(jv, xcolv) ==0 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!	
!!*****************************
!!*    Indicator one/one      *
!!*****************************
!else if (tvq==8) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                if ( x(iv, xcolv)==1 .and. x(k0, xcolv) ==1 ) then
!                    dtfv = dtfv + 1.0
!                endif !( x(iv, xcolv)==1 .and. x(k0, xcolv) ==1 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                if ( x(k0, xcolv)==1 .and. x(jv, xcolv) ==1 ) then
!                    dtfv = dtfv + 1.0
!                endif !( x(k0, xcolv)==1)*(x(jv, xcolv) ==1 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!*****************************
!!*    Indicator x(i)=x(j)    *
!!*****************************
!else if (tvq==9) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                if ( x(iv, xcolv) == x(k0, xcolv) ) then
!                    dtfv = dtfv + 1.0
!                endif !( x(iv, xcolv) == x(k0, xcolv) )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                if ( x(k0, xcolv)==x(jv, xcolv) ) then
!                    dtfv = dtfv + 1.0
!                endif !( x(k0, xcolv)==x(jv, xcolv) )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!*****************************
!!*    Indicator x(i)/=x(j)   *
!!*****************************
!else if (tvq==10) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                if ( x(iv, xcolv) .ne. x(k0, xcolv) ) then
!                    dtfv = dtfv + 1.0
!                endif ! ( x(iv, xcolv) .ne. x(k0, xcolv) )               
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                if ( x(k0, xcolv) .ne. x(jv, xcolv) ) then
!                    dtfv = dtfv + 1.0
!                endif !( x(k0, xcolv) .ne. x(jv, xcolv) )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!*****************************
!!*    Product x(i)*x(k)      *
!!*****************************
!else if (tvq==11) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)*x(k0, xcolv) )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)*x(jv, xcolv) )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!*****************************
!!*     x(j)      *
!!*****************************
!else if (tvq==12) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + x(jv, xcolv)
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!        if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + x(k0, xcolv)
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!   enddo !k0
!
!!******************************
!!* Differential Homophily 1/1 *
!!******************************
!else if (tvq==21) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==1 )*( x(k0, xcolv) ==1)
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==1)*(x(jv, xcolv) ==1 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 1/2 *
!!******************************
!else if (tvq==22) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==1 )*( x(k0, xcolv) ==2 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==1)*(x(jv, xcolv) ==2 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 1/3 *
!!******************************
!else if (tvq==23) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==1 )*( x(k0, xcolv) ==3 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==1)*(x(jv, xcolv) ==3 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 1/4 *
!!******************************
!else if (tvq==24) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==1 )*( x(k0, xcolv) ==4 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==1)*(x(jv, xcolv) ==4 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 1/5 *
!!******************************
!else if (tvq==25) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==1 )*( x(k0, xcolv) ==5 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==1)*(x(jv, xcolv) ==5 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 2/1 *
!!******************************
!else if (tvq==31) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==2 )*( x(k0, xcolv) ==1 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==2)*(x(jv, xcolv) ==1 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 2/2 *
!!******************************
!else if (tvq==32) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==2 )*( x(k0, xcolv) ==2 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==2)*(x(jv, xcolv) ==2 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 2/3 *
!!******************************
!else if (tvq==33) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==2 )*( x(k0, xcolv) ==3 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==2)*(x(jv, xcolv) ==3 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 2/4 *
!!******************************
!else if (tvq==34) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==2 )*( x(k0, xcolv) ==4 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==2)*(x(jv, xcolv) ==4 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 2/5 *
!!******************************
!else if (tvq==35) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==2 )*( x(k0, xcolv) ==5 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==2)*(x(jv, xcolv) ==5 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 3/1 *
!!******************************
!else if (tvq==41) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv= dtfv + ( x(iv, xcolv)==3 )*( x(k0, xcolv) ==1 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==3)*(x(jv, xcolv) ==1 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 3/2 *
!!******************************
!else if (tvq==42) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==3 )*( x(k0, xcolv) ==2 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==3)*(x(jv, xcolv) ==2 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 3/3 *
!!******************************
!else if (tvq==43) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==3 )*( x(k0, xcolv) ==3 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==3)*(x(jv, xcolv) ==3 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 3/4 *
!!******************************
!else if (tvq==44) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==3 )*( x(k0, xcolv) ==4 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==3)*(x(jv, xcolv) ==4 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 3/5 *
!!******************************
!else if (tvq==45) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==3 )*( x(k0, xcolv) ==5 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==3)*(x(jv, xcolv) ==5 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 4/1 *
!!******************************
!else if (tvq==51) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==4 )*( x(k0, xcolv) ==1 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==4)*(x(jv, xcolv) ==1 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 4/2 *
!!******************************
!else if (tvq==52) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==4 )*( x(k0, xcolv) ==2 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==4)*(x(jv, xcolv) ==2 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 4/3 *
!!******************************
!else if (tvq==53) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==4 )*( x(k0, xcolv) ==3 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==4)*(x(jv, xcolv) ==3 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 4/4 *
!!******************************
!else if (tvq==54) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==4 )*( x(k0, xcolv) ==4 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==4)*(x(jv, xcolv) ==4 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 4/5 *
!!******************************
!else if (tvq==55) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==4 )*( x(k0, xcolv) ==5 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==4)*(x(jv, xcolv) ==5 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 5/1 *
!!******************************
!else if (tvq==61) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==5 )*( x(k0, xcolv) ==1 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==5)*(x(jv, xcolv) ==1 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 5/2 *
!!******************************
!else if (tvq==62) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==5 )*( x(k0, xcolv) ==2 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==5)*(x(jv, xcolv) ==2 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 5/3 *
!!******************************
!else if (tvq==63) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==5 )*( x(k0, xcolv) ==3 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==5)*(x(jv, xcolv) ==3 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 5/4 *
!!******************************
!else if (tvq==64) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==5 )*( x(k0, xcolv) ==4 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==5)*(x(jv, xcolv) ==4 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 5/5 *
!!******************************
!else if (tvq==65) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==5 )*( x(k0, xcolv) ==5 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==5)*(x(jv, xcolv) ==5 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!else
!	dtfv=0.0
!endif
!
!
!!*********************
!!*    cyclic triangles   *
!!*********************
!elseif (tvq==100) then
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            dtfv = dtfv + 1.0*mu(jv,k0)*mu(k0,iv)
!        endif !(k0 .ne. iv)
! 
!    enddo !k0
!	
else
	dtfvq=0.0
endif

!dtfv = dtfv/(nv*(nv-1)*(nv-2))
!dtfv = dtfv/(nv*(nv-1))
!dtfv = dtfv/10000.0

return
end function dtfvq



!************************************
!* subroutine: psi_mf
!* version: 1.0
!* date: 06/22/2014
!* last modified: 
!* author: Angelo Mele
!* description: computes log-constant
!*              approdximation with mean field
!*              
!************************************
!subroutine psi_mf(n,q,pp,mu,tmu,theta,x,dt, tol, psi)
subroutine psi_mf(n,q,pp,mu,theta,x,dt, tol, psi)

implicit none

!input variables
integer, intent(in)::n
integer, intent(in)::q
integer, intent(in)::pp(2)
real(8), intent(inout), dimension(n,n)::mu
!real(8), intent(in), dimension(pp(2))::tmu
real(8), intent(in), dimension(pp(2))::theta
real(8), intent(in), dimension(n,q)::x
integer, intent(in), dimension(pp(2),2)::dt
real(8), intent(in)::tol

!output variables
real(8), intent(out):: psi

!working variables
integer i,j,k,kk  !counters
real(8) eps
real(8), dimension(pp(2))::dttemp
real(8) muij_new
real(8) dtfu,dtfv
!real(8), dimension(n,n)::mu
real(8), dimension(pp(2))::tmu
!real(8), dimension(pp(2))::scaletheta
real(8) potential 
real(8) In
real(8) logar

! initialize
!mu = mux
psi = 0.0
eps = 10.0
i = 1
j = 1
dttemp(:) = 0.0

! initialize  psi
tmu(:)=0.0
call sstat(pp,n,q,mu,x,dt,tmu)
!scaletheta(:) = 1.0
!scaletheta(pp(1)+1:pp(2)) = 0.5
potential = dot_product(theta,tmu)
call entropy(n,mu,In)
psi = potential + In

! loop over matrix mu to update
do while (eps>tol) 
            if(i/=j) then
                ! direct utility update
                do kk=1,pp(1)
                    dttemp(kk) = 2.0*dtfu(n,q,x,dt(kk,1),i,j,dt(kk,2)) 
                enddo
                ! indirect utility
                do kk=pp(1)+1,pp(2)
                    dttemp(kk) = dtfv(n,q,mu,x,dt(kk,1),i,j,dt(kk,2)) 
                enddo
   !print*, 'i=',i,'j=',j
                ! update mu(i,j)
                muij_new = 1.0/exp(-dot_product(theta,dttemp))
                !print*, 'mij_new', muij_new
                eps = dot_product(theta,dttemp)*(muij_new-mu(i,j))
                !print*, eps
!                eps = eps - muij_new*log(muij_new) - (1-muij_new)*log(1-muij_new)
                eps = eps - muij_new*logar(muij_new) - (1-muij_new)*logar(1-muij_new)
                !print*, eps
 !               eps = eps + mu(i,j)*log(mu(i,j)) + (1-mu(i,j))*log(1-mu(i,j))
                eps = eps + mu(i,j)*logar(mu(i,j)) + (1-mu(i,j))*logar(1-mu(i,j))
                !print*, 'eps', eps
    
                ! update log-constant estimate        
                psi = psi + eps   
                !print*, 'psi', psi
                mu(i,j) = muij_new
                mu(j,i) = mu(i,j)
            endif !(i/=j)
            
            if (i==n .and. j==n) then 
                i = 1
                j = 1
            elseif (i<n .and. j==n) then
                j = 1
                i = i+1
            elseif (i<n .and. j<n) then
                j = j+1
            endif
enddo !while loop


psi = psi/(n*n)

end subroutine psi_mf





!************************************
!* subroutine: psi_mf2
!* version: 1.0
!* date: 07/06/2016
!* last modified: 07/13/2016
!* author: Angelo Mele
!* description: computes log-constant
!*              approdximation with mean field
!*              improved version      
!************************************
!subroutine psi_mf(n,q,pp,mu,tmu,theta,x,dt, tol, psi)
subroutine psi_mf2(n,q,pp,mu,theta,x,dt, tol, maxiterations, psi)

implicit none

!input variables
integer, intent(in)::n
integer, intent(in)::q
integer, intent(in)::pp(2)
real(8), intent(inout), dimension(n,n)::mu
!real(8), intent(in), dimension(pp(2))::tmu
real(8), intent(in), dimension(pp(2))::theta
real(8), intent(in), dimension(n,q)::x
integer, intent(in), dimension(pp(2),2)::dt
real(8), intent(in)::tol
integer, intent(in)::maxiterations

!output variables
real(8), intent(out):: psi

!working variables
real(8) psinew
integer t,i,j,k,kk  !counters
real(8) eps
real(8), dimension(pp(2))::dttemp
!real(8) muij_new
real(8) dtfu,dtfv
!real(8), dimension(n,n)::mu
real(8), dimension(pp(2))::tmu
!real(8), dimension(pp(2))::scaletheta
real(8) potential 
real(8) In
real(8) logar

! initialize
!mu = mux
psi = 0.0
eps = 10.0
dttemp(:) = 0.0

! initialize  psi
tmu(:)=0.0
call sstat(pp,n,q,mu,x,dt,tmu)
!scaletheta(:) = 1.0
!scaletheta(pp(1)+1:pp(2)) = 0.5
potential = dot_product(theta,tmu)
call entropy(n,mu,In)
psi = potential + In


! loop over matrix mu to update
do t = 2,maxiterations
    do i=1,n
        do j=i,n
            if(i/=j) then
                ! direct utility update
                do kk=1,pp(1)
                    dttemp(kk) = 2.0*dtfu(n,q,x,dt(kk,1),i,j,dt(kk,2)) 
                enddo
                ! indirect utility
                do kk=pp(1)+1,pp(2)
                    dttemp(kk) = dtfv(n,q,mu,x,dt(kk,1),i,j,dt(kk,2)) 
                enddo
                ! update mu(i,j)
                mu(i,j) = 1.0/exp(-dot_product(theta,dttemp))
                mu(j,i) = mu(i,j)
            endif !(i/=j)
        enddo ! j
    enddo ! i 
    
    ! compute psi again
    call sstat(pp,n,q,mu,x,dt,tmu)
    potential = dot_product(theta,tmu)
    call entropy(n,mu,In)
    psinew = potential + In
    
    ! check convergence
    eps = psinew - psi
    
    ! if converged, stop
    if (eps < tol) then
        exit 
    endif
    ! if did not converge, update psi and update mu again
    psi = psinew
enddo ! t

! divide by n^2
psi = psi/(n*n)

end subroutine psi_mf2




!************************************
!* subroutine: sstat
!* version: 2.0
!* date: 06/22/2014
!* last modified: 
!* author: Angelo MeleF
!* description: computes the sufficient stats
!*              given a vector of change statistics
!*              
!************************************
subroutine sstat(ps,ns,qs,mu,x,dts,ts)

! SSTAT: functions that computes ALL sufficient statistis
!        given a vector of change statistics
!
! INPUT
! ps   = vector with number of suff stats (ppu,ppm,ppv) (INTEGER(3))
! ns   = number of players (INTEGER)
! qs   = number of variables in xxx (columns of xxx), INTEGER
! mu    = matrix mu with the network (real, dimension(ns,ns))
! x    = matrix with exogenous variables (real, dimesion(ns,qs))
! dts  = vector containing the specification (INTEGER, DIMENSION(ppp,2)
!
! OUTPUT
! ts   = vector of sufficient statistics


implicit none

! inpur vbls
integer, dimension(2), intent(in):: ps
integer, intent(in):: ns
integer, intent(in):: qs
real(8), intent(in):: mu(ns,ns)
real(8), dimension(ns,qs), intent(in):: x
integer, dimension(ps(2),2), intent(in):: dts

! output vbls
real(8), dimension(ps(2)), intent(out):: ts
! working vbls
integer i,j,j1,k
real(8) dtfu,dtfvq


ts(1:ps(2))=0.0

 
! statistics for u
do i=1,ns
	do j = i,ns
        do k=1,ps(1)
    		ts(k) = ts(k) + 2.0*mu(i,j)*dtfu(ns,qs,x,dts(k,1),i,j,dts(k,2))
        enddo !k
	enddo !j
enddo !i

! statistics for v
do i=1,ns
	do j = i,ns
        
        do k=ps(1)+1,ps(2)
			ts(k) = ts(k) + mu(i,j)*dtfvq(ns,qs,mu,x,dts(k,1),i,j,dts(k,2))
        enddo !k
	enddo !j
enddo !i

!ts(:) = 10.0
end subroutine sstat




!************************************
!* function: logar
!* version: 1.0
!* date: 09/28/2015
!* last modified: 09/28/2015
!* author: Angelo Mele
!* description: computes logs, preventing underflow
!*              
!************************************
function logar(number)

! LOGAR : functions that computes log of a number
!         avoiding underflow for small numbers
!
! INPUT
! number = real number 

implicit none

real(8) number
real(8) logar

if (number < 1.0E-30) then
    logar = log(1.0E-30)
else
    logar = log(number)
endif
end function


