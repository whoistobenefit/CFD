integer,parameter:: mfi=101,mr=151,ma=max(mfi,mr),mi=1
real, dimension(1:mfi,1:mr)::tet,tets,tet1,tet2,vi,vis,vi1,vi2,ft,ft1,ft2,fts,ufi,ur
real, dimension(mi:ma)::a,b,c,d,e
character(100):: sim ,name

!namelist/CHANNEL/ Re,Pr, dt, v_trub, dfi, dr, sumft, sumvi, sumtet

!----------------
itk = 400
ish = 100
it = 0

dt=0.001
Re=10.
Pr=0.1
Gr=0.

!channel geometry
alf=30.
r1=1.
r2=2.

sumft = 0.
sumvi = 0.
sumtet = 0.

dr=1./(mr-1)
dr2=dr*dr
             
dfi=alf*r2*2.*3.1415/(mfi-1)/360.
dfi2=dfi*dfi

!temperature I.C.    
tet(1,2:mr-1)=0.
tet(mfi,2:mr-1)=1.
tet(2:mfi-1,1)=0.
tets=tet

!I.C. for velocities
ur(2:mfi-1,1)=1.
ufi(2:mfi-1,1)=0.

!stream-function I.C.
do i=1,mfi
    fi=(i-1)*dfi
    ft(i,1)=-fi*r2          
enddo
ft(mfi,2:mr)=ft(mfi,1)
ft(1,1:mr)=0.
fts=ft


!MAIN CYCLE
do while (it<=itk)

    !--------------STREAM FUNCTION
    !sweep for y1
    do i=2, mfi-1
        do j=2,mr-1
            r=r2-(j-1)*dr
            a(j)=1./dr2
            c(j)=1./dr2
            b(j)=1./dt+2./dr2
            d(j)=fts(i,j)/dt-vi(i,j)+(fts(i,j+1)-fts(i,j-1))/2./dr/r
        enddo      
        a(1)=0.
        b(1)=1.
        c(1)=0.
        d(1)=fts(i,1)
        a(mr)=1.
        b(mr)=1.
        c(mr)=0.
        d(mr)=0.
        call Tom(1,mr,a,b,c,d,e,mi,ma)
        ft1(i,1:mr)=e(1:mr)
    enddo    
    
    !sweep for x1       
    do j=2,mr-1
        do i=2,mfi-1  
            r=r2-(j-1)*dr
            a(i)=1./dfi2/r**2
            c(i)=1./dfi2/r**2
            b(i)=1./dt+2./dfi2/r**2
            d(i)=ft1(i,j)/dt
        enddo 
        a(1)=0.
        b(1)=1.
        c(1)=0.
        d(1)=0.
        a(mfi)=0.
        b(mfi)=1.
        c(mfi)=0.
        d(mfi)=fts(mfi,j)
        call Tom(1,mfi,a,b,c,d,e,mi,ma)
        ft(1:mfi,j)=e(1:mfi)
    enddo 
    
    !local B.C. for stream function
    ft(1:mfi,mr)=ft(1:mfi,mr-1)
    ft2=ft-fts   
            
    !B.C. for velocities
    do i=2, mfi-1
        do j=2, mr-1
        r=r2-(j-1)*dr
        ufi(i,j)=(ft(i,j+1)-ft(i,j-1))/2./dr
        ur(i,j)=-(ft(i+1,j)-ft(i-1,j))/r/2./dfi
        enddo
    enddo
    do i=2,mfi-1
       ufi(i,mr)=ufi(i,mr-1)
       ur(i,mr)=ur(i,mr-1)*(1+dr/r1)
    enddo
    
    
    !--------------VORTICITY
    !sweep for y1
    do i=2, mfi-1
        do j=2,mr-1
            r=r2-(j-1)*dr
            aur=abs(ur(i,j))
            a(j)=(aur+ur(i,j))/dr/2.+1./Re/dr2
            c(j)=(aur-ur(i,j))/dr/2.+1./Re/dr2
            b(j)=1./dt+2./Re/dr2+aur/dr            
            d(j)=vis(i,j)/dt+(vis(i,j+1)-vis(i,j-1))/2./dr/r/Re
        enddo 
        a(mr)=1.
        b(mr)=1.
        c(mr)=0.
        d(mr)=0.
        a(1)=0.
        b(1)=1.
        c(1)=0.
        d(1)=0.
        call Tom(1,mr,a,b,c,d,e,mi,ma)
        vi1(i,1:mr)=e(1:mr)
    enddo
    
    !sweep for x1
    do j=2, mr-1
        do i=2,mfi-1
            r=r2-(j-1)*dr
            aufi=abs(ufi(i,j))
            a(i)=(aufi+ufi(i,j))/dfi/r/2.+1./dfi2/r**2/Re
            c(i)=(aufi-ufi(i,j))/dfi/r/2.+1./dfi2/r**2/Re
            b(i)=2./dfi2/Re/r**2+1./dt+aufi/dfi/r
            d(i)=vi1(i,j)/dt+(tet(i+1,j)-tet(i-1,j))/2./dfi*Gr/Re**2
        enddo 
        a(1)=0
        b(1)=1.
        c(1)=0.
        d(1)=2.*(ft(2,j)-ft(1,j))/dfi2/r**2
        a(mfi)=0.
        b(mfi)=1.
        c(mfi)=0.
        d(mfi)=2.*(ft(mfi-1,j)-ft(mfi,j))/dfi2/r**2
        call Tom(1,mfi,a,b,c,d,e,mi,ma)
        vi(1:mfi,j)=e(1:mfi)
    enddo
        
    !local B.C. for vorticity
    vi(2:mfi-1,mr)=vi(2:mfi-1,mr-1)
    vi2=vi-vis 


    !--------------TEMPERATURE
    !sweep for y1
    do ll=1,3
        do i=2, mfi-1
            do j=2,mr-1
                r=r2-(j-1)*dr
                aur=abs(ur(i,j))
                a(j)=(aur+ur(i,j))/dr/2.+1/dr2/Re/Pr
                c(j)=(aur-ur(i,j))/dr/2.+1/dr2/Re/Pr
                b(j)=1./dt+2./Re/Pr/dr2+aur/dr 
                d(j)=tets(i,j)/dt+(tets(i,j+1)-tets(i,j-1))/2./dr/r/Re/Pr
            enddo 
            a(mr)=1.
            b(mr)=1.
            c(mr)=0.
            d(mr)=0. 
            a(1)=0.
            b(1)=1.
            c(1)=0.
            d(1)=0. 
            call Tom(1,mr,a,b,c,d,e,mi,ma)
            tet1(i,1:mr)=e(1:mr)
        enddo
        
        !sweep for x1
        do j=2, mr-1
            do i=2,mfi-1
                r=r2-(j-1)*dr
                aufi=abs(ufi(i,j))
                a(i)=(aufi+ufi(i,j))/dfi/r/2.+1./dfi2/r**2/Re/Pr
                c(i)=(aufi-ufi(i,j))/dfi/r/2.+1./dfi2/r**2/Re/Pr
                b(i)=2./dfi2/r**2/Re/Pr+1./dt+aufi/dfi/r
                d(i)=tet1(i,j)/dt
            enddo
            a(1)=0.
            b(1)=1.
            c(1)=1.
            d(1)=0.
            a(mfi)=0.
            b(mfi)=1.
            c(mfi)=0.
            d(mfi)=1.!tet(mfi,j) 
            call Tom(1,mfi,a,b,c,d,e,mi,ma)
            tet(1:mfi,j)=e(1:mfi)
        enddo 
       
       !local B.C. for temperature
        tet(2:mfi-1,mr)=tet(2:mfi-1,mr-1)
        tet(1,2:mr-1)=tet(2,2:mr-1)
        tet2=tet-tets
        tets=tet
    enddo

    vis=vi
    fts=ft
    tets=tet
 
    !CONVERGENCE CHECK
    do i=2, mfi-1
        do j=2, mr-1
            sumft= sumft+abs(ft2(i,j))/dt
            sumvi= sumvi+abs(vi2(i,j))/dt
            sumtet= sumtet+abs(tet2(i,j))/dt
        enddo
    enddo
    
    if (it==itk-1) then
        write(*, '(I7,3es14.4)') it, sumft, sumvi, sumtet
    endif
    
    it = it+1
enddo !while


!--------------------- Write to files--------------------
! sim=Adjustl(sim)
open(23,file=trim('name1')//' izo_tem_ft_ufi_ur_vi '//' .dat')
do i=1,mfi
    do j=1,mr-1
        fi=(i-1)*dfi
        r=r2-(j-1)*dr
		write(23,'(7es14.4)') fi,r, tet (i,j), ft (i,j), ufi(i,j), ur(i,j) ,vi(i,j)
    enddo
enddo
    
open(24,file=trim('name2')//' ur '//' .dat')
do i=1,mfi
    fi=(i-1)*dfi
    write(24,'(6es14.4)') fi, ur(i,2), ur (i,mr/3), ur(i,mr/2), ur(i,mr/4), ur(i,mr) 
enddo 

end Program


!sweep method
Subroutine Tom(i0,iN,a,b,c,d,e,mi,ma)
    integer i0,iN,mi,ma,i
    real,dimension (mi:ma):: a,b,c,d,e,alf,bet
    alf(i0)=c(i0)/b(i0); bet(i0)=d(i0)/b(i0)
    do i=i0+1,iN-1
        alf(i)=c(i)/(b(i)-a(i)*alf(i-1))
    enddo
    
    do i=i0+1,iN
        bet(i)=(d(i)+a(i)*bet(i-1))/(b(i)-a(i)*alf(i-1))
    enddo
    
    e(iN)=bet(iN)
    do i=iN-1,i0,-1
        e(i)=alf(i)*e(i+1)+bet(i)
    enddo
end
