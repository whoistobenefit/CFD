integer,parameter:: mx=101,my=101,ma=max(mx,my),mi=1
real, dimension(1:mx,1:my)::tet,tets,tet1,tet2,vi,vis,vi1,vi2,ft,ft1,ft2,fts,ux,uy
real, dimension(mi:ma)::a,b,c,d,e
character(100):: sim ,name

!namelist/CHANNEL/ Re,Pr, dt, dl_trub, dx, dy, sumft, sumvi, sumtet

!-------------
itk = 400
ish = 100
it = 0

dt=0.0005
Re=100.
Pr=100.
Gr=0.

dl_trub = 1.

sumft = 0.
sumvi = 0.
sumtet = 0.
           
dy=1./(my-1)
dy2=dy*dy
dx=dl_trub/(mx-1)
dx2=dx*dx

!channel geometry
y0=3./10.
j0=nint(y0/dy)+1
y0=(j0-1)*dy
      
y1=4./10.
j1=nint(y1/dy)+1
y1=(j1-1)*dy
     
y2=6./10.
j2=nint(y2/dy)+1
y2=(j2-1)*dy
     
y3=7./10.
j3=nint(y3/dy)+1
y3=(j3-1)*dy

!temperature I.C.       
tet(1,2:j0)=0.
tet(1,j0+1:j1)=1.
tet(1,j1+1:my-1)=0.
tet(2:mx-1,my)=0.
tets=tet

!I.C. for velocities
ux(1,2:j0)=0.
ux(1,j0+1:j1)=1.
ux(1,j1+1:j2)=0.
ux(1,j2+1:j3)=1.
ux(1,j3+1:my-1)=0.
       
!stream-function I.C.
ft(1,1:j0)=0.
do j=j0+1,j1
    y=(j-1)*dy
    ft(1,j)=y-y0
enddo    
ft(1,j1+1:j2)=ft(1,j1)
do j=j2+1,j3
    y=(j-1)*dy
    ft(1,j)=y-y2+ft(1,j2)
enddo    
ft(1,j3+1:my)=ft(1,j3)
ft(2:mx,my)=ft(1,my)
ft(2:mx,1)=0.
fts=ft
  
   
!MAIN CYCLE
do while (it<=itk)
    
    !--------------STREAM FUNCTION
    !sweep for x1
    do j=2, my-1
        do i=2,mx-1       
            a(i)=1/dx2
            c(i)=1/dx2
            b(i)=1/dt+2./dx2
            d(i)=fts(i,j)/dt-vi(i,j)
        enddo
        a(1)=0.
        b(1)=1.
        c(1)=0.
        d(1)=fts(1,j)
        a(mx)=1.
        b(mx)=1.
        c(mx)=0.
        d(mx)=0.
        call Tom(1,mx,a,b,c,d,e,mi,ma)
        ft1(1:mx,j)=e(1:mx)
    enddo
           
    !sweep for y1
    do i=2, mx-1
        do j=2,my-1
            a(j)=1/dy2
            c(j)=1/dy2
            b(j)=1/dt+2/dy2
            d(j)=ft1(i,j)/dt
        enddo
        a(1)=0.
        b(1)=1.
        c(1)=0.
        d(1)=fts(i,1)
        a(my)=0.
        b(my)=1.
        c(my)=0.
        d(my)=fts(i,my)
        call Tom(1,my,a,b,c,d,e,mi,ma)
        ft(i,1:my)=e(1:my)
    enddo
    
    !local B.C. for stream function
    ft(mx,2:my-1)=ft(mx-1,2:my-1)
    ft2=ft-fts 
    
    !B.C. for velocities
    do i=2, mx-1
        do j=2, my-1
            ux(i,j)=(ft(i,j+1)-ft(i,j-1))/2/dy
            uy(i,j)=-(ft(i+1,j)-ft(i-1,j))/2/dx
        enddo
    enddo
    ux(mx,2:my-1)=ux(mx-1,2:my-1)
    uy(mx,2:my-1)=uy(mx-1,2:my-1)
    uy(1,j0+1:j1)=uy(2,j0+1:j1)
    uy(1,j2+1:j3)=uy(2,j2+1:j3)
   
    
    !--------------VORTICITY
    !sweep for x1
    do j=2, my-1
        do i=2,mx-1
            aux=abs(ux(i,j))
            a(i)=(aux+ux(i,j))/2/dx+1/Re/dx2
            c(i)=(aux-ux(i,j))/2/dx+1./Re/dx2
            b(i)=1/dt+2/Re/dx2+aux/dx            
            d(i)=vis(i,j)/dt-(tet(i+1,j)-tet(i-1,j))/2/dx*Gr/Re/Re
        enddo !i
        
        if( j>=2 .and. j<=j0) then
            a(1)=0.
            b(1)=1.
            c(1)=0.
            d(1)=2.*(ft(2,j)-ft(1,j))/dx2
        endif
        if (j>=j0+1 .and. j<=j1) then
            a(1)=0.
            b(1)=1.
            c(1)=0.
            d(1)=0.
        endif
        if (j>=j1+1 .and. j<=j2) then
            a(1)=0.
            b(1)=1.
            c(1)=0.
            d(1)=2.*(ft(2,j)-ft(1,j))/dx2
        endif
        if (j>=j2+1 .and. j<=j3) then
            a(1)=0.
            b(1)=1.
            c(1)=0.
            d(1)=0.
        endif
        if (j>=j3+1 .and. j<=my-1) then
            a(1)=0.
            b(1)=1.
            c(1)=0.
            d(1)=2.*(ft(2,j)-ft(1,j))/dx2
        endif
        a(mx)=1.
        b(mx)=1.
        c(mx)=0.
        d(mx)=0.
        call Tom(1,mx,a,b,c,d,e,mi,ma)
        vi1(1:mx,j)=e(1:mx)
    enddo
    
    !sweep for y1
    do i=2, mx-1
        do j=2,my-1
            auy=abs(uy(i,j))
            a(j)=(auy+uy(i,j))/2/dy+1/Re/dy2
            c(j)=(auy-uy(i,j))/2/dy+1/dy2/Re
            b(j)=2./(dy2*re)+1./dt+auy/dy
            d(j)=vi1(i,j)/dt
        enddo !i             
        a(1)=0.
        b(1)=1.
        c(1)=0.
        d(1)=2*(ft(i,2)-ft(i,1))/dy2             
        a(my)=0.
        b(my)=1.
        c(my)=0.
        d(my)=2*(ft(i,my-1)-ft(i,my))/dy2
        call Tom(1,my,a,b,c,d,e,mi,ma)
        vi(i,1:my)=e(1:my)
    enddo
        
    !local B.C. for vorticity        
    vi(mx,2:my-1)=vi(mx-1,2:my-1)
    vi(1,2:j0)=2.*(ft(2,2:j0)-ft(1,2:j0))/dx2  
    vi(1,j1+1:j2)=2.*(ft(2,j1+1:j2)-ft(1,j1+1:j2))/dx2
    vi(1,j3:my-1)=2.*(ft(2,j3:my-1)-ft(1,j3:my-1))/dx2
    vi2=vi-vis
       
    vis=vi
    fts=ft
       
    !--------------TEMPERATURE
    !sweep for x1              
    do ll=1,3
        do j=2, my-1
            do i=2,mx-1
                aux=abs(ux(i,j))
                a(i)=(aux+ux(i,j))/2/dx+1/Re/Pr/dx2
                c(i)=(aux-ux(i,j))/2/dx+1./Re/Pr/dx2
                b(i)=1/dt+2/Re/Pr/dx2+aux/dx
                d(i)=tets(i,j)/dt
            enddo !i
            if( j>=2 .and. j<=j0) then
                a(1)=0.
                b(1)=1.
                c(1)=1.
                d(1)=0.
            endif      
            if (j>=j0+1 .and. j<=j1) then
                a(1)=0.
                b(1)=1.
                c(1)=0.
                d(1)=1.
            endif      
            if (j>=j1+1 .and. j<=j2) then
                a(1)=0.
                b(1)=1.
                c(1)=1.
                d(1)=0.
            endif      
            if (j>=j2+1 .and. j<=j3) then
                a(1)=0.
                b(1)=1.
                c(1)=0.
                d(1)=0.
            endif      
            if (j>=j3+1 .and. j<=my-1) then
                a(1)=0.
                b(1)=1.
                c(1)=1.
                d(1)=0.
            endif 
            a(mx)=1.
            b(mx)=1.
            c(mx)=0.
            d(mx)=0.
            call Tom(1,mx,a,b,c,d,e,mi,ma)
            tet1(1:mx,j)=e(1:mx)
        enddo
    
        !sweep for y1
        do i=2, mx-1
            do j=2,my-1
                auy=abs(uy(i,j))
                a(j)=(auy+uy(i,j))/2/dy+1/Re/Pr/dy2
                c(j)=(auy-uy(i,j))/2/dy+1/dy2/Re/Pr
                b(j)=2./dy2/re/Pr+1./dt+auy/dy
                d(j)=tet1(i,j)/dt
            enddo !i
            a(1)=0.
            b(1)=1.
            c(1)=1.
            d(1)=0. 
            a(my)=1.
            b(my)=1.
            c(my)=0.
            d(my)=0.
            call Tom(1,my,a,b,c,d,e,mi,ma)
            tet(i,1:my)=e(1:my)
        enddo 
        
        !local B.C. for temperature
        tet(mx,2:my-1)=tet(mx-1,2:my-1)
        tet(2:mx-1,1)=tet(2:mx-1,2)
        tet(1,j0+1:j1)=tet(2,j0+1:j1)
        tet(1,j1+1:j2)=tet(2,j1+1:j2)
        tet(1,j3+1:my-1)=tet(2,j3+1:my-1)
        tet2=tet-tets
        tets=tet
    enddo

    !CONVERGENCE CHECK
    do i = 2, mx - 1
        do j = 2, my - 1
            sumft = sumft + abs(ft2(i, j)) / dt
            sumvi = sumvi + abs(vi2(i, j)) / dt
    		sumtet = sumtet + abs(tet2(i, j)) / dt
        enddo
    enddo
    
    if (it==itk-1) then
        write(*, '(I7,3es14.4)') it, sumft, sumvi, sumtet
    endif

    it = it+1
enddo !while


!--------------------- Write to files--------------------
open(23,file=trim('name1')//'grid_psi_temp'//' .dat')
    do i=1,mx
        do j=1,my
		    x=(i-1)*dx
            y=(j-1)*dy
		    write(23,'(7es14.4)') x,y,  ft (i,j), tet(i,j)
        enddo
    enddo
            
            
    open(52,file=trim('name2')//'grid_vels_vort'//' .dat')
    do i=1,mx
        do j=1,my
		    x=(i-1)*dx
            y=(j-1)*dy
            write(52,'(7es14.4)') x,y, ux(i,j), uy(i,j), vi(i,j)
        enddo
    enddo  

end Program


!sweep function
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


