integer,parameter:: mr=201,mf=101,ma=max(mr,mf),mi=1
real, dimension(mr,mf)::tet,tet2,tets,tet1,vi,vi2,vis,vi1,ft,ft1,ft2,fts,ur,uf
real, dimension(mi:ma)::a,b,c,d,e
character(100):: sim ,name

namelist/Kanal_plosk/ Re,Pr, dt, dl_trub, dr, df, sumft, sumvi, sumtet

!-------------
itk = 400
ish = 100
it = 0

dt=0.0002
Re=10.
Pr=1.0

!channel geometry
r0=1.
r1=2.
pi=atan(1.)*4

dr=(r1-r0)/(mr-1)
df=pi /(mf-1)
dr2=dr*dr
df2=df*df

!temperature I.C.  
tet(mr,2:mf-1)=0.
tet(1,2:mf-1)=1.
tet(2:mr-1,1)=0.
tets=tet

!stream-function I.C.
do i=1,mr
    r=r0+(i-1)*dr
    ft(i,1)=log(r/r0)
enddo
ft(mr,2:mf)=ft(mr,1)
fts=ft

!velocities I.C.
do i=2,mr-1
    r=r0+(i-1)*dr
    uf(i,1)=1./r
enddo
    


!MAIN CYCLE          
do while (it<=itk)
    
    !--------------STREAM FUNCTION
    !sweep for x1
    do j=2, mf-1
        do i=2,mr-1
            r=r0+(i-1)*dr
            a(i)=1/dr2
            c(i)=1/dr2
            b(i)=1/dt+2/dr2
            d(i)=fts(i,j)/dt-vi(i,j)+1/r*(ft(i+1,j)-ft(i-1,j))/2/dr
        enddo !i
        a(1)=0.
        b(1)=1.
        c(1)=0.
        d(1)=fts(1,j)
        a(mr)=0
        b(mr)=1.
        c(mr)=0.
        d(mr)=fts(mr,j)
        call Tom(1,mr,a,b,c,d,e,mi,ma)
        ft1(1:mr,j)=e(1:mr)
    enddo

    !sweep for y1
    do i=2, mr-1
        do j=2,mf-1
            r=r0+(i-1)*dr
            a(j)=1/df2/r/r
            c(j)=1/df2/r/r
            b(j)=1/dt+2/df2/r/r
            d(j)=ft1(i,j)/dt
        enddo !i
        a(1)=0.
        b(1)=1.
        c(1)=0.
        d(1)=fts(i,1)
        a(mf)=1.
        b(mf)=1.
        c(mf)=0.
        d(mf)=0.
        call Tom(1,mf,a,b,c,d,e,mi,ma)
        ft(i,1:mf)=e(1:mf)
    enddo
    
    !local B.C. for stream function
    ft(2:mr-1,mf)=ft(2:mr-1,mf-1)
    ft2=ft-fts

    !B.C. for velocities
    do i=2, mr-1
        do j=2, mf-1
            r=r0+(i-1)*dr
            ur(i,j)=-(ft(i,j+1)-ft(i,j-1))/2/df/r
            uf(i,j)=(ft(i+1,j)-ft(i-1,j))/2/dr
        enddo
    enddo
    ur(2:mr-1,mf)=ur(2:mr-1,mf-1)
    uf(2:mr-1,mf)=uf(2:mr-1,mf-1)
    ur(2:mr-1,1)=ur(2:mr-1,2)

    !--------------VORTICITY
    !sweep for x1
    do j=2, mf-1
        do i=2,mr-1
            r=r0+(i-1)*dr
            aur=abs(ur(i,j))
            a(i)=(aur+ur(i,j))/2/dr+1/Re/dr2
            c(i)=(aur-ur(i,j))/2/dr+1./Re/dr2
            b(i)=1/dt+2/Re/dr2+aur/dr
            d(i)=vis(i,j)/dt+1/r/Re*(vis(i+1,j)-vis(i-1,j))/2/dr
        enddo !i
        a(1)=0.
        b(1)=1.
        c(1)=0.
        d(1)=2*(ft(2,j)-ft(1,j))/dr2
        a(mr)=1.
        b(mr)=1.
        c(mr)=0.
        d(mr)=2*(ft(mr-1,j)-ft(mr,j))/dr2
        call Tom(1,mr,a,b,c,d,e,mi,ma)
        vi1(1:mr,j)=e(1:mr)
    enddo
    
    !sweep for y1
    do i=2, mr-1
        do j=2,mf-1
            r=r0+(i-1)*dr
            auf=abs(uf(i,j))
            a(j)=(auf+uf(i,j))/2/df/r+1/Re/df2/r/r
            c(j)=(auf-uf(i,j))/2/df/r+1/df2/Re/r/r
            b(j)=2./(df2*re)/r/r+1./dt+auf/df/r
            d(j)=vi1(i,j)/dt
        enddo !i
        a(1)=0
        b(1)=1.
        c(1)=0
        d(1)=0.
        a(mf)=1.
        b(mf)=1.
        c(mf)=0.
        d(mf)=0.
        call Tom(1,mf,a,b,c,d,e,mi,ma)
        vi(i,1:mf)=e(1:mf)
    enddo
    
    !local B.C. for vorticity
    vi(2:mr-1,mf)=vi(2:mr-1,mf-1)
    do j=2,mf-1
        vi(1,j)=2*(ft(2,j)-ft(1,j))/dr2
        vi(mr,j)=2*(ft(mr-1,j)-ft(mr,j))/dr2
    enddo
    vi2=vi-vis
    
    
    !--------------TEMPERATURE
    do ll=1,3
        !sweep for x1
        do j=2, mf-1
            do i=2,mr-1
                r=r0+(i-1)*dr
                aur=abs(ur(i,j))
                a(i)=(aur+ur(i,j))/2/dr+1/Re/Pr/dr2
                c(i)=(aur-ur(i,j))/2/dr+1./Re/Pr/dr2
                b(i)=1/dt+2/Re/Pr/dr2+aur/dr
                d(i)=tets(i,j)/dt +1/r/Re*(tets(i+1,j)-tets(i-1,j))/2/dr/Pr
            enddo !iЗк
            a(1)=0.
            b(1)=1.
            c(1)=0.
            d(1)=1.
            a(mr)=0.
            b(mr)=1.
            c(mr)=0.
            d(mr)=0.
            call Tom(1,mr,a,b,c,d,e,mi,ma)    
            tet1(1:mr,j)=e(1:mr)
        enddo

        !sweep for y1
        do i=2, mr-1
            do j=2,mf-1
                r=r0+(i-1)*dr
                auf=abs(uf(i,j))
                a(j)=(auf+uf(i,j))/2/df/r+1/Re/Pr/df2/r/r
                c(j)=(auf-uf(i,j))/2/df/r+1/df2/Re/Pr/r/r
                b(j)=2./df2/re/Pr/r/r+1./dt+auf/df/r
                d(j)=tet1(i,j)/dt
            enddo !i
            a(1)=0
            b(1)=1.
            c(1)=0
            d(1)=0.
            a(mf)=1.
            b(mf)=1.
            c(mf)=0.
            d(mf)=0.
            call Tom(1,mf,a,b,c,d,e,mi,ma)
            tet(i,1:mf)=e(1:mf)
        enddo

        !local B.C. for temperature
        tet(mr,1:mf)=0.
        tet(1,1:mf)=1.
        tet(2:mr-1,1)=0.
        tet(2:mr-1,mf)=tet(2:mr-1,mf-1)
        tet2=tet-tets
        tets=tet
    enddo

    vis=vi
    fts=ft
    tets=tet

    sumft=0.
    sumvi=0.
    sumtet=0.
    
    !CONVERGENCE CHECK
    do i=2, mr-1
        do j=2, mf-1
            sumft= sumft+abs(ft2(i,j))/dt
            sumvi= sumvi+abs(vi2(i,j))/dt
            sumtet= sumtet+abs(tet2(i,j))/dt
        enddo
    enddo
    
    it = it+1
enddo !while


!--------------------- Write to files--------------------
open(23,file=trim('name1')//' izo_tem_ft_ur_uf_vi '//' .dat')
do i=2,mr-1
    do j=2,mf-1
        r=r0+(i-1)*dr
        f=(j-1)*df
        x=r*cos(f)
        y=r*sin(f)
        write(23,'(7es14.4)') x,y, tet (i,j), ft (i,j), ur(i,j), uf(i,j) ,vi(i,j)
    enddo
enddo

open(24,file=trim('name2')//' ur '//' .dat')
do j=1,mf
    f=(j-1)*df
    write(24,'(6es14.4)') f, ur(2,j), ur (4,j), ur(8,j), ur(16,j), ur(mr,j)
enddo
    
open(25,file=trim('name3')//' tet '//' .dat')
do j=1,mf
    f=(j-1)*df
    write(25,'(6es14.4)') f, tet(2,j), tet (4,j), tet(8,j), tet(16,j), tet(mr,j)
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