
integer,parameter:: mx = 100, my = 150, ma = max(mx, my), mi = 1
real, dimension(1:mx, 1 : my)::tet, tets, tet1, tet2, vi, vis, vi1, vi2, ft, ft1, ft2, fts, ux, uy
real, dimension(mi:ma)::a, b, c, d, e
character(100)::sim, name

namelist / HEAT / Re, Pr, dt, Hz, Dl, dx, dy, sumft, sumvi, sumtet


!-------------
itk = 400
ish = 100
it = 0

dt = 0.0025 
Re = 20.
Pr = 1.0

!channel geometry
Hz = 5. 
Dl = 3.
x1 = 2.
y1 = 1.

sumft = 0.
sumvi = 0.
sumtet = 0.

dy = Hz / (my - 1)
dy2 = dy*dy

dx = Dl / (mx - 1)
dx2 = dx*dx

i1 = nint(x1 / dx) + 1
j1 = nint(y1 / dy) + 1

! temperature I.C.
tet(1, 2:j1) = 0.
tet(mx, 1:my) = 1.
tets = tet

! longitudinal velocity I.C.
ux(1, 2:j1-1) = 1.

! stream-function I.C.
do j = 1, j1
    y = (j - 1)*dy
    ft(1, j) = y
enddo
ft(i1+1:mx-1, my) = ft(i1+1:mx-1, my-1)


!MAIN CYCLE
do while (it <= itk)

    !--------------STREAM FUNCTION
    !sweep for x1
    do j = 2, j1
        do i = 2, mx-1
            a(i) = 1 / dx2
            c(i) = 1 / dx2
            b(i) = 1 / dt + 2. / dx2
            d(i) = fts(i, j) / dt - vi(i, j)
        enddo !i
        a(1) = 0.
        b(1) = 1.
        c(1) = 0.
        d(1) = fts(1, j)
        a(mx) = 0.
        b(mx) = 1.
        c(mx) = 0.
        d(mx) = 0.
        call Tom(1, mx, a, b, c, d, e, mi, ma)
        ft1(1:mx, j) = e(1:mx)
    enddo
	   
    !sweep for x2
    do j = j1+1, my-1
        do i = i1+1, mx - 1
            a(i) = 1 / dx2
	        c(i) = 1 / dx2
	        b(i) = 1 / dt + 2. / dx2
	        d(i) = fts(i, j) / dt - vi(i, j)
        enddo !i
	    a(i1) = 0.
        b(i1) = 1.
        c(i1) = 0.
        d(i1) = fts(1, j1)
        a(mx) = 0.
        b(mx) = 1.
        c(mx) = 0.
        d(mx) = 0.
	    call Tom(1, mx, a, b, c, d, e, mi, ma)
		ft1(1:mx, j) = e(1:mx)
	enddo

    !sweep for y1
    do i = 2, i1
        do j = 2, j1-1
            a(j) = 1 / dy2
			c(j) = 1 / dy2
            b(j) = 1 / dt + 2 / dy2
            d(j) = ft1(i, j) / dt
		enddo !i
        a(1) = 0.
        b(1) = 1.
        c(1) = 0.
        d(1) = 0.
	    a(j1) = 0.
        b(j1) = 1.
        c(j1) = 0.
        d(j1) = fts(1, j1)
        call Tom(1, j1, a, b, c, d, e, mi, ma)
        ft(i, 1:j1) = e(1:j1)
	enddo

    !sweep for y2
    do i = i1 + 1, mx-1
        do j = 2 , my-1
            a(j) = 1 / dy2
            c(j) = 1 / dy2
	        b(j) = 1 / dt + 2 / dy2
            d(j) = ft1(i, j) / dt
		enddo !i
        a(1) = 0.
        b(1) = 1.
        c(1) = 0.
        d(1) = 0.
        a(my) =1.
        b(my) = 1.
        c(my) = 0.
        d(my) = 0
        call Tom(1, my, a, b, c, d, e, mi, ma)
        ft(i, 1:my) = e(1:my)
    enddo !i
    
    ! local B.C.
    ft(i1:mx-1, my) = ft(i1:mx-1, my-1)
    ft2 = ft - fts
    
    ! B.C. for velocities
    do i = 2, mx - 1
        do j = 2, my - 1
            ux(i, j) = (ft(i, j + 1) - ft(i, j - 1)) / 2 / dy
            uy(i, j) = -(ft(i + 1, j) - ft(i - 1, j)) / 2 / dx
    	enddo
    enddo
    ux(i1:mx-1, my) = ux(i1:mx-1, my-1)
    uy(i1:mx-1, my) = uy(i1:mx-1, my-1)
    
    
    !--------------VORTICITY
    !sweep for x1
    do j = 2, j1
        do i = 2, mx - 1
            aux = abs(ux(i, j))							
            a(i) = (aux + ux(i, j)) / 2 / dx + 1 / Re / dx2
            c(i) = (aux - ux(i, j)) / 2 / dx + 1. / Re / dx2
            b(i) = 1 / dt + 2 / Re / dx2 + aux / dx
            d(i) = vis(i, j) / dt
        enddo !i
    	a(1) = 0.
    	b(1) = 1.
    	c(1) = 0.
    	d(1)= 0.
    	a(mx) = 0.
    	b(mx) = 1.
    	c(mx)= 0.
    	d(mx)= 2*(ft(mx-1, j) - ft(mx, j)) / dx2
    	call Tom(1, mx, a, b, c, d, e, mi, ma)
        vi1(1:mx, j) = e(1:mx)
    enddo
    			
    !sweep for x2			
    do j = j1 + 1, my
        do i = i1 + 1, mx - 1
            aux = abs(ux(i, j))
    		a(i) = (aux + ux(i, j)) / 2 / dx + 1 / Re / dx2
    		c(i) = (aux - ux(i, j)) / 2 / dx + 1. / Re / dx2
    		b(i) = 1 / dt + 2 / Re / dx2 + aux / dx
    		d(i) = vis(i, j) / dt
        enddo !i
        a(i1) = 0.
        b(i1) = 1.
        c(i1) = 0.
        d(i1) = 2*(ft(i1+1, j) - ft(i1, j)) / dx2
    	a(mx) = 0.
    	b(mx) = 1.
    	c(mx) = 0.
    	d(mx) = 2*(ft(mx-1, j) - ft(mx, j)) / dx2
    	call Tom(1, mx, a, b, c, d, e, mi, ma)
        vi1(1:mx, j) = e(1:mx)
    enddo
    
    
    !sweep for y1
    do i = 2, i1
        do j = 2, j1 - 1
            auy = abs(uy(i, j))
    		a(j) = (auy + uy(i, j)) / 2 / dy + 1 / Re / dy2
            c(j) = (auy - uy(i, j)) / 2 / dy + 1 / dy2 / Re
            b(j) = 2. / (dy2*re) + 1. / dt + auy / dy
            d(j) = vi1(i, j) / dt
    	enddo !i
    	a(1) = 0.
    	b(1) = 1.
    	c(1) = 0.
    	d(1) = 2*(ft(2, j) - ft(1, j)) / dy2
    	a(j1) = 0.
    	b(j1) = 1.
    	c(j1) = 0.
    	d(j1) = 2*(ft(j1-1, j) - ft(j1, j)) / dy2
    	call Tom(1, my, a, b, c, d, e, mi, ma)
        vi(i, 1:j1) = e(1:j1)
    enddo
    
    !sweep for y2					 
    do i = i1 + 1, mx - 1
        do j = 2, my
            auy = abs(uy(i, j))
            a(j) = (auy + uy(i, j)) / 2 / dy + 1 / Re / dy2
            c(j) = (auy - uy(i, j)) / 2 / dy + 1 / dy2 / Re
    		b(j) = 2. / (dy2*re) + 1. / dt + auy / dy
            d(j) = vi1(i, j) / dt
    	enddo !i
    	a(1) = 0.
    	b(1) = 1.
    	c(1) = 0.
    	d(1) = 2*(ft(2, j) - ft(1, j)) / dy2
        a(my) = 1.
        b(my) = 1.
        c(my) = 0.
        d(my) = 0
        call Tom(1, my, a, b, c, d, e, mi, ma)
        vi(i, 1:my) = e(1:my)
    enddo
    						 
    ! local B.C. for vorticity
    vi(i1:mx-1, my) = vi(i1:mx-1, my-1)
    vi2 = vi - vis
    
    	
    !--------------TEMPERATURE
    !sweep for x1
    do ll = 1, 3
    	do j = 2, j1
    		do i = 2, mx - 1 
    		    aux = abs(ux(i, j))
    		    a(i) = (aux + ux(i, j)) / 2 / dx + 1 / Re / Pr / dx2
    			c(i) = (aux - ux(i, j)) / 2 / dx + 1. / Re / Pr / dx2
    			b(i) = 1 / dt + 2 / Re / Pr / dx2 + aux / dx
    			d(i) = tets(i, j) / dt
    		enddo !i
    	    a(1) = 0.
    	    b(1) = 1.
    	    c(1) = 0.
    	    d(1) = 0.
    	    a(mx) = 0.
    	    b(mx) = 1.
    	    c(mx) = 0.
    	    d(mx) = 1.
    	    call Tom(1, mx, a, b, c, d, e, mi, ma)
    	    tet1(1:mx, j) = e(1:mx)
        enddo
    
    	!sweep for x2
        do j = j1 + 1, my
            do i = i1 + 1, mx - 1
                aux = abs(ux(i, j))
    			a(i) = (aux + ux(i, j)) / 2 / dx + 1 / Re / Pr / dx2
                c(i) = (aux - ux(i, j)) / 2 / dx + 1. / Re / Pr / dx2
    			b(i) = 1 / dt + 2 / Re / Pr / dx2 + aux / dx
    			d(i) = tets(i, j) / dt
    		enddo !i
    		a(i1) = 0.
    		b(i1) = 1.
    		c(i1) = 1.
    		d(i1) = 0.
    		a(mx) = 0.
    		b(mx) = 1.
    		c(mx) = 0.
    		d(mx) = 1.
            call Tom(1, mx, a, b, c, d, e, mi, ma)
    		tet1(1:mx, j) = e(1:mx)
    	enddo
    	
    	
    	!sweep for y1
        do i = 2, i1
            do j = 2, j1 - 1
    		    auy = abs(uy(i, j))
                a(j) = (auy + uy(i, j)) / 2 / dy + 1 / Re / Pr / dy2
    		    c(j) = (auy - uy(i, j)) / 2 / dy + 1 / dy2 / Re / Pr
    		    b(j) = 2. / dy2 / re / Pr + 1. / dt + auy / dy
    		    d(j) = tet1(i, j) / dt
            enddo !i
    		a(1) = 0.
    		b(1) = 1.
    		c(1) = 1.
    		d(1) = 0.
    		a(j1) = 1.
    		b(j1) = 1.
    		c(j1) = 0.
    		d(j1) = 0.
    		call Tom(1, my, a, b, c, d, e, mi, ma)
    	    tet(i, 1:j1) = e(1:j1)
    	enddo
    
        !sweep for y2
        do i = i1 + 1, mx - 1
            do j = 2, my 
                auy = abs(uy(i, j))
    		    a(j) = (auy + uy(i, j)) / 2 / dy + 1 / Re / Pr / dy2
    		    c(j) = (auy - uy(i, j)) / 2 / dy + 1 / dy2 / Re / Pr
    		    b(j) = 2. / dy2 / re / Pr + 1. / dt + auy / dy
    		    d(j) = tet1(i, j) / dt
    		enddo !i
    	    a(1) = 0.
    	    b(1) = 1.
    	    c(1) = 1.
    	    d(1) = 0.
    		a(my) = 1.
    		b(my) = 1.
    		c(my) = 0.
    		d(my) = 0
    		call Tom(1, my, a, b, c, d, e, mi, ma)
    		tet(i, 1:my) = e(1:my)
    	enddo
    
    
        ! local B.C. for temperature
        tet(mx, 2:my-1) = tet(mx - 1, 2:my-1)
    	tet(2:mx-1, 1) = tet(2:mx-1, 2)
    	tet2 = tet - tets
    	tets = tet
    enddo
    

    vis = vi
    fts = ft
    tets = tet
    
    
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
		    write(23,'(7es14.4)') x,y,  ft (i,j)!, tet(i,j)
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


																