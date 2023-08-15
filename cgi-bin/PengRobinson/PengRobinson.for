c     Cálculo da Pressão de Vapor de Compostos Puros 
c     Equação de Peng-Robinson

c     Declaração das variáveis usadas no programa

      implicit none
      character *200 filename
	integer maxit, options, opcao
	parameter (maxit=1000)
      double precision a,b,at,bt,rg,p,pc,pvap,te,tc,w,fl,fv,zl,zv,alfa,
     *a0,a1,a2,a3,erro,k,z1,z2,z3,aq,bq,cq,delta,tol,x(maxit),
     *f(maxit),df(maxit),tr,pr,da,hrl,hrv,hgi,hl,hv,acp,bcp,ccp,dcp,
     *sgi,srl,srv,sl,sv,vl,vv
      integer it,i,n,novo
      character(len=20) :: data
      
c     Dados da substƒncia

c      write(*,*) 'Entre com Tc (K), Pc (bar), w'
c      read(*,*) tc,pc,w
      call getarg(1,data)
      read(data,*) tc
      call getarg(2,data)
      read(data,*) pc
      pc = pc*1.d5
      call getarg(3,data)
      read(data,*) w

c      write(*,*) 'Digite as constantes para o calculo de Cp (J/mol.K)'
c      read(*,*) acp,bcp,ccp,dcp
      call getarg(4,data)
      read(data,*) acp
      call getarg(5,data)
      read(data,*) bcp
      call getarg(6,data)
      read(data,*) ccp
      call getarg(7,data)
      read(data,*) dcp

       
c      tc = 304.1
c      	pc = 73.8*1.d5
c	w = .239

c	acp = 19.8jkkk
c	bcp = .07344
c	ccp = -5.602d-5
c	dcp = 1.715d-8

c	write(*,*) 'Entre com T(K) e P(bar) de referencia'
c	read(*,*) tr,pr
      call getarg(8,data)
      read(data,*) tr
      call getarg(9,data)
      read(data,*) pr
	pr=pr*1.d5

c	tr = 298.15
c	pr = 1.d5
	
c      write(*,*) 'Para calcular dados de equil¡brio digite 1'
c      write(*,*) 'Para calcular dados de em qualquer P e T digite 2'
c      read(*,*) opcao

      call getarg(10,data)
      read(data,*) te

      call getarg(11,data)
      read(data,*) p
      p = p*1.d5

      call getarg(12,data)
      read(data,*) opcao

      options = 0
      if ( options == 1 ) then
         write(*,*) " -------------------------------------------- <br>"
         write(*,*) " -------------------------------------------- <br>"
         write(*,*) " -------------------------------------------- <br>"
         write(*,*) " Dados de entrada: <br>"
         write(*,*) " tc = ", tc, "<br>"
         write(*,*) " pc = ", pc, "<br>"
         write(*,*) " w = ", w, "<br>"
         write(*,*) " acp = ", acp, "<br>"
         write(*,*) " bcp = ", bcp, "<br>"
         write(*,*) " ccp = ", ccp, "<br>"
         write(*,*) " dcp = ", dcp, "<br>"
         write(*,*) " tr = ", tr, "<br>"
         write(*,*) " pr = ", pr, "<br>"
         write(*,*) " te = ", te, "<br>"
         write(*,*) " p = ", p, "<br>"
         write(*,*) " tipo de calculo (opcao) = ", opcao, "<br>"
         write(*,*) " -------------------------------------------- <br>"
         write(*,*) " -------------------------------------------- <br>"
         write(*,*) " -------------------------------------------- <br>"
      end if

      if (opcao.eq.1) then

c1000  write(*,*) 'C lculo das propriedades de l¡quido e vapor saturados'

	if (te.gt.tc) then
	write(*,*) 'Atencao!!! Temperatura acima da critica - nao ha
     *pressao de vapor. Escolha outra temperatura'
      stop
	end if

5     rg = 8.314
      k = .37464 + 1.54226*w - .26992*(w**2.d0)
	alfa = (1.d0 + k*(1.d0-(te/tc)**.5))**2.d0
	at = .45724*(rg**2.d0)*(tc**2.d0)*alfa/pc
	bt = .0778*rg*tc/pc
      
      it = 0

10    it = it + 1
c      write (*,*) 'Iteracao numero ',it
      a = at*p/((rg**2.d0)*(te**2.d0))
	b = bt*p/(rg*te)

c     Resolucao da equacao de estado cubica (z**3+a2*z**2+a1*z+a0)

c     Coeficientes

      a3 = 1.d0
	a2 = -(1.d0-b)
	a1 = a-(3.d0*(b**2.d0))-2*b
	a0 = -(a*b-(b**2.d0)-(b**3.d0))
	
      tol = 1.d-9
	n = 1000
	
	x(1) = 1.d0
	f(1) = (a3*x(1)**3)+(a2*x(1)**2)+a1*x(1)+a0
	df(1) = (3.d0*a3*x(1)**2)+2.d0*a2*x(1)+a1

c     Método de Newton - cálculo da primeira raiz

	do i=1,n

	  x(i+1) = x(i) - (f(i)/df(i))
        f(i+1) = (a3*x(i+1)**3)+(a2*x(i+1)**2)+a1*x(i+1)+a0
	  df(i+1) = (3.d0*a3*x(i+1)**2)+2.d0*a2*x(i+1)+a1

        if (dabs(f(i+1)).lt.tol) goto 100
	  
      enddo

      if (i.eq.n) then
	write(*,*) 'Numero de iteracoes maximo - calculo nao converge'
	stop
	endif

100   z1 = x(i+1)
c      write(*,*) 'z1 = ',z1

c     Cálculo das outras 2 raízes (equação de segundo grau)

	aq = a3
	bq = a3*z1 + a2
	cq = a3*(z1**2) + a2*z1 + a1
	delta = bq**2 - 4*aq*cq
	
	if (dabs(delta).lt.1.d-15) then
	  delta = 0
	  z2 = z1
	  z3 = z1
	  write(*,*) 'Ponto Critico'
	  goto 300
	endif 
	
	if (delta.lt.0) then
	  z2 = z1
	  z3 = z1
c	  write(*,*) 'REGIAO DE 1 FASE - MELHORE O CHUTE INICIAL'
	  if (p.lt.pc) then 
	    p = p + 1.d5
	    goto 5
c	  else
c	    p = p - 1.d5
c         goto 5
	  endif
	endif
	

	z2 = (-bq - delta**.5)/(2*aq)
	z3 = (-bq + delta**.5)/(2*aq)


c      write(*,*) 'z2 = ',z2
c	write(*,*) 'z3 = ',z3
c	write(*,*)

c     Identificação da menor e da maior raiz

	zl = dmin1(z1,z2)
	zl = dmin1(zl,z3)

	zv = dmax1(z1,z2)
	zv = dmax1(zv,z3)

	vl = zl*rg*te/p
	vv = zv*rg*te/p

c     Cálculo das fugacidades - devem ser iguais

      fv = p*exp((zv-1.d0)-dlog(zv-b)-(a/(2.d0*(2.d0**.5)*b))*
     *dlog((zv+(1.d0+2.d0**.5)*b)/(zv+(1.d0-2.d0**.5)*b)))


      fl = p*exp((zl-1.d0)-dlog(zl-b)-(a/(2.d0*(2.d0**.5)*b))*
     *dlog((zl+(1.d0+2.d0**.5)*b)/(zl+(1.d0-2.d0**.5)*b)))
      
      erro = dabs((fl/fv)-1.d0)

	if (erro.ge.1.d-4) then

	p = (p*(fl/fv)+p)/2.d0	

c	write (*,*) 'Nova estimativa de Pvap: ',p*1.d-5,' Pa <br>'
	
	goto 10
	
	endif

c	Calculo das entalpias

	da = -.45724*(rg**2)*(tc**2)*k*((alfa/(te*tc))**.5)/pc

	hgi = acp*(te-tr) + (bcp/2)*(te**2 - tr**2) + 
     *(ccp/3)*(te**3 - tr**3) + (dcp/4)*(te**4 - tr**4)
     
      write(*,*) 'H(gas ideal)=',hgi,'J/mol <br>'

	hrl = rg*te*(zl-1.)+((te*da-at)/(2*(2**.5)*bt))*
     *dlog((zl+(1.+(2**.5))*b)/(zl+(1.-(2**.5))*b))

	hl = hgi + hrl

	write(*,*) 'Hr(liquido) =',hrl,'J/mol <br>'
	write(*,*) 'H(liquido) =',hl,'J/mol <br>'

	hrv = rg*te*(zv-1.)+((te*da-at)/(2*(2**.5)*bt))*
     *dlog((zv+(1.+(2**.5))*b)/(zv+(1.-(2**.5))*b))

	hv = hgi + hrv

      write(*,*) 'Hr(vapor) =',hrv,'J/mol <br>'
	write(*,*) 'H(vapor) =',hv,'J/mol <br>'
	write(*,*) "<br>"

c	Calculo das entropias

	sgi = acp*dlog(te/tr) + bcp*(te - tr) + (ccp/2)*(te**2 - tr**2)+
     *(dcp/3)*(te**3 - tr**3) - rg*dlog(p/pr)

	srl = rg*dlog(zl-b)+da/(2*(2**.5)*bt)*
     *dlog((zl+(1.+(2**.5))*b)/(zl+(1.-(2**.5))*b))

	sl = sgi + srl
	
	write(*,*) 'S(gas ideal) =',sgi,'J/mol.K <br>'

	write(*,*) 'Sr(liquido) =',srl,'J/mol.K <br>'
	write(*,*) 'S(liquido) =',sl,'J/mol.K <br>'
	write(*,*) "<br>"

	srv = rg*dlog(zv-b)+da/(2*(2**.5)*bt)*
     *dlog((zv+(1.+(2**.5))*b)/(zv+(1.-(2**.5))*b))

	sv = sgi + srv


      write(*,*) 'Sr(vapor) =',srv,'J/mol.K <br>'
	write(*,*) 'S(vapor) =',sv,'J/mol.K <br>'
	write(*,*) "<br>"

c     Resultado final      

	write(*,*) 'zv = ',zv, "<br>"
	write(*,*) 'zl = ',zl, "<br>"
      write(*,*) "<br>"

      write(*,*) 'vv = ',vv,'m3/mol <br>'
	write(*,*) 'vl = ',vl,'m3/mol <br>'
      write(*,*) "<br>"

      write(*,*) 'fv = ',fv*1.d-5,'bar <br>'
      write (*,*) 'fl = ',fl*1.d-5,'bar <br>'
      write(*,*) "<br>"

300	pvap = p 

	write(*,*) 'A pressao de vapor e de ',p*1.d-5, 'bar <br>'

c	write(*,*) 'Se quiser fazer novo calculo digite 1'
c	write(*,*) 'Se quiser calcular dados em outra P e T digite 2'
      stop
	
c opcao
      end if
      
      

      if (opcao.eq.2) then

c2000  write(*,*) 'Cálculo das propriedades da substância'

c      write(*,*) 'Digite a pressÆo (bar) e a temperatura (K)'
c      read (*,*) p,te
c      p = p*1.d5
      

      rg = 8.314
      k = .37464 + 1.54226*w - .26992*(w**2.d0)
	alfa = (1.d0 + k*(1.d0-(te/tc)**.5))**2.d0
	at = .45724*(rg**2.d0)*(tc**2.d0)*alfa/pc
	bt = .0778*rg*tc/pc
	
	a = at*p/((rg**2.d0)*(te**2.d0))
	b = bt*p/(rg*te)

c     Resolucao da equacao de estado cubica (z**3+a2*z**2+a1*z+a0)

c     Coeficientes

      a3 = 1.d0
	a2 = -(1.d0-b)
	a1 = a-(3.d0*(b**2.d0))-2*b
	a0 = -(a*b-(b**2.d0)-(b**3.d0))

      tol = 1.d-9
	n = 1000

	x(1) = 1.d0
	f(1) = (a3*x(1)**3)+(a2*x(1)**2)+a1*x(1)+a0
	df(1) = (3.d0*a3*x(1)**2)+2.d0*a2*x(1)+a1

c     Método de Newton - cálculo da primeira raiz

	do i=1,n

	  x(i+1) = x(i) - (f(i)/df(i))
        f(i+1) = (a3*x(i+1)**3)+(a2*x(i+1)**2)+a1*x(i+1)+a0
	  df(i+1) = (3.d0*a3*x(i+1)**2)+2.d0*a2*x(i+1)+a1

        if (dabs(f(i+1)).lt.tol) goto 400

      enddo

      if (i.eq.n) then
	write(*,*) 'Numero de iteracoes maximo - calculo nao converge'
	stop
	endif

400   z1 = x(i+1)
c      write(*,*) 'z1 = ',z1

c     Cálculo das outras 2 raízes (equação de segundo grau)

	aq = a3
	bq = a3*z1 + a2
	cq = a3*(z1**2) + a2*z1 + a1
	delta = bq**2 - 4*aq*cq

	if (dabs(delta).lt.1.d-15) then
	  delta = 0
	  z2 = z1
	  z3 = z1
	  write(*,*) 'Ponto Critico'
	  goto 600
	endif

	if (delta.lt.0) then
	  z2 = z1
	  z3 = z1
          goto 600
	endif


	z2 = (-bq - delta**.5)/(2*aq)
	z3 = (-bq + delta**.5)/(2*aq)


c      write(*,*) 'z2 = ',z2
c	write(*,*) 'z3 = ',z3
c	write(*,*)

c     Identificação da menor e da maior raiz

600	zl = dmin1(z1,z2)
	zl = dmin1(zl,z3)

	zv = dmax1(z1,z2)
	zv = dmax1(zv,z3)

	write(*,*) 'zv = ',zv, "<br>"
	write(*,*) 'zl = ',zl , "<br>"
	write(*,*) "<br>"
	
	vl = zl*rg*te/p
	vv = zv*rg*te/p

        write(*,*) 'vv = ',vv,'m3/mol <br>'
	write(*,*) 'vl = ',vl,'m3/mol <br>'
	write(*,*) "<br>"

c     Cálculo das fugacidades

      fv = p*exp((zv-1.d0)-dlog(zv-b)-(a/(2.d0*(2.d0**.5)*b))*
     *dlog((zv+(1.d0+2.d0**.5)*b)/(zv+(1.d0-2.d0**.5)*b)))

      write(*,*) 'fv = ',fv*1.d-5,'bar <br>'

      fl = p*exp((zl-1.d0)-dlog(zl-b)-(a/(2.d0*(2.d0**.5)*b))*
     *dlog((zl+(1.d0+2.d0**.5)*b)/(zl+(1.d0-2.d0**.5)*b)))

      write (*,*) 'fl = ',fl*1.d-5,'bar <br>'
	write(*,*) "<br>"

c	Calculo das entalpias

	da = -.45724*(rg**2)*(tc**2)*k*((alfa/(te*tc))**.5)/pc

	hgi = acp*(te-tr) + (bcp/2)*(te**2 - tr**2) +
     *(ccp/3)*(te**3 - tr**3) + (dcp/4)*(te**4 - tr**4)

 	hrl = rg*te*(zl-1.)+((te*da-at)/(2*(2**.5)*bt))*
     *dlog((zl+(1.+(2**.5))*b)/(zl+(1.-(2**.5))*b))

	hl = hgi + hrl
	
	write(*,*) 'H(gas ideal) =',hgi,'J/mol <br>'

	write(*,*) 'Hr(liquido) =',hrl,'J/mol <br>'
	write(*,*) 'H(liquido) =',hl,'J/mol <br>'
	write(*,*) "<br>"

	hrv = rg*te*(zv-1.)+((te*da-at)/(2*(2**.5)*bt))*
     *dlog((zv+(1.+(2**.5))*b)/(zv+(1.-(2**.5))*b))
     
	hv = hgi + hrv

	write(*,*) 'Hr(vapor) =',hrv,'J/mol <br>'
	write(*,*) 'H(vapor) =',hv,'J/mol <br>'
	write(*,*) "<br>"



c	Calculo das entropias

	sgi = acp*dlog(te/tr) + bcp*(te - tr) + (ccp/2)*(te**2 - tr**2)+
     *(dcp/3)*(te**3 - tr**3) - rg*dlog(p/pr)

      srl = rg*dlog(zl-b)+da/(2*(2**.5)*bt)*
     *dlog((zl+(1.+(2**.5))*b)/(zl+(1.-(2**.5))*b))

	sl = sgi + srl
	
	write(*,*) 'S(gas ideal) =',sgi,'J/mol.K <br>'

	write(*,*) 'Sr(liquido) =',srl,'J/mol.K <br>'
	write(*,*) 'S(liquido) =',sl,'J/mol.K <br>'
	write(*,*) "<br>"

	srv = rg*dlog(zv-b)+da/(2*(2**.5)*bt)*
     *dlog((zv+(1.+(2**.5))*b)/(zv+(1.-(2**.5))*b))

	sv = sgi + srv

	write(*,*) 'Sr(vapor) =',srv,'J/mol.K <br>'
	write(*,*) 'S(vapor) =',sv,'J/mol.K <br>'
	write(*,*) "<br>"

c      read (*,*)
      
c      write(*,*) 'Se quiser fazer novo calculo em outra P e T digite 1'
c      write(*,*) 'Se quiser calcular outra pressao de vapor digite 2'
c        read(*,*) novo
c	if (novo.eq.1) goto 2000
c	if (novo.eq.2) goto 1000

      end if

	
	
	
	
	
	
	
	
      end




