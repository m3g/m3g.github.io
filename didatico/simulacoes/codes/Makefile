all : aleatorio ativ29 atividade18 min1 min2 randomsearch2 randomsearch sim1 sim2 simplex \
      kineticmodel \
      modelab
aleatorio : aleatorio.f90
	gfortran -o aleatorio     aleatorio.f90
ativ29 : ativ29.f90
	gfortran -o ativ29        ativ29.f90
atividade18 : atividade18.f90
	gfortran -o atividade18   atividade18.f90
min1 : min1.f90
	gfortran -o min1          min1.f90
min2 : min2.f90
	gfortran -o min2          min2.f90
randomsearch2 : randomsearch2.f90
	gfortran -o randomsearch2 randomsearch2.f90
randomsearch : randomsearch.f90
	gfortran -o randomsearch  randomsearch.f90
sim1 : sim1.f90
	gfortran -o sim1          sim1.f90
sim2 : sim2.f90
	gfortran -o sim2          sim2.f90
simplex : simplex.f90
	gfortran -o simplex       simplex.f90
kineticmodel : kineticmodel.f90
	gfortran -o kineticmodel kineticmodel.f90
modelab : modelab.f90 NelderMeadMinimizer.f90
	gfortran -o modelab modelab.f90 NelderMeadMinimizer.f90
