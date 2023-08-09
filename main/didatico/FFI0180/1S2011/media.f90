
real :: rels, nota1, nota2, media

do
  read(*,*) rels, nota1, nota2
  media = (nota1 + nota2)/2.
  write(*,"( f5.2,tr2,f5.1 )") media, 0.04999999 + 0.4*rels + 0.6*media
end do

end
