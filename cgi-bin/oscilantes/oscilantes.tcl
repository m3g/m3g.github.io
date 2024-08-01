#!/usr/bin/tclsh

set basehtml "http://leandro.iqm.unicamp.br/leandro/html/didatico/qf632/oscilantes"
set filehtml "$basehtml/files"
set files "/home/leandro/public_html/leandro/html/didatico/qf632/oscilantes/files"

puts "Content-type: text/html"
puts ""
puts "<head>"
puts "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\">"
puts "<link rel=\"stylesheet\" type=\"text/css\" href=$basehtml/estilo.css>"
puts "<title>Reações Oscilantes</title>"
puts "</head>"
puts "<body>"
puts "<table width=80% align=center><tr><td align=center>"
puts "<br><h2>Resultado da simulação:</h2><br>"

catch {
set data [read stdin $env(CONTENT_LENGTH)]
set data [ split $data "\n" ]
set i 0
set ivar 0
set ival 0
foreach arg $data {
  incr i
  set val($i) $arg
  if { [ string first "name=" $val($i) ] > 0 } {
    set ivar $i
  }
  if { $i == [ expr $ivar + 2 ] } {
    incr ival
    set value($ival) $val($i)
  }

}

if { $value(11) > 100. } { 
  puts " ERRO: O tempo deve ser no máximo 100.<br><br> "
  puts "<a href=$basehtml>\[voltar\]</a>"
  exit
}
for { set i 1 } { $i <= 11 } { incr i } {
  if { $value($i) < 0. } {
    puts " ERRO: Nenhum dos parâmetros pode ser negativo. <br><br>"
    puts "<a href=$basehtml>\[voltar\]</a>"
    exit
  }
} 

set rand [ expr int(1000*rand()) ]

catch { 
exec -ignorestderr ./oscilantes $value(1) $value(2) $value(3) $value(4) $value(5) $value(6) $value(7) $value(8) $value(9) $value(10) $value(11) $files/dados_$rand.dat
} msg
puts $msg

set GREP [ catch { exec grep ERRO $files/dados_$rand.dat } ERROR ]
if { $GREP == 0 } { 
  puts "$ERROR <br>"
  puts "Tente com novos parâmetros iniciais <br> (geralmente concentrações ou constantes de velocidade menores)<br>"
  puts "<a href=$basehtml>\[voltar\]</a>"
  exit
} 

exec -ignorestderr sed -e s/DADOSRAND/dados_$rand/ plot.gnp > $files/plot_$rand.gnp 

exec -ignorestderr rm -rf $files/dados_$rand.jpg
exec -ignorestderr gnuplot $files/plot_$rand.gnp 2>@ stderr
puts "Arquivo de dados: <a href=$filehtml/dados_$rand.dat>dados_$rand.dat</a><br>"
puts "<img src=$filehtml/dados_$rand.jpg>" 

puts "<br><br>"
puts "<b>Dados:</b> <br>"
puts "\[A\]<sub>0</sub> = $value(1) &nbsp;&nbsp;&nbsp;"
puts "\[B\]<sub>0</sub> = $value(2) &nbsp;&nbsp;&nbsp;"
puts "\[C\]<sub>0</sub> = $value(3) &nbsp;&nbsp;&nbsp;"
puts "\[D\]<sub>0</sub> = $value(4)<br>"
puts "k<sub>1</sub> = $value(5) &nbsp;&nbsp;&nbsp;"
puts "k<sub>-1</sub> = $value(6) &nbsp;&nbsp;&nbsp;"
puts "k<sub>2</sub> = $value(7) &nbsp;&nbsp;&nbsp;"
puts "k<sub>-2</sub> = $value(8) &nbsp;&nbsp;&nbsp;"
puts "k<sub>3</sub> = $value(9) &nbsp;&nbsp;&nbsp;"
puts "k<sub>-3</sub> = $value(10)<br>"
puts "Tempo = $value(11)<br><br>"
puts "<a href=http://leandro.iqm.unicamp.br/leandro/html/didatico/qf632/oscilantes>\[voltar\]</a>"

} msg
puts $msg
puts "</td></tr></table>"
puts "</body>"
