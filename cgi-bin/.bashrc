# .bashrc

# User specific aliases and functions

# Source global definitions

if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

#Perguntar antes de apagar, mover, etc:
alias rm='rm -i'
alias cp='cp -i'
alias mv='mv -i'
alias mkdir='mkdir -p'
alias ls='ls -hF --color=never'
alias update='echo "Update system..."; sudo apt-get update; sudo apt-get dist-upgrade'

#Fix wacom
alias wacom=`xsetwacom set "Wacom Bamboo Pen Pen stylus" Mode Relative`


#Move to trash instead of deleting
alias rm='trash-put'

# Don't ask if I want to add keys all the time and export X
alias ssh='ssh -X -o StrictHostKeyChecking=no'

# Prompt de comando
#PS1="\e[0;33m\u@\h:\w%\e[m " 
PS1="\u@\h:\w% "

# PATH
PATH=/sbin:/usr/sbin:/usr/local/bin:/bin:/usr/bin
PATH=$PATH:/usr/X11R6/bin
PATH=$PATH:/home/leandro/programs/lovoalign
PATH=$PATH:/home/leandro/programs/mdlovofit/bin
PATH=$PATH:/home/leandro/programs/packmol
PATH=$PATH:/home/leandro/programs/NAMD_2.9_Linux-x86_64
PATH=$PATH:/home/leandro/programs/molden
PATH=$PATH:/home/leandro/programs/vmd
PATH=$PATH:/usr/local/grace/bin
PATH=$PATH:/home/leandro/programs/mdanalysis/bin
PATH=$PATH:/home/leandro/programs/g03
PATH=$PATH:/home/leandro/programs/gv4
PATH=$PATH:/home/leandro/programs/voidoo
PATH=$PATH:/home/leandro/programs/atd_scripts
PATH=$PATH:/home/leandro/Documents/hf/gamess/gamess

# Path dos SCRIPTS
PATH=$PATH:/home/leandro/programs/scripts

# Path to toyprob ampl
PATH=$PATH:/home/leandro/Documents/termobio/algencan/bin/ampl

# Path to sasref
PATH=$PATH:/opt/atsas/bin

# Path to coot
PATH=$PATH:/home/leandro/programs/coot/bin

# Atualizar o PATH
export PATH

# Fix backspace
if [ `hostname` == "taz" ]; then
    stty erase '^?' >& /dev/null
   #stty erase ^H
   alias fixbackspace=`stty erase ^? >& /dev/null`
fi

# For namd on cuda
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/leandro/programs/NAMD_2.9_Linux-x86_64-multicore-CUDA

# For amber ptraj
PATH=$PATH:/home/leandro/programs/amber/amber11/bin 
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/leandro/programs/amber/amber11/bin
export LD_LIBRARY_PATH

