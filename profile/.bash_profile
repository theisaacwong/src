# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi

# User specific environment and startup programs

PATH=$PATH:$HOME/bin

export PATH

# RT #101835
umask u=rwx,g=rwx,o=rx


alias snakesub='mkdir -p log; snakemake --ri --jobname "{rulename}.{jobid}" --drmaa " -V -cwd -j y -o ./log -e ./log -l h_rt={resources.hrs}:00:00 -l mfree={resources.mem}G -pe serial {threads} -w n -S /bin/bash" -w 60'
alias stats='qstat -u $( whoami )'

if [ "$color_prompt" = yes ]; then
	PS1='\[\033[01;32m\]\u@\h:\[\033[01;34m\]\W\[\033[00m\]$ '
	#PS1='[\u@\h \W]\$ '
else
    PS1='[\u@\h \W]\$ '
fi

export DRMAA_LIBRARY_PATH=/opt/sge/drmaa/lib/lx-amd64/libdrmaa.so.1.0

. /etc/profile.d/modules.sh
if test ! -z $MODULESHOME; then
   module load modules modules-init modules-gs/prod modules-eichler/prod
fi

alias ls='ls --color=auto'

alias l='ls -lhrt'
alias lss='less -S'
alias awk9="awk '{print \$9}'"
alias t2n="tr '\t' '\n'"
alias cpgen="cp /net/eichler/vol28/home/iwong1/nobackups/scripts/qsub/template.qsub ."
alias R="R --no-save"
alias loadr="module load R/4.3.2"
alias collate="Rscript --vanilla ~/nobackups/scripts/R/collateTables.R "

header() {
	head -1 $1 | t2n | awk '{print NR " " $0}'
}

dim() {
	NCOL=$(head -1 $1 | grep -Po "\t" | wc -l)
	NROW=$(wc -l $1 | grep -Po [1-9]*(?= ))
	echo $NROW  $NCOL
}

powerReader() { column -t -s$'\t' "$@" | less -S; }
pr() { column -t -s$'\t' "$@" | less -S; }

alias cdautism="cd /net/eichler/vol28/projects/autism_genome_assembly/nobackups/"
alias cdlra="cd /net/eichler/vol28/projects/long_read_archive/nobackups/"
alias cdont="cd /net/eichler/vol28/projects/ont_sequencing/backups/"
alias cdhome="cd /net/eichler/vol28/home/iwong1/nobackups/"
alias cdprojects="cd /net/eichler/vol28/projects/"
alias qsta="qstat -u iwong1"

source /net/eichler/vol28/home/iwong1/nobackups/scripts/bash/acd_func.sh
PS1='\[\033[01;32m\]\u@\h:\[\033[01;34m\]\W\[\033[00m\]$ '

export MODSW=/net/eichler/vol28/software/modules-sw/
export MODREP=/net/eichler/vol28/7200/software/modules-repo/${MODULES_REL}

export APPY="/net/eichler/vol28/home/iwong1/nobackups/apptainer/"

LS_COLORS='rs=0:di=01;34:ln=01;36:mh=00:pi=40;33:so=01;35:do=01;35:bd=40;33;01:cd=40;33;01:or=40;31;01:mi=00:su=37;41:sg=30;43:ca=30;41:tw=30;42:ow=34;42:st=37;44:ex=01;32:*.tar=01;31:*.tgz=01;31:*.arc=01;31:*.arj=01;31:*.taz=01;31:*.lha=01;31:*.lz4=01;31:*.lzh=01;31:*.lzma=01;31:*.tlz=01;31:*.txz=01;31:*.tzo=01;31:*.t7z=01;31:*.zip=01;31:*.z=01;31:*.dz=01;31:*.gz=01;31:*.lrz=01;31:*.lz=01;31:*.lzo=01;31:*.xz=01;31:*.zst=01;31:*.tzst=01;31:*.bz2=01;31:*.bz=01;31:*.tbz=01;31:*.tbz2=01;31:*.tz=01;31:*.deb=01;31:*.rpm=01;31:*.jar=01;31:*.war=01;31:*.ear=01;31:*.sar=01;31:*.rar=01;31:*.alz=01;31:*.ace=01;31:*.zoo=01;31:*.cpio=01;31:*.7z=01;31:*.rz=01;31:*.cab=01;31:*.wim=01;31:*.swm=01;31:*.dwm=01;31:*.esd=01;31:*.jpg=01;35:*.jpeg=01;35:*.mjpg=01;35:*.mjpeg=01;35:*.gif=01;35:*.bmp=01;35:*.pbm=01;35:*.pgm=01;35:*.ppm=01;35:*.tga=01;35:*.xbm=01;35:*.xpm=01;35:*.tif=01;35:*.tiff=01;35:*.png=01;35:*.svg=01;35:*.svgz=01;35:*.mng=01;35:*.pcx=01;35:*.mov=01;35:*.mpg=01;35:*.mpeg=01;35:*.m2v=01;35:*.mkv=01;35:*.webm=01;35:*.webp=01;35:*.ogm=01;35:*.mp4=01;35:*.m4v=01;35:*.mp4v=01;35:*.vob=01;35:*.qt=01;35:*.nuv=01;35:*.wmv=01;35:*.asf=01;35:*.rm=01;35:*.rmvb=01;35:*.flc=01;35:*.avi=01;35:*.fli=01;35:*.flv=01;35:*.gl=01;35:*.dl=01;35:*.xcf=01;35:*.xwd=01;35:*.yuv=01;35:*.cgm=01;35:*.emf=01;35:*.ogv=01;35:*.ogx=01;35:*.aac=00;36:*.au=00;36:*.flac=00;36:*.m4a=00;36:*.mid=00;36:*.midi=00;36:*.mka=00;36:*.mp3=00;36:*.mpc=00;36:*.ogg=00;36:*.ra=00;36:*.wav=00;36:*.oga=00;36:*.opus=00;36:*.spx=00;36:*.xspf=00;36:'
export LS_COLORS

app() {
        apptainer exec -C --bind /net:/net ${1} /bin/bash
}

appy() {
	apptainer exec -C --bind /net:/net /net/eichler/vol28/home/iwong1/nobackups/apptainer/${1}.sif /bin/bash
}

logr() {
	ls -lh | awk '{print $9}' | xargs -I {} sh -c "echo '************ {} ************'; cat {}; echo '\n\n______________________________________________________________________________________________________________________________________________\n----------------------------------------------------------------------------------------------------------------------------------------------\n\n'" | lss
}

trimlogs() {
	ls -lhd */ | grep iwong1 | head -n -1 | awk '{print $9}' | xargs rm -rf
}

alias cpgen="cp /net/eichler/vol28/home/iwong1/nobackups/src/snake/runsnake ."

# Usage: diskusage vol28
alias diskusage='/net/eichler/vol28/7200/software/tools/storage/diskusage'

cp ~/.bash_profile ~/nobackups/src/profile/


