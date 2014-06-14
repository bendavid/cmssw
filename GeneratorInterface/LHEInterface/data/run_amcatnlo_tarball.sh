#!/bin/bash

echo "   ______________________________________     "
echo "         Running Madgraph5                    "
echo "   ______________________________________     "

repo=${1}
echo "%MSG-MG5 repository = $repo"

name=${2} 
echo "%MSG-MG5 gridpack = $name"

nevt=${3}
echo "%MSG-MG5 number of events requested = $nevt"

rnum=${4}
echo "%MSG-MG5 random seed used for the run = $rnum"

LHEWORKDIR=`pwd`

if [[ -d madevent ]]
    then
    echo 'madevent directory found'
    echo 'Setting up the environment'
    rm -rf madevent
fi
mkdir madevent; cd madevent

# retrieve the wanted gridpack from the official repository 
fn-fileget -c `cmsGetFnConnect frontier://smallfiles` ${repo}/${name}_tarball.tar.gz 

#check the structure of the tarball
tar xzf ${name}_tarball.tar.gz ; rm -f ${name}_tarball.tar.gz ;

cd mgbasedir/${name}

#make sure lhapdf points to local cmssw installation area
echo "lhapdf = `echo "$LHAPATH/../../../full/bin/lhapdf-config"`" >> ./Cards/amcatnlo_configuration.txt

echo "done" > runscript.dat
echo "set nevents $nevt" >> runscript.dat
echo "set iseed $rnum" >> runscript.dat
echo "done" >> runscript.dat

domadpsin=0
if [ -f ./Cards/madspin_card.dat ] ;then
  #set random seed for madspin
  rnum2=$(($rnum+1000000))
  echo "$(echo `echo "set seed $rnum2"` | cat - ./Cards/madspin_card.dat)" > ./Cards/madspin_card.dat
  domadpsin=1
fi

#generate events
cat runscript.dat | ./bin/generate_events -ox -n $name

if [ "$domadspin" -gt "0" ] ; then 
    mv ./Events/${name}_decayed_1/events.lhe.gz $LHEWORKDIR/${name}_final.lhe.gz
else
    mv ./Events/${name}/events.lhe.gz $LHEWORKDIR/${name}_final.lhe.gz
fi

cd $LHEWORKDIR
gzip -d ${name}_final.lhe.gz
#cp ${name}_final.lhe ${name}_final.lhe.bak

ls -l
echo

exit 0

