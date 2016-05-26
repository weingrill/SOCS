#!/bin/bash
# store first command line parameter
FILE=$1 #m_67_rot_nw_20141202.save
# get filename without extension
filename=${FILE%.*}
#target contains filename without extension and _date
target=${filename%?????????}

echo "fetch files from stella"
rsync -e ssh -r -v  -t sro@stella:/stella/home/www/uploads ~

# create xml files
cd /z/operator/uploads

echo "generate the xml file"
java -Djava.ext.dirs=/usr/share/java/ -cp /z/operator/java stella.jview.JTargetMaker\$Recreate proposal.create weingrill/save/$filename.save

echo "move the xml file to /submit/"
if [ -e "*.xml" ]
then
	mv *.xml weingrill/submit/
else
	echo "no XML file found"
fi

echo "copy to Tenerife"
scp -p ~/uploads/weingrill/submit/$target.xml stella@wifsip:stella/master1/testing/targets/
#java stella.util.SchedulerAccess Stella1MasterMind -add $target.xml STELLA1

echo "activate target in observing queue"
ssh stella@wifsip java -Djava.ext.dirs=/usr/lib/j2se/ext stella.util.SchedulerAccess Stella1MasterMind -add stella/master1/testing/targets/$target.xml STELLA1
