#!/bin/bash
# xml file has been copied to sro@stella:/stella/home/www/uploads/weingrill/save/
echo "fetch files from stella"
rsync -e ssh -r -v  -t sro@stella:/stella/home/www/uploads ~

# create xml files
cd /z/operator/uploads

#generate the xml file
#java -Djava.ext.dirs=/usr/share/java/ -cp /z/operator/java stella.jview.JTargetMaker\$Recreate proposal.create weingrill/save/$1

echo "move the xml file to /submit/"
if [ -e "*.xml" ]
then
        mv *.xml weingrill/submit/
else
        echo "no XML file found"
fi

echo "copy to Tenerife"
scp -p ~/uploads/weingrill/save/$1 stella@wifsip:stella/master1/testing/targets/
#java stella.util.SchedulerAccess Stella1MasterMind -add $1 STELLA1
ssh stella@wifsip java -Djava.ext.dirs=/usr/lib/j2se/ext stella.util.SchedulerAccess Stella1MasterMind -add stella/master1/testing/targets/$1 STELLA1
