#!/bin/sh

currentDir=$PWD

# First, get a list of subfolders containing the .wav files

rm dirlist

cd $1
ls -d * > ../dirlist
cd ../

# Define path to some files (sonobat, config, submit)

fileConfigTemp=$currentDir/SonobatSever20170125.conf
sonobatFile=$currentDir/SonobatSever20170125
submitFile=$currentDir/submit.sh

configFileName=SonobatSever20170125.conf


# use loop to get the subfolder name from dirlist
while IFS='' read -r line || [[ -n "$line" ]]; do
	echo $line
	fullDir=$currentDir/$1/$line

	# define where the job will be running
	workDir=$currentDir/work_$1/$line
	mkdir -p $workDir


	# define entries of TargetDirectory and LogDirectory 
	targetDirEntry='TargetDirectory="'$fullDir'"'
	logDirEntry='LogDirectory="'$workDir'"'


	# use bash commands to generate the config file for each subfolder and store and temp.config
	sed -n 1,14p $fileConfigTemp > temp.config
	echo $logDirEntry >> temp.config
	sed -n 16,17p $fileConfigTemp >> temp.config
	echo $targetDirEntry >> temp.config
	sed -n 19,28p $fileConfigTemp >> temp.config
	

	# copy/move files needed for computation/submitting to workDir
	mv temp.config $workDir/$configFileName
	cp $sonobatFile $workDir/
	cp $submitFile $workDir/

	
	# submit jobs to compute node
	cd $workDir
	msub -N $1/${line} $submitFile


	# move back 
	cd $currentDir
	
done < ./dirlist

