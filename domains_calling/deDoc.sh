#!/bin/bash
FOLDER_PATH=$1

FILES=${FOLDER_PATH}/formatted_maps/graph/*
for f in $FILES
do
	echo $f
	time /usr/bin/java -jar /home/magnitov/Software/deDoc-1.0.0/deDoc.jar $f
done

rm ${FOLDER_PATH}/formatted_maps/graph/*uncontinued*
mkdir -p ${FOLDER_PATH}/TADs/deDoc/
mv ${FOLDER_PATH}/formatted_maps/graph/*deDoc* ${FOLDER_PATH}/TADs/deDoc/

