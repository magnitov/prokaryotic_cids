FOLDER_PATH=$1
REPLICA=$2
RESOLUTION=$3

THR=$((RESOLUTION*3))

for IS in $((RESOLUTION*5)) $((RESOLUTION*6)) $((RESOLUTION*7)) $((RESOLUTION*8)) $((RESOLUTION*9)) $((RESOLUTION*10)) $((RESOLUTION*11)) $((RESOLUTION*12)) $((RESOLUTION*13)) $((RESOLUTION*14)) $((RESOLUTION*15))
do
	for IDS in $((RESOLUTION*2)) $((RESOLUTION*4)) $((RESOLUTION*6)) $((RESOLUTION*8)) $((RESOLUTION*10))
	do
		IS_2=$((IS+1))
		IDS_2=$((IDS+1))
		perl scripts/matrix2insulation.pl -i ${FOLDER_PATH}/${REPLICA}.tab -is ${IS} -ids ${IDS} -im mean -nt 0.1 -bmoe 2 -v
	done
done
