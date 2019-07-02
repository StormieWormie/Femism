

#awk '{print "./DATA/"$1}' $1 | xargs mkdir
#awk '{print $2" "$3" "$4}' $1 | xargs $2
LOCATION="$3"
mkdir $LOCATION
NUMBER=$(< $2 wc -l)
while read line
do


    echo "$line" | awk '{print $1"\tout of'"$NUMBER"'"}'
    echo "$line" | awk '{print "'$LOCATION'/"$1}' | xargs mkdir
    #echo "$line" | awk '{print "'$LOCATION'/"$1"/ "$2" "$3" "$4}' #| xargs $1
    echo "$line" | awk '{printf "'$LOCATION'/"$0"\n"}' | xargs $1
done < $2
