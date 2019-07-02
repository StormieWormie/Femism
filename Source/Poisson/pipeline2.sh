

#awk '{print "./DATA/"$1}' $1 | xargs mkdir
#awk '{print $2" "$3" "$4}' $1 | xargs $2
LOCATION="$3"
mkdir $LOCATION
while read line
do
    echo "newline"
    echo "$line"
    echo "$line" | awk '{print "'$LOCATION'/"$1}' | xargs mkdir
    echo "$line" | awk '{print "'$LOCATION'/"$1"/ "$2" "$3" "$4}' | xargs $1
done < $2