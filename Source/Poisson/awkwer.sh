LOCATION="$1"
echo $LOCATION

while read line
do
    #echo "newline"
    #echo "$line"
    #echo "$line" | awk '{print "../../Data/Poissondata2/"$1}'
    echo "$line" | awk '{print "'$LOCATION'"}'
done < parameters.txt