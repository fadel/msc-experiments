# stop.sh
#
# Generate proper stop words list from the 'stop.txt' file.


# Original source: http://snowball.tartarus.org/algorithms/english/stop.txt
# NOTE: in our experiments, stop.txt has been modified to include the last stop
# words (stop.txt is included).

sed 's/|.*//g' <stop.txt \
    | sed 's/ \+//g' \
    | sed '/^$/d' >words.txt
