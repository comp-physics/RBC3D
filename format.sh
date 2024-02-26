echo "formatting example cases and common directory"

echo "$PWD"
fprettify ./examples -r --indent 2
fprettify ./common -r --indent 2

echo "done"
