#!/bin/csh

echo "file $1 trig $2 pt $3 frame $4"

nohup root -b << EOF
.x fitlambda.C($1,$2,$3,$4)
.q
EOF

