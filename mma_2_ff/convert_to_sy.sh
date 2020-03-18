#!/bin/bash
#==================================================================================
#    FireFly - Reconstructing rational functions and polynomial over finite fields.
#    Copyright (C) 2020  Jonas Klappert and Fabian Lange
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#==================================================================================

# Script to convert a Mathematica expression list to parsable format
file=$1
suff="_c"
file_c="$file$suff"

echo "Converting $file to $file_c."

cp $file $file_c

#tr -d '' < $file >  $file_c

sed -i 's/"//g' $file_c
sed -i 's/{//g' $file_c
sed -i 's/}/;/g' $file_c
sed -i 's/\\//g' $file_c
sed -i 's/,/;\n/g' $file_c
