#
# Sources:
# https://stackoverflow.com/questions/50951339/gnuplot-plotting-series-based-on-label-in-third-column
#

set datafile separator ','
set key autotitle columnhead

set terminal png
set output 'errors.png'

delimiter = "#"    # some character that does not appear in category name
categories = ""

stats "errors.csv" using (categories = categories." ".delimiter.strcol(2).delimiter) nooutput

unique_categories = ""
do for [cat in categories] {
   if (strstrt (unique_categories, cat) == 0) {
      unique_categories = unique_categories." ".cat
   }
}


set logscale y
set logscale x

set xlabel "dt"
set ylabel "error"
set title "Euler-Maruyama and Milstein Error vs Step Size"

plot for [cat in unique_categories] "errors.csv" using 1:(delimiter.strcol(2).delimiter eq cat ? $3 : NaN) title cat[2:strlen(cat)-1]

