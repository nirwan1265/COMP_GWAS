for famfile in *.fam
do
famname=${famfile%.*}
txtfile="${famname}.txt"
if [[ -f "$txtfile" ]]; then

  awk 'FNR==NR{a[NR]=$1;next}{$6=a[FNR];print}' "$txtfile" "$famfile" > "${famname}_new.fam"
  mv "${famname}_new.fam" "$famfile"
fi
done
