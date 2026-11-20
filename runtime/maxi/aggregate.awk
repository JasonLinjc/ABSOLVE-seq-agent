BEGIN { FS=OFS="\t" }
NR==FNR { seq[$1]=$2; next }
FNR%4==1 { if(substr($0,1,1)!="@") exit 2 }
FNR%4==2 {
    if(length($0)<48) exit 3
    b=substr($0,1,16); t=substr($0,40,9); a=substr($0,49,32)
    print b,t,((t in seq) && a==seq[t]) ? 1 : 0
    n++
}
END {if (!n || FNR!=4*n) exit 4}
