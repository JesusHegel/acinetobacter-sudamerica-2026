BEGIN{FS=OFS="\t"; V=500}
FNR==NR && FILENAME==ARGV[1]{len[$1"|"$2]=$3; next}
FILENAME==ARGV[2]{is[$1"|"$2]=is[$1"|"$2]" "$3","$4; next}
{
  acc=$1; ctg=$2; ini=$3; fin=$4; hebra=$5; gen=$6
  k=acc"|"ctg; L=len[k]
  estado="ausente"
  if(hebra=="+"){ margen=ini-1; zi=ini-V; zf=ini }
  else          { margen=L-fin;  zi=fin;   zf=fin+V }
  n=split(is[k],a," ")
  for(i=1;i<=n;i++){
    if(a[i]=="") continue
    split(a[i],b,","); s=b[1]; e=b[2]
    if(s>e){t=s; s=e; e=t}
    if(e>=zi && s<=zf) estado="presente"
  }
  if(estado=="ausente" && margen<V) estado="no_determinable"
  print acc, gen, ctg, hebra, L, margen, estado
}
