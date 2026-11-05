/hdd/yanhuan/envs/tum/bin/samtools view -F 0x900 "$1" \
| awk '
{
  cigar=$6;
  nm=0;
  aln=0;


  for(i=12;i<=NF;i++){
    if($i ~ /^NM:i:/){
      split($i,a,":");
      nm=a[3];
    }
  }

  tmp=cigar;
  while(match(tmp, /[0-9]+[M=X]/)){
    aln += substr(tmp, RSTART, RLENGTH-1);
    tmp = substr(tmp, RSTART+RLENGTH);
  }

  if(aln > 45 && (1 - nm/aln) > 0.95)
    print $1
}' \
| sort -u
