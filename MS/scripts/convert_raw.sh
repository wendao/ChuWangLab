
work_dir=$(pwd)

for raw in $*
do
 echo "cd /home/wendao/RawConverter/Release/; env WINEPREFIX=~/.wine32 wine RawConverter.exe ${work_dir}/${raw} --mgf --mzXML --select_mono_prec &> ${work_dir}/${raw}.log; cd $work_dir"
done | parallel -j 4

